

#include	"CGBL_dmesh.h"
#include	"CreateMesh.h"

#include	"ImageIO.h"
#include	"Timer.h"

#include	<queue>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// compass directions
static int dxs[8] = {1, 1, 0, -1, -1, -1,  0,  1};
static int dys[8] = {0, 1, 1,  1,  0, -1, -1, -1};






/* --------------------------------------------------------------- */
/* class PQElm --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Priority queue element - track costs for cutting corners

class PQElm {
public:
	int		to;		// from node zero to here...
	double	cost;	// and the cost to get here
public:
	PQElm( int to, double cost ) : to(to), cost(cost) {};

	bool operator < ( const PQElm &rhs ) const
		{return cost > rhs.cost;};	// priority low if cost high
};

/* --------------------------------------------------------------- */
/* class Grf ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Graph node - nodes are perimeter corners

class Grf{
public:
	int		x, y,	// vertex coordinates
			back;	// node we came from
	double	cost;	// cost to get here from node zero
public:
	Grf( const vertex &corner, int back, double cost )
	: x(corner.x), y(corner.y), back(back), cost(cost) {};

	double Dist( const Grf& rhs ) const;
};


double Grf::Dist( const Grf& rhs ) const
{
	double	dx = x - rhs.x,
			dy = y - rhs.y;

	return sqrt( dx*dx + dy*dy );
}

/* --------------------------------------------------------------- */
/* SegPointDist -------------------------------------------------- */
/* --------------------------------------------------------------- */

static double SegPointDist(
	const Grf	&v0,
	const Grf	&v1,
	const Grf	&v2 )
{
	return SegPointDist( v0.x, v0.y, v1.x, v1.y, v2.x, v2.y );
}

/* --------------------------------------------------------------- */
/* ShortcutCost -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Assume a line from entry s to entry e of corner vertices.
// The cost of this shortcut is the sum of all perpendicular
// displacements of excluded intermediate vertices.
//
static double ShortcutCost( const vector<Grf> &v, int s, int e )
{
	double	sum = 0.0;

	for( int i = s + 1; i < e; ++i )
		sum += SegPointDist( v[s], v[e], v[i] );

	return sum;
}

/* --------------------------------------------------------------- */
/* SetEdgeLimits ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SetEdgeLimits(
	double		&lmin,
	double		&lmax,
	const IBox	&B,
	FILE*		flog )
{
	int		bx = B.R - B.L, by = B.T - B.B;

// Prefer triangle edges in range [MNL, 2xMNL] and generally
// prefer a max:min ratio of about 2:1.

	lmin = GBL.mch.MNL,
	lmax = 2 * lmin;

// If region skinnier than typical triangle...

	if( bx < lmin || by < lmin ) {

		// set trial values
		lmin = fmin( bx, by );
		lmax = fmax( bx, by );

		// push lmin down more if lmax small
		lmin = fmin( lmin, lmax / 2.0 );

		// split longer dimension
		lmax = lmax / 3.0;

		// push lmax up more to restore 2:1
		lmax = fmax( lmax, 2.0 * lmin );

		fprintf( flog,
		"Reducing minl from [%d %d] to [%f %f].\n",
		GBL.mch.MNL, 2*GBL.mch.MNL, lmin, lmax );
	}
}

/* --------------------------------------------------------------- */
/* MakeMap ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Makes a bitmap (ones where there are pts) that adds a one pixel
// wide border of zeros all around the points. This map will serve
// as the coordinate system for mesh building. MeshCreate must
// offset control points to real coords before exiting.
//
static void MakeMap(
	vector<uint8>		&map,
	int					&w,
	int					&h,
	const IBox			&B,
	const vector<Point>	&pts )
{
	w = B.R - B.L + 3,
	h = B.T - B.B + 3;

	map.resize( w * h, 0 );

	int	npts = pts.size();

	for( int i = 0; i < npts; ++i ) {

		int	x = int(pts[i].x) - B.L + 1;
		int y = int(pts[i].y) - B.B + 1;

		map[x + w*y] = 1;
	}
}

/* --------------------------------------------------------------- */
/* FindLeftEdge -------------------------------------------------- */
/* --------------------------------------------------------------- */

// If we can find the left edge of the bounds set (x0,y0)
// and return true.
//
static bool FindLeftEdge(
	int					&x0,
	int					&y0,
	const vector<uint8>	&map,
	int					w,
	int					h,
	FILE*				flog )
{
// We position ourselves at horizizontal midline (1,h/2)
// and advance to the right until finding a point.

	x0 = 1;
	y0 = h/2;

	for( ; x0 < w; ++x0 ) {

		if( map[x0 + w*y0] )
			return true;
	}

// Failed.
// This could happen if the overlap cut the connected region
// into two parts, neither of which intersects the horizonal
// midline. Return failure.

	fprintf( flog, "Can't find starting point??\n" );

#ifdef CREATEMESH_WRITE_DEBUG_IMAGES
	Raster8ToTif8( "Bogus3.tif", &map[0], w, h );
#endif

	return false;
}

/* --------------------------------------------------------------- */
/* GetBoundaryPoints --------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if successfully find a contiguous boundary.
//
static bool GetBoundaryPoints(
	vector<vertex>		&bndry,
	const vector<uint8>	&map,
	int					w,
	int					h,
	int					x0,
	int					y0,
	FILE*				flog )
{
// Starting at x0, y0 and working CCW, collect all points on the
// region boundary.

	int		x = x0, y = y0, dir = 0;

	do {

		//fprintf( flog, "-- %6d %6d %3d\n", x, y, dir );

		int	op = (dir+4)%8;	// sweep starts opposite to dir
		int	k;

		// Sweep dir CCW until pointing at non-zero neighbor.
		// Note that the sweep always starts behind us and is
		// always on our exterior flank.

		for( k = 1; k <= 8; ++k ) {

			int	j	= (op+k)%8,
				dx	= dxs[j],
				dy	= dys[j];

			if( map[x+dx + w*(y+dy)] ) {

				bndry.push_back( vertex( x, y, j ) );

				x	= x + dx;
				y	= y + dy;
				dir	= j;

				goto next_bndry_pt;
			}
		}

		// Failed

		fprintf( flog, "Input was an isolated pixel??\n" );

#ifdef CREATEMESH_WRITE_DEBUG_IMAGES
		Raster8ToTif8( "Bogus4.tif", &map[0], w, h );
#endif

		return false;

next_bndry_pt:;
	} while( x != x0 || y != y0 );

	return true;
}

/* --------------------------------------------------------------- */
/* FindCorners --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void FindCorners(
	vector<vertex>			&corners,
	const vector<vertex>	&bndry,
	FILE*					flog )
{
// Compress the ordered boundary points into a list of corners,
// that is, the set of points at which the direction changes.

	int		nbps = bndry.size();

	for( int i = 0; i < nbps; ++i ) {

		int	iprev = (i+nbps-1)%nbps;

		if( bndry[i].dir != bndry[iprev].dir ) {

			vertex	v = bndry[i];

			v.orig = i;
			corners.push_back( v );
		}
	}

	int	ncorn = corners.size();

// print corner count

	fprintf( flog,
	"Got %d corners from %d boundary points.---\n", ncorn, nbps );

// table of corners

	fprintf( flog, " Vtx\t Point\t    Vx\t    Vy\n" );

	for( int i = 0; i < ncorn; ++i ) {

		const vertex&	C = corners[i];

		fprintf( flog, "%4d\t%6d\t%6d\t%6d\n", i, C.orig, C.x, C.y );
	}

	fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* SplitLongSides ------------------------------------------------ */
/* --------------------------------------------------------------- */

// If any edges are too long, split them.
//
// Note that a bit of oversplitting actually helps
// the corner cutting operation that follows this.
//
static void SplitLongSides(
	vector<vertex>	&corners,
	double			lmin,
	FILE*			flog )
{
	int	ncorn = corners.size();

	for( int i = 0; i < ncorn; ++i ) {

		int		inext	= (i+1)%ncorn;
		double	dx		= corners[inext].x - corners[i].x;
		double	dy		= corners[inext].y - corners[i].y;
		double	len		= sqrt( dx*dx + dy*dy );

		// len > lmin/2 ?

		if( 2.0 * len > lmin ) {

			int	nins = int(2.0*len/lmin) + 1;

			for( int j = 1; j <= nins; ++j ) {

				double	frac = double(j)/(nins+1);
				vertex	nv(	int(corners[i].x + dx*frac),
							int(corners[i].y + dy*frac) );

				corners.insert( corners.begin() + i + j, nv );

				fprintf( flog, "Inserted at (%d %d).\n", nv.x, nv.y );
			}

			i		+= nins;
			ncorn	+= nins;
		}
	}

// print corner count

	fprintf( flog,
	"\nSplitting long sides gives %d corners.---\n", ncorn );

// table of corners

	fprintf( flog, " Vtx\t    Vx\t    Vy\n" );

	for( int i = 0; i < ncorn; ++i ) {

		const vertex&	C = corners[i];

		fprintf( flog, "%4d\t%6d\t%6d\n", i, C.x, C.y );
	}

	fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* CutCorners ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Here we shrink the contour a bit and reduce the segment count
// by excluding a few outlier points (cutting corners). All line
// segments have an intrinsic 'cost' value, which is their length.
// The cost for cutting corners is the sum of the perpendicular
// distances of the cut corners from the new shortcut segment.
//
static void CutCorners(
	vector<lineseg>			&edges,
	double					&lmin,
	double					lmax,
	const vector<vertex>	&corners,
	FILE*					flog )
{
	int		ncorn		= corners.size();
	bool	big_print	= false;

// Attempt up to ten times until lmin adjusted small enough.

	for( int II=1; II<=10 && edges.size() < 3; ++II, lmin *= 0.8 ) {

		fprintf( flog,
		"Pass %d: Length range [%f, %f].\n\n", II, lmin, lmax );

		edges.clear();

		/* ---------- */
		/* Init graph */
		/* ---------- */

		vector<Grf>	graph;

		graph.reserve( ncorn + 1 );

		// Set all back pointers to -1, and all costs to infinity

		for( int i = 0; i < ncorn; ++i )
			graph[i] = Grf( corners[i], -1, double(BIG) );

		// Final entry is duplicate of 1st--the home position

		graph[ncorn] = graph[0];

		/* --------------------------------------------- */
		/* Build low cost path from first to last corner */
		/* --------------------------------------------- */

		// Priority queues, which work like heaps, make sure
		// the topmost element has the highest 'priority'. In
		// our case, the PQElm::comparison function assigns
		// highest priority to the lowest 'cost' element.
		//
		// We will build a path by enqueing low cost hops and
		// dequeing them again to see how we might make further
		// hops from there, until we get home. The priority
		// queue device makes sure we consider lower cost hops
		// ahead of others, which is a fine guess to make, but
		// this does not guarantee an optimal overall path. We
		// do not require optimality, though, just a good path.

		{
			/* ---------------------------------------- */
			/* Initialize queue of intermediate corners */
			/* ---------------------------------------- */

			priority_queue<PQElm>	q;

			q.push( PQElm( 0, 0.0 ) );
			graph[0].cost = 0.0;

			/* ------------------------------- */
			/* Process corners until back home */
			/* ------------------------------- */

			while( !q.empty() ) {

				/* ------------------------ */
				/* Dequeue an endpoint T.to */
				/* ------------------------ */

				PQElm	T = q.top();

				q.pop();

				fprintf( flog,
				"Node 0 to %4d; xy=(%4d %4d); cost=%f\n",
				T.to, graph[T.to].x, graph[T.to].y, T.cost );

				/* --------- */
				/* Home yet? */
				/* --------- */

				if( T.to == ncorn )
					break;

				/* ----------------------------------- */
				/* Is it a lower cost way to get here? */
				/* ----------------------------------- */

				if( T.cost > graph[T.to].cost )
					continue;

				/* ------- */
				/* Enqueue */
				/* ------- */

				// Starting from the dequeued node [s], consider
				// paths to all other downstream nodes [j]. If the
				// length is reasonable AND if the cost to get to
				// [j] is lower than the current cost at [j], then
				// update node [j] AND enter [j] in the queue to
				// consider how to proceed from there.

				int	s = T.to;  // the source node

				for( int j = s + 1; j < ncorn + 1; ++j ) {

					double	len = graph[j].Dist( graph[s] );

					if( big_print ) {
						fprintf( flog,
						"Target %3d; xy=(%4d %4d); len=%f\n",
						j, graph[j].x, graph[j].y, len );
					}

					/* ------------------ */
					/* Length reasonable? */
					/* ------------------ */

					if( lmin <= len && len <= lmax ) {

						double	cost =
						graph[s].cost + ShortcutCost( graph, s, j );

						if( big_print ) {
							fprintf( flog,
							"Cost: src=%f; tot=%f; trg=%f\n",
							graph[s].cost, cost, graph[j].cost );
						}

						/* ----------- */
						/* Lower cost? */
						/* ----------- */

						if( cost < graph[j].cost ) {

							graph[j].cost = cost;
							graph[j].back = s;

							if( big_print ) {
								fprintf( flog,
								"Pushing node %3d;"
								" cost=%f; back=%3d\n",
								j, cost, s );
							}

							/* --------------------- */
							/* Start again from here */
							/* --------------------- */

							q.push( PQElm( j, cost ) );
						}
					}
				}
			}
		}

		/* ------------------- */
		/* All the way around? */
		/* ------------------- */

		if( graph[ncorn].cost == BIG ) {

			fprintf( flog,
			"No path to final vertex??"
			" Will reduce lmin and try again.\n" );
			continue;
		}

		fprintf( flog, "\n" );

		/* ------------------- */
		/* Collect final edges */
		/* ------------------- */

		// Deduce connectivity by backtracking from final vertex
		// and creating the edge list in reverse order (inserting
		// new entries at front). This keeps polygon in same order
		// as the vertices from which it was made.

		for( int i = ncorn; i > 0; i = graph[i].back ) {

			const Grf&	prev	= graph[graph[i].back];
			double		d		= graph[i].Dist( prev );

			fprintf( flog,
			"Edge from (%4d %4d) to (%4d %4d); len %7.2f\n",
			prev.x, prev.y, graph[i].x, graph[i].y, d );

			edges.insert( edges.begin(),
			lineseg( prev.x, prev.y, graph[i].x, graph[i].y ) );
		}
	}

	fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* UncrossEdges -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Here we fix tangles of a specific type wherein we are looking
// for crossings among an ordered list [0..ne) of edges. Assume
// that the first crossing pair is [i], [k]; with k > i. Further
// assume that only these two are crossed, while any intervening
// segments [i+1..k-1] are locally ordered without crossings. In
// other words, there is a locally ordered sequence of vertices
// between [i] and [k] suffering a global order reversal.
//
// The vertices needing reversal include the tip of [i] through
// the tail of [k], since again, we assume proper ordering of all
// vertices outside the pathological section. We will reverse the
// subsection by swapping pairs of vertices symmetrically placed
// about the subsection midpoint.
//
static void UncrossEdges( vector<lineseg> &edges, FILE* flog )
{
	int	ne = edges.size(), err;

	do {

		/* ------------------- */
		/* Print current edges */
		/* ------------------- */

		for( int i = 0; i < ne; ++i ) {

			const lineseg&	L = edges[i];

			fprintf( flog,
			"Pgon edge %3d: (%4d %4d) (%4d %4d).\n",
			i, L.v[0].x, L.v[0].y, L.v[1].x, L.v[1].y );
		}

		/* ---------------------------- */
		/* Look at all edge pairs {i,k} */
		/* ---------------------------- */

		err = false;

		for( int i = 0; i < ne - 1 && !err; ++i ) {

			for( int k = i + 1; k < ne && !err; ++k ) {

				const lineseg&	Li = edges[i];
				const lineseg&	Lk = edges[k];

				if( OpenSegsCross(
						Li.v[0], Li.v[1],
						Lk.v[0], Lk.v[1] ) ) {

					fprintf( flog,
					"Edges %d and %d cross; (%4d %4d)-(%4d %4d)"
					" and (%4d %4d)-(%4d %4d).\n", i, k,
					Li.v[0].x, Li.v[0].y, Li.v[1].x, Li.v[1].y,
					Lk.v[0].x, Lk.v[0].y, Lk.v[1].x, Lk.v[1].y );

					err = true;

					/* ------ */
					/* Repair */
					/* ------ */

					// It's cumbersome to address vertices as
					// tails and tips of edges, so temporarily
					// copy them to a simple vertex list. Note
					// that getting tails indeed gets them all
					// since the edges describe a closed loop.
					// Every tip is someone else's tail.

					vector<vertex>	v( ne );

					for( int ie = 0; ie < ne; ++ie )
						v[ie] = edges[ie].v[0];

					// Reverse inclusive range [i+1..k]
					// using symmetric pairwise swaps.

					int	nv = k - (i+1) + 1,	// num verts
						md = nv / 2;		// midpoint

					for( int j = 0; j < md; ++j ) {

						int		a	= (i+1+j)%ne,
								b	= (k-j)%ne;
						vertex	t;

						t		= v[a];
						v[a]	= v[b];
						v[b]	= t;
					}

					// Copy all vertices back to edges

					for( int ie = 0; ie < ne; ++ie ) {

						edges[ie].v[0] = v[ie];
						edges[ie].v[1] = v[(ie+1)%ne];
					}
				}
			}
		}

	} while( err );
}

/* --------------------------------------------------------------- */
/* ReverseWhole -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Here we address the worry that multiple rounds of section
// reversal in UncrossEdges() may have left the whole figure
// reversed. We can detect that using a negative area test.
//
static void ReverseWhole( vector<lineseg> &edges, FILE* flog )
{
	double	area = AreaOfPolygon( edges );

	fprintf( flog, "After edge untangling area %.2f\n", area );

	if( area < 0.0 ) {

		// reverse every edge

		fprintf( flog, "Reversing...\n" );

		int	ne = edges.size();

		for( int i = 0; i < ne; ++i ) {

			lineseg&	L = edges[i];
			vertex		t;

			t		= L.v[0];
			L.v[0]	= L.v[1];
			L.v[1]	= t;
		}
	}
}

/* --------------------------------------------------------------- */
/* SetInternalVertices ------------------------------------------- */
/* --------------------------------------------------------------- */

// @@@ Note:
// The only criteria we apply for these points is that they are
// on or within the boundary and not too close to a corner or to
// each other. However, these may not make useful triangles if
// too close to an edge. Deeper in the interior may be better.
//
static void SetInternalVertices(
	vector<vertex>			&vinside,
	vector<uint8>			&map,
	int						w,
	int						h,
	const vector<lineseg>	&edges,
	int						rclear,
	FILE*					flog )
{
/* --------------------------------------------- */
/* Clear map in neighborhood of existing corners */
/* --------------------------------------------- */

	int	ne = edges.size(),
		np = w * h;

	for( int i = 0; i < ne; ++i )
		RemoveFromMap( map, w, h, edges[i].v[0], rclear );

/* ----------------------------------- */
/* Find well-separated internal points */
/* ----------------------------------- */

	for( int i = 0; i < np; ++i ) {

		if( map[i] ) {

			/* -------------------------------------------- */
			/* Clear neighborhood around prospective vertex */
			/* -------------------------------------------- */

			int	y = i / w;
			int	x = i - w * y;

			fprintf( flog,
			"Add internal vertex at (%4d %4d).\n", x, y );

			vertex	newv( x, y );

			RemoveFromMap( map, w, h, newv, rclear );

			/* -------------------------- */
			/* Apply simple 'inside' test */
			/* -------------------------- */

			// Make guaranteed outside point and count edge
			// crossings from there to proposed inside point.

			vertex	outside( -10, y );

			int	m = CountCrossings( edges, newv, outside );

			if( m & 1 )
				vinside.push_back( newv );
			else if( m ) {

				fprintf( flog,
				"Warning: Supposedly internal point is outside"
				" - %d crossings. Continuing.\n", m );
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* ListEdgesMatlab ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void ListEdgesMatlab(
	const vector<lineseg>	&edges,
	FILE*					flog )
{
	int	ne = edges.size();

	fprintf( flog,
	"\n--------- Start Matlab format ------------------\n" );

	for( int i = 0; i < ne; ++i ) {

		const lineseg&	L = edges[i];

		fprintf( flog,
		"x=[%d %d]; y=[%d %d]; plot(x,y); hold on;\n",
		L.v[0].x, L.v[1].x, L.v[0].y, L.v[1].y );
	}
}

/* --------------------------------------------------------------- */
/* LongestEdge --------------------------------------------------- */
/* --------------------------------------------------------------- */

static int LongestEdge( const vector<lineseg> &edges )
{
	double	Dmax	= -1.0;
	int		ne		= edges.size(),
			Imax	= -1;

	for( int i = 0; i < ne; ++i ) {

		double	d = edges[i].LenSqr();

		if( d > Dmax ) {

			Dmax = d;
			Imax = i;
		}
	}

	return Imax;
}

/* --------------------------------------------------------------- */
/* UniqueVerts --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Combined list of unique edge and internal vertices

class UVert {
public:
	vertex	v;
	int		type,
			indx;
public:
	UVert()	{};

	UVert( const vertex& v, int type, int indx )
	: v(v), type(type), indx(indx) {};
};


// Form a list of unique vertices with their types and source
// indices for use by the BestVertex() search.
//
// Return uv.size();
//
static int UniqueVerts(
	vector<UVert>			&uv,
	const vector<lineseg>	&edges,
	const vector<vertex>	&vinside )
{
	int	n = edges.size();

// First add the edge.v0, which are guaranteed unique

	for( int i = 0; i < n; ++i )
		uv.push_back( UVert( edges[i].v[0], 0, i ) );

// Next add unique edge.v1
//
// The v1's are unique from other v1's, but are they
// different from listed v0's?

	for( int i = 0; i < n; ++i ) {

		const vertex&	v = edges[i].v[1];

		// unique?

		for( int j = 0; j < n; ++j ) {

			if( v == edges[j].v[0] )
				goto next_i;
		}

		uv.push_back( UVert( v, 1, i ) );

next_i:;
	}

// Finally, add internal verts which are always unique

	n = vinside.size();

	for( int i = 0; i < n; ++i )
		uv.push_back( UVert( vinside[i], -1, i ) );

// Kill excess space and return count

	uv.resize( n = uv.size() );

	return n;
}

/* --------------------------------------------------------------- */
/* BestVertex ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set the type and list-index for the best vertex to complete
// triangle currently having side va -> vb.
//
// type = {-2=unset, -1=internal, 0=edge-v0, 1=edge-v1}.
// indx = index into vinside[] or edges[] according to type.
//
static void BestVertex(
	int						&type,
	int						&indx,
	const vertex			&va,
	const vertex			&vb,
	const vector<lineseg>	&edges,
	const vector<vertex>	&vinside,
	FILE*					flog )
{
	type	= -2;
	indx	= -1;

	vertex			vm( (va.x+vb.x)/2, (va.y+vb.y)/2 ); // midpoint
	vector<UVert>	uv;
	double			Dbest = BIG;
	int				nu;

	nu = UniqueVerts( uv, edges, vinside );

	for( int i = 0; i < nu; ++i ) {

		const vertex&	vc = uv[i].v;

		if( vc == va || vc == vb )
			continue;

		double	D = vm.DistSqr( vc ),
				A = AreaOfTriangle( va, vb, vc );
		int		L = LeftSide( va, vb, vc );

		fprintf( flog,
		"#%3d (%4d %4d); dist=%12.2f; left=%d; area=%11.2f; ",
		i, vc.x, vc.y, D, L, A );

		// require vc on interior side of va->vb
		if( !L ) {
			printf( "rjct: not L\n" );
			continue;
		}

		// require small...
		if( D >= Dbest ) {
			printf( "rjct: big D\n" );
			continue;
		}

		// ...but not too small
		if( A <= GBL.mch.MTA ) {
			printf( "rjct: sml A\n" );
			continue;
		}

		// don't cross any remaining edges
		if( AnyCrossing( edges, va, vc ) ) {
			printf( "rjct: crs va\n" );
			continue;
		}

		// ditto
		if( AnyCrossing( edges, vb, vc ) ) {
			printf( "rjct: crs vb\n" );
			continue;
		}

		// don't enclose other vertices
		if( AnyInside( va, vb, vc, edges, vinside ) ) {
			printf( "rjct: any inside\n" );
			continue;
		}

		Dbest	= D;
		type	= uv[i].type;
		indx	= uv[i].indx;

		printf( "keep: *\n" );
	}
}

/* --------------------------------------------------------------- */
/* AddTriangle --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void AddTriangle(
	vector<triangle>	&tri,
	vector<vertex>		&ctl,
	const vertex		&va,
	const vertex		&vb,
	const vertex		&vc )
{
	triangle		t;
	vector<vertex>	V( 3 );
	int				nc = ctl.size();

	V[0] = va;	// makes vertices indexable
	V[1] = vb;
	V[2] = vc;

// For each triangle vertex, see if matches existing control point
// and either refer to that or add the vertex as new control point.

	for( int iv = 0; iv < 3; ++iv ) {

		for( int ic = 0; ic < nc; ++ic ) {

			if( ctl[ic] == V[iv] ) {
				// matches--use it
				t.v[iv] = ic;
				goto next_iv;
			}
		}

		// no match--append
		ctl.push_back( V[iv] );
		t.v[iv] = nc++;

next_iv:;
	}

	tri.push_back( t );
}

/* --------------------------------------------------------------- */
/* UpdatePolygon ------------------------------------------------- */
/* --------------------------------------------------------------- */

// We've already deleted a full edge va->vb from the figure, but
// we now edit the edge list to keep the polygon continuous and
// closed. While new directed triangle edges run va->vb->vc->va,
// when we remove this triangle, the remaining polygon contour
// runs ...va->vc->vb... So these are the segments we consider
// for addition back into the figure.
//
// As we add new segments we must be careful to avoid overlaps
// with existing edges. IsSubseg() checks if a given edge is
// entirely contained in the 'test' segment, and if so, it clips
// 'edge' out and pushes any non-overlapped parts of 'test' back
// to the prospective addition list (StackOfSegs). In this case,
// we remove the contained edge from the figure.
//
// Segments with no overlaps are added to the polygon.
//
static void UpdatePolygon(
	vector<lineseg>	&edges,
	const vertex	&va,
	const vertex	&vb,
	const vertex	&vc )
{
// Create a stack of prospective new segments to add.
// Initialize stack with va->vc and vc->vb.

	stack<lineseg>	StackOfSegs;

	StackOfSegs.push( lineseg( va, vc ) );
	StackOfSegs.push( lineseg( vc, vb ) );

// Examine segments until none left

	for( ; StackOfSegs.size() > 0; ) {

		// Fetch a prospective segment

		lineseg	test = StackOfSegs.top();
		StackOfSegs.pop();

		// Check for overlap with each existing edge

		for( int j = 0; j < edges.size(); ++j ) {

			if( IsSubseg( StackOfSegs, edges[j], test ) ) {

				// edge[j] is wholly contained in 'test' so
				// we remove edge[j] ourselves. InSubseg()
				// has already added the leftover pieces
				// to the stack for the next pass.

				edges.erase( edges.begin() + j );
				goto next_seg;
			}
		}

		// No overlaps -- add to figure

		edges.push_back( test );

next_seg:;
	}
}

/* --------------------------------------------------------------- */
/* OffsetControlPoints ------------------------------------------- */
/* --------------------------------------------------------------- */

static void OffsetControlPoints(
	vector<vertex>	&ctl,
	const IBox		&B )
{
	int	dx = B.L - 1, dy = B.B - 1;
	int	nc = ctl.size();

	for( int i = 0; i < nc; ++i ) {

		ctl[i].x += dx;
		ctl[i].y += dy;
	}
}

/* --------------------------------------------------------------- */
/* MeshGetBounds ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MeshGetBounds(
	IBox				&B,
	const vector<Point>	&pts,
	FILE*				flog )
{
	BBoxFromPoints( B, pts );

	fprintf( flog,
	"Region size is [%d %d] in x, [%d %d] in y.\n",
	B.L, B.R, B.B, B.T );
}

/* --------------------------------------------------------------- */
/* MeshMakeSingleTri --------------------------------------------- */
/* --------------------------------------------------------------- */

// Return one large triangle (and its three control points)
// on the given bounding box.
//
void MeshMakeSingleTri(
	vector<triangle>	&tri,
	vector<vertex>		&ctl,
	const IBox			&B,
	FILE*				flog )
{
	tri.resize( 1 );
	ctl.resize( 3 );

// set one triangle

	triangle	t;

	t.v[0] = 0;
	t.v[1] = 1;
	t.v[2] = 2;

	tri[0] = t;

// and its three control points

	if( B.R - B.L > B.T - B.B ) {
		// horizontal
		ctl[0] = vertex( B.L, B.B );
		ctl[1] = vertex( B.R, B.B );
		ctl[2] = vertex( (B.L + B.R) / 2, B.T );
	}
	else {
		// vertical
		ctl[0] = vertex( B.L, B.B );
		ctl[1] = vertex( B.R, (B.B + B.T) / 2 );
		ctl[2] = vertex( B.L, B.T );
	}

	fprintf( flog, "Computed simple triangle.\n" );
}

/* --------------------------------------------------------------- */
/* MeshCreate ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Create an approximate bounding polygon for given points.
// We prefer a few pixels of error (perhaps 5 pixels) over
// having too fine a fragmentation of edges. Divide polygon
// into triangles.
//
// Fill in vector of triangles and control points.
//
// Return value = {OK=0, error=1+}.
//
// Notes
// -----
// For the mesh building step, we are at liberty to make any
// triangles we want, and we choose to make triangle vertices
// only at integer pixel coordinates. This helps make various
// calculations about whether a point is on a line or on this
// or that side of a line an exact calculation. That explains
// the peculiar choice of integer (x,y) in class vertex.
//
int MeshCreate(
	vector<triangle>	&tri,
	vector<vertex>		&ctl,
	const vector<Point>	&pts,
	const IBox			&B,
	FILE*				flog )
{
	clock_t	t0 = StartTiming();

	tri.clear();
	ctl.clear();

	int	npts = pts.size();

/* ---------------------- */
/* Set edge length limits */
/* ---------------------- */

	double	lmin, lmax;

	SetEdgeLimits( lmin, lmax, B, flog );

/* ---------------------------------------- */
/* Make bitmap of ones where we have points */
/* ---------------------------------------- */

	vector<uint8>	map;
	int				w, h;

	MakeMap( map, w, h, B, pts );

/* -------------------------- */
/* Find left edge of boundary */
/* -------------------------- */

	int		x0, y0;

	if( !FindLeftEdge( x0, y0, map, w, h, flog ) )
		return 3;

/* --------------------------------- */
/* Get boundary points, then corners */
/* --------------------------------- */

	vector<vertex>	corners;

	{
		vector<vertex>	bndry;

		if( !GetBoundaryPoints( bndry, map, w, h, x0, y0, flog ) )
			return 4;

		FindCorners( corners, bndry, flog );
	}

/* ------------------- */
/* Split long segments */
/* ------------------- */

	SplitLongSides( corners, lmin, flog );

/* --------------------------------------- */
/* Reduce segment count by cutting corners */
/* --------------------------------------- */

	vector<lineseg>	edges;

	CutCorners( edges, lmin, lmax, corners, flog );

/* ---------------------------------- */
/* Repair any segment crossing errors */
/* ---------------------------------- */

	UncrossEdges( edges, flog );

	ReverseWhole( edges, flog );

/* -------------------------------- */
/* Create list of internal vertices */
/* -------------------------------- */

	vector<vertex>	vinside;

	SetInternalVertices( vinside, map, w, h,
		edges, int((lmin+lmax)/2), flog );

/* ----------------------------- */
/* Divide polygon into triangles */
/* ----------------------------- */

// Take the longest remaining edge, one at a time. Find the
// best triangle for that edge, and then remove the triangle
// from the figure. Repeat until no edges remain. The result
// will be a list of control points, and referring triangles.

	double	Aexpect = 0.0;

	for( ; edges.size() > 0; ) {

		/* ---------- */
		/* Check area */
		/* ---------- */

		// At the end of each pass, we calculate the expected
		// area for the next pass which is: Alast - Atri.
		// We skip expected area test for 1st pass.

		double	area = AreaOfPolygon( edges );

		if( area < 0.0 ||
//			(tri.size() && fabs( area - Aexpect ) > 0.001) ) {
			(tri.size() && fabs( area - Aexpect ) > 1.000) ) {

			fprintf( flog,
			"Internal error! Area %f < 0, or not expected %f\n",
			area, Aexpect );

			return 5;
		}

		/* -------------------- */
		/* Report current state */
		/* -------------------- */

		ListEdgesMatlab( edges, flog );

		fprintf( flog, "\nEdges %ld; Area %f\n", edges.size(), area );

		/* ------------------- */
		/* Remove longest edge */
		/* ------------------- */

		int which = LongestEdge( edges );

		vertex	va( edges[which].v[0] ),
				vb( edges[which].v[1] );

		edges.erase( edges.begin() + which );

		fprintf( flog,
		"\nWorking on edge %d; (%d %d) -> (%d %d).\n",
		which, va.x, va.y, vb.x, vb.y );

		/* --------------------------------------- */
		/* Seek best triangle vertex for this edge */
		/* --------------------------------------- */

		int type, indx;

		BestVertex( type, indx, va, vb, edges, vinside, flog );

		/* ---------------- */
		/* If none found... */
		/* ---------------- */

		if( type == -2 ) {

			fprintf( flog, "\nNo legal triangle at all??"
			" %ld triangles so far, area limit %d\n",
			tri.size(), GBL.mch.MTA );

			if( tri.size() > 0 ) {

				// We'll call it good if we've got
				// at least one triangle already.

				break;
			}

			// There are no triangles to be found, so...
			// we'll make one up from the bounding box.

			fprintf( flog, "STAT: Fall back to single triangle.\n" );

			MeshMakeSingleTri( tri, ctl, B, flog );

			goto exit;
		}

		/* --------------------------------- */
		/* Construct triangle va->vb->vc->va */
		/* --------------------------------- */

		vertex	vc = (type < 0 ? vinside[indx] : edges[indx].v[type]);
		double	Atri = AreaOfTriangle( va, vb, vc );

		fprintf( flog,
		"Triangle (%d %d) (%d %d) (%d %d); area %f\n",
		va.x, va.y, vb.x, vb.y, vc.x, vc.y, Atri );

		AddTriangle( tri, ctl, va, vb, vc );

		/* -------------------- */
		/* Update list of edges */
		/* -------------------- */

		UpdatePolygon( edges, va, vb, vc );

		/* -------------- */
		/* Update vinside */
		/* -------------- */

		// If the edge came from the internal list,
		// it's internal no longer, so remove it.

		if( type == -1 )
			vinside.erase( vinside.begin() + indx );

		/* -------------------- */
		/* Update expected area */
		/* -------------------- */

		Aexpect = area - Atri;
	}

/* --------------------------------- */
/* Convert back to input coordinates */
/* --------------------------------- */

	OffsetControlPoints( ctl, B );

/* -------------- */
/* Report success */
/* -------------- */

exit:
	fprintf( flog,
	"\nSTAT: From %d pts, got %ld triangles, %ld control points.\n",
	npts, tri.size(), ctl.size() );

	StopTiming( flog, "MeshCreate", t0 );

	return 0;
}

/* --------------------------------------------------------------- */
/* MeshCreateX --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Divide bounding box B into regular array of triangles.
// Fill in vector of triangles and control points.
//
// Return 0 always (no error).
//
int MeshCreateX(
	vector<triangle>	&tri,
	vector<vertex>		&ctl,
	const vector<Point>	&pts,
	const IBox			&B,
	FILE*				flog )
{
	clock_t	t0 = StartTiming();

	tri.clear();
	ctl.clear();

/* --------------- */
/* Regular spacing */
/* --------------- */

	double	Dx	= GBL.mch.MNL,
			Dy	= GBL.mch.MNL;
	int		Lx	= B.R - B.L,
			Ly	= B.T - B.B,
			Nx,
			Ny;

	if( !Dx || Dx >= Lx ) {
		Nx = 1;
		Dx = Lx;
	}
	else {
		Nx = int(ceil( Lx / Dx ));
		Dx = Lx / Nx;
	}

	if( !Dy || Dy >= Ly ) {
		Ny = 1;
		Dy = Ly;
	}
	else {
		Ny = int(ceil( Ly / Dy ));
		Dy = Ly / Ny;
	}

/* ----------------------------------------- */
/* Reduce tri count while Atri < GBL.mch.MTA */
/* ----------------------------------------- */

	if( GBL.A.z != GBL.B.z ) {

		while( Nx*Ny > 1 && (Lx*Ly) / (Nx*Ny * 2) < GBL.mch.MTA ) {

			if( Nx >= Ny )
				Dx = Lx / --Nx;
			else
				Dy = Ly / --Ny;
		}
	}

/* ----------------------- */
/* Report basic grid specs */
/* ----------------------- */

	fprintf( flog, "Lx Dx Nx: %5d %8.2f %3d\n", Lx, Dx, Nx );
	fprintf( flog, "Ly Dy Ny: %5d %8.2f %3d\n", Ly, Dy, Ny );

/* ---------- */
/* Create ctl */
/* ---------- */

	for( int iy = 0; iy <= Ny; ++iy ) {

		int	y = (iy < Ny ? B.B + int(iy*Dy) : B.T);

		for( int ix = 0; ix <= Nx; ++ix ) {

			int	x = (ix < Nx ? B.L + int(ix*Dx) : B.R);

			ctl.push_back( vertex( x, y ) );
		}
	}

/* ---------- */
/* Create tri */
/* ---------- */

	int	w = Nx + 1;

	for( int iy = 0; iy < Ny; ++iy ) {

		for( int ix = 0; ix < Nx; ++ix ) {

			triangle	t;

			t.v[0] = ix + w * iy;
			t.v[1] = t.v[0] + 1;
			t.v[2] = t.v[0] + w;
			tri.push_back( t );

			t.v[0] = t.v[2];
			t.v[2] = t.v[0] + 1;
			tri.push_back( t );
		}
	}

/* ---------------- */
/* Remove empty tri */
/* ---------------- */

	const double	occ = 0.30;

	int	ntri = tri.size(),
		npnt = pts.size();

// map pts into their triangles

	vector<int>	in( ntri, 0 );

	for( int i = 0; i < npnt; ++i ) {

		vertex	v( int(pts[i].x), int(pts[i].y) );

		for( int j = 0; j < ntri; ++j ) {

			const triangle& T = tri[j];

			if( InTriangle(
				ctl[T.v[0]], ctl[T.v[1]], ctl[T.v[2]], v ) ) {

				++in[j];
				break;
			}
		}
	}

// remove tri with low occupancy

	for( int i = ntri - 1; i >= 0; --i ) {

		const triangle&	T = tri[i];

		if( !in[i] ||
			in[i] <= occ * AreaOfTriangle(
			ctl[T.v[0]], ctl[T.v[1]], ctl[T.v[2]] ) ) {

			tri.erase( tri.begin() + i );
			--ntri;
		}
	}

/* ----------------------- */
/* Remove unreferenced ctl */
/* ----------------------- */

	if( ntri < in.size() ) {

		fprintf( flog,
		"\nOf %ld triangles, %d were above %3d%% occupancy.\n",
		in.size(), ntri, int(occ*100.0) );

		for( int i = ctl.size() - 1; i >= 0; --i ) {

			for( int j = 0; j < ntri; ++j ) {

				const triangle&	T = tri[j];

				if( T.v[0] == i || T.v[1] == i || T.v[2] == i )
					goto next_i;
			}

			// not ref'd

			ctl.erase( ctl.begin() + i );

			for( int j = 0; j < ntri; ++j ) {

				triangle&	T = tri[j];

				if( T.v[0] > i ) --T.v[0];
				if( T.v[1] > i ) --T.v[1];
				if( T.v[2] > i ) --T.v[2];
			}
next_i:;
		}
	}

/* ------------ */
/* Final report */
/* ------------ */

	fprintf( flog,
	"\nSTAT: From %ld pts, got %ld triangles, %ld control points.\n",
	pts.size(), tri.size(), ctl.size() );

	StopTiming( flog, "MeshCreate", t0 );

	return 0;
}


