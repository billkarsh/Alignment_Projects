

#include	"Geometry.h"
#include	"Maths.h"
#include	"TAffine.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// compute (v1-O) x (v2-O)
//
#define	CROSS( O, v1, v2 )										\
	((v1.x - O.x)*(v2.y - O.y) - (v2.x - O.x)*(v1.y - O.y))






/* --------------------------------------------------------------- */
/* SegPointDist -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return closest approach of point 2 to any point
// on closed line segment from point 0 to point 1.
//
double SegPointDist(
	int		x0,
	int		y0,
	int		x1,
	int		y1,
	int		x2,
	int		y2 )
{
// make (x0,y0) the origin;
// form vector A from 0 to 1;
// form vector B from 0 to 2.

	x1	-= x0;
	x2	-= x0;
	y1	-= y0;
	y2	-= y0;

	int	BdotA = x1*x2 + y1*y2;

	if( BdotA <= 0 ) {

		// B closest to origin

		return sqrt( x2*x2 + y2*y2 );
	}

	int	AdotA = x1*x1 + y1*y1;

	if( BdotA >= AdotA ) {

		// B closest to A

		x2	-= x1;
		y2	-= y1;

		return sqrt( x2*x2 + y2*y2 );
	}

// otherwise, dist is perpendicular component: |AxB|/|A|

	return abs( double(x1*y2 - x2*y1) ) / sqrt( AdotA );
}

/* --------------------------------------------------------------- */
/* BBoxFromPoints ------------------------------------------------ */
/* --------------------------------------------------------------- */

void BBoxFromPoints( IBox &B, const vector<Point> &pts )
{
	double	xlo, xhi, ylo, yhi;
	int		npts = pts.size();

	xlo = xhi = pts[0].x;
	ylo = yhi = pts[0].y;

	for( int i = 1; i < npts; ++i ) {

		double	t;

		if( (t = pts[i].x) < xlo )
			xlo = t;
		else if( t > xhi )
			xhi = t;

		if( (t = pts[i].y) < ylo )
			ylo = t;
		else if( t > yhi )
			yhi = t;
	}

	B.L	= (int)floor( xlo );
	B.R	= (int)ceil( xhi );

	B.B	= (int)floor( ylo );
	B.T	= (int)ceil( yhi );
}


void BBoxFromPoints( DBox &B, const vector<Point> &pts )
{
	int	npts = pts.size();

	B.L = B.R = pts[0].x;
	B.B = B.T = pts[0].y;

	for( int i = 1; i < npts; ++i ) {

		double	t;

		if( (t = pts[i].x) < B.L )
			B.L = t;
		else if( t > B.R )
			B.R = t;

		if( (t = pts[i].y) < B.B )
			B.B = t;
		else if( t > B.T )
			B.T = t;
	}

	B.L	= floor( B.L );
	B.R	= ceil( B.R );

	B.B	= floor( B.B );
	B.T	= ceil( B.T );
}

/* --------------------------------------------------------------- */
/* BoxesFromShifts ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Suppose two images, each starting at (0,0);
// Image 1 is w1xh1 and image 2 is w2xh2.
// Now let image 1 be shifted by adding (x,y) to its coordinates.
// Calculate the overlap region in each image's own coordinates.
//
// Note: Box bottom & top are exchanged w.r.t. most conventions...
// Here, bottom has lower y-value, so that (t - b) is positive.
//
void BoxesFromShifts(
	IBox	&B1,
	IBox	&B2,
	int		w1,
	int		h1,
	int		w2,
	int		h2,
	int		x,
	int		y )
{
// after shifting, bbox for image2 will be
// [x, x+w-1] in x, and [y, y+h-1] in y

	B2.L = max( 0, x );			// greater of left edges
	B2.R = min( w2-1, x+w1-1 );	// least of the rights
	B2.B = max( 0, y );			// greater of the bottoms
	B2.T = min( h2-1, y+h1-1 );	// lesser of the tops
	B2.R = max( B2.R, B2.L );	// non-negative width
	B2.T = max( B2.T, B2.B );	// non-negative height

// bounding box relative to image 1 is simply offset

	B1.L = B2.L - x;
	B1.R = B2.R - x;
	B1.B = B2.B - y;
	B1.T = B2.T - y;
}

/* --------------------------------------------------------------- */
/* TightestBBox -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Try all integer angle rotations [-45,45] to determine the
// orientation of the point set having smallest bbox.
//
// Return that bbox and corresponding angle.
//
int TightestBBox( DBox &B, const vector<Point> &pts )
{
	int	np = pts.size();

/* --------------------- */
/* Create region outline */
/* --------------------- */

// For each y-row get the min and max x-coord.

	vector<Point>	outline;

	BBoxFromPoints( B, pts );

	{
		int				ny = int(B.T - B.B) + 1;
		vector<double>	minx( ny, B.R + 1 );
		vector<double>	maxx( ny, B.L - 1 );

		for( int i = 0; i < np; ++i ) {

			int	iy = int(floor( pts[i].y - B.B ));

			minx[iy] = fmin( minx[iy], pts[i].x );
			maxx[iy] = fmax( maxx[iy], pts[i].x );
		}

		for( int iy = 0; iy < ny; ++iy ) {

			if( minx[iy] <= maxx[iy] ) {

				outline.push_back( Point( minx[iy], iy + B.B ) );
				outline.push_back( Point( maxx[iy], iy + B.B ) );
			}
		}
	}

/* ------------------------------------------- */
/* Outline only useful if fewer pts than whole */
/* ------------------------------------------- */

	const vector<Point>	*ppointset = &outline;

	if( outline.size() >= np ) {

		ppointset = &pts;
		outline.clear();
	}

/* --------------- */
/* Find best angle */
/* --------------- */

	double	best_area	= (B.R - B.L) * (B.T - B.B);
	int		best_angle	= 0;

	for( int angle = -45; angle <= 45; ++angle ) {

		if( angle == 0 )
			continue;

		vector<Point>	P = *ppointset;
		TAffine			T;
		DBox			box;
		double			area;

		T.NUSetRot( angle*PI/180 );
		T.Apply_R_Part( P );
		BBoxFromPoints( box, P );

		area = (box.R - box.L) * (box.T - box.B);

		if( area < best_area ) {

			B			= box;
			best_area	= area;
			best_angle	= angle;
		}
	}

	return best_angle;
}

/* --------------------------------------------------------------- */
/* Propagate ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Modify given vector of doubles (v) as follows...
// Starting at v[first], take all connected points whose values
// fall between tmin and tmax, and set those values to set_to.
//
// Return count of qualified points (region area).
//
int Propagate(
	vector<double>	&v,
	int				w,
	int				h,
	int				first,
	double			tmin,
	double			tmax,
	double			set_to )
{
	stack<int>	st;
	int			cnt = 0;

	st.push( first );

	while( !st.empty() ) {

		int		j = st.top();

		st.pop();

		if( tmin < v[j] && v[j] < tmax ) {

			int		y = j / w;
			int		x = j - w * y;

			v[j] = set_to;	// remove from further consideration
			++cnt;

			// push the four neighbors
			if( x - 1 >= 0 )	st.push( j - 1 );
			if( x + 1 <  w )	st.push( j + 1 );
			if( y - 1 >= 0 )	st.push( j - w );
			if( y + 1 <  h )	st.push( j + w );
		}
	}

	return cnt;
}


// Fill vector of points (plist), by examining values (v)...
// Starting at v[first], take all connected points whose values
// exceed thresh and add them to (plist).
//
// Note that this function will modify (v) by setting the values
// for qualified points to the value set_to. This is a necessary
// consequence of bookkeeping.
//
// Return count of qualified points (region area).
//
int Propagate(
	vector<Point>	&plist,
	vector<double>	&v,
	int				w,
	int				h,
	int				first,
	double			thresh,
	double			set_to )
{
	stack<int>	st;
	int			cnt = 0;

	plist.clear();
	plist.reserve( w * h );

	st.push( first );

	while( !st.empty() ) {

		int		j = st.top();

		st.pop();

		if( v[j] > thresh ) {

			int		y = j / w;
			int		x = j - w * y;

			plist.push_back( Point( x, y ) );

			v[j] = set_to;	// remove from further consideration
			++cnt;

			// push the four neighbors
			if( x - 1 >= 0 )	st.push( j - 1 );
			if( x + 1 <  w )	st.push( j + 1 );
			if( y - 1 >= 0 )	st.push( j - w );
			if( y + 1 <  h )	st.push( j + w );
		}
	}

	return cnt;
}

/* --------------------------------------------------------------- */
/* MapBlobRng ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given image <I> and empty <map>, create connected object
// by placing a one into map at each pixel in I whose value
// is in range [tmin, tmax]. Starts at seed point <first>.
//
// Return object pixel count (area).
//
int MapBlobRng(
	vector<uint8>			&map,
	const vector<double>	&I,
	int						w,
	int						h,
	int						first,
	double					tmin,
	double					tmax )
{
	stack<int>	st;
	int			cnt = 0;

	st.push( first );

	while( !st.empty() ) {

		int		j = st.top();

		st.pop();

		if( !map[j] && tmin <= I[j] && I[j] <= tmax ) {

			int		y = j / w;
			int		x = j - w * y;

			map[j] = 1;
			++cnt;

			// push the four neighbors
			if( x - 1 >= 0 )	st.push( j - 1 );
			if( x + 1 <  w )	st.push( j + 1 );
			if( y - 1 >= 0 )	st.push( j - w );
			if( y + 1 <  h )	st.push( j + w );
		}
	}

	return cnt;
}


int MapBlobRng(
	vector<uint8>	&map,
	const uint8		*I,
	int				w,
	int				h,
	int				first,
	int				tmin,
	int				tmax )
{
	stack<int>	st;
	int			cnt = 0;

	st.push( first );

	while( !st.empty() ) {

		int		j = st.top();

		st.pop();

		if( !map[j] && tmin <= I[j] && I[j] <= tmax ) {

			int		y = j / w;
			int		x = j - w * y;

			map[j] = 1;
			++cnt;

			// push the four neighbors
			if( x - 1 >= 0 )	st.push( j - 1 );
			if( x + 1 <  w )	st.push( j + 1 );
			if( y - 1 >= 0 )	st.push( j - w );
			if( y + 1 <  h )	st.push( j + w );
		}
	}

	return cnt;
}

/* --------------------------------------------------------------- */
/* MapBlobVar ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given image <I> and empty <map>, create connected object
// by placing a one into map at each pixel in I whose value
// is within <tol> of the average of its <size> x <size>
// neighbors. Starts at seed <first>.
//
// Return object pixel count (area).
//
int MapBlobVar(
	vector<uint8>	&map,
	const uint8		*I,
	int				w,
	int				h,
	int				first,
	int				size,
	int				tol )
{
	stack<int>	st;
	int			cnt = 0;

	size /= 2;	// halfwidth

	st.push( first );

	while( !st.empty() ) {

		int	j = st.top();

		st.pop();

		if( map[j] )
			continue;

		int	y0 = j / w,
			x0 = j - w * y0,
			A = -I[j],
			n = -1,
			L, R, B, T;

		L = max( x0 - size, 0 );
		R = min( x0 + size, w - 1 );
		B = max( y0 - size, 0 );
		T = min( y0 + size, h - 1 );

		for( int y = B; y <= T; ++y ) {
			for( int x = L; x <= R; ++x ) {
				A += I[x+w*y];
				++n;
			}
		}

		if( iabs( I[j] - A/n ) <= tol ) {

			map[j] = 1;
			++cnt;

			// push the four neighbors
			if( x0 - 1 >= 0 )	st.push( j - 1 );
			if( x0 + 1 <  w )	st.push( j + 1 );
			if( y0 - 1 >= 0 )	st.push( j - w );
			if( y0 + 1 <  h )	st.push( j + w );
		}
	}

	return cnt;
}

/* --------------------------------------------------------------- */
/* DilateMap1Pix ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Dilate bitmap <map> 1 pixel.
//
// Version is sloppy in that it does not act on 1-pixel border.
//
void DilateMap1Pix( vector<uint8> &map, int w, int h )
{
	vector<uint8>	org = map;
	int				N	= w * h - w - 1;

	for( int i = w + 1; i < N; ++i ) {

		map[i] |= org[i - 1];
		map[i] |= org[i + 1];
		map[i] |= org[i - w];
		map[i] |= org[i + w];
	}
}

/* --------------------------------------------------------------- */
/* ErodeMap1Pix -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Erode bitmap <map> 1 pixel.
//
// Version is sloppy in that it does not act on 1-pixel border.
//
void ErodeMap1Pix( vector<uint8> &map, int w, int h )
{
	vector<uint8>	org = map;
	int				N	= w * h - w - 1;

	for( int i = w + 1; i < N; ++i ) {

		map[i] &= org[i - 1];
		map[i] &= org[i + 1];
		map[i] &= org[i - w];
		map[i] &= org[i + w];
	}
}

/* --------------------------------------------------------------- */
/* PixelListFromPolygon ------------------------------------------ */
/* --------------------------------------------------------------- */

// Compile a pixel list from a polygon.
// Polygon is defined by the points in vector Pgon,
// with the end point wrapping around.
//
void PixelListFromPolygon(
	vector<Point>			&Plist,
	const vector<Point>		&Pgon )
{
	double	xmin, ymin, xmax, ymax;
	int		i, x, y, np = Pgon.size();

	xmin = xmax = Pgon[0].x;
	ymin = ymax = Pgon[0].y;

	for( i = 1; i < np; ++i ) {

		double	t;

		if( (t = Pgon[i].x) < xmin )
			xmin = t;
		else if( t > xmax )
			xmax = t;

		if( (t = Pgon[i].y) < ymin )
			ymin = t;
		else if( t > ymax )
			ymax = t;
	}

	Plist.clear();
	Plist.reserve( (size_t)((xmax - xmin) * (ymax - ymin)) );

	for( x = int(xmin-1); x <= xmax+1; ++x ) {

		for( y = int(ymin-1); y <= ymax+1; ++y ) {

			// is (x,y) inside?
			// Draw a line to infinity, count the crossings to see

			int		nint = 0;

			for( i = 0; i < np; ++i ) {

				// line goes from Pgon[i] to Pgon[i2]
				int		i2 = (i+1)%np;

				// ignore horizontal lines
				if( Pgon[i].y == Pgon[i2].y )
					continue;

				// entirely above; ignore it
				if( Pgon[i].y >= y && Pgon[i2].y >= y )
					continue;

				// entirely below; ignore it
				if( Pgon[i].y < y && Pgon[i2].y < y )
					continue;

				// otherwise find x intercept
				double	xx = Pgon[i].x +
							(y-Pgon[i].y)/(Pgon[i2].y-Pgon[i].y)
							*(Pgon[i2].x-Pgon[i].x);

				if( xx > x )
					++nint;
			}

			if( nint & 1 )
				Plist.push_back( Point( x, y ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* ImageFromValuesAndPoints -------------------------------------- */
/* --------------------------------------------------------------- */

// Create image from point list.
// Image's new origin is (0,0).
//
void ImageFromValuesAndPoints(
	vector<double>			&I,
	int						w,
	int						h,
	const vector<double>	&v,
	const vector<Point>		&plist,
	int						xmin,
	int						ymin )
{
// Size and zero destination image

	I.assign( w * h, 0.0 );

// Paint with points

	int	np = plist.size();

	for( int i = 0; i < np; ++i ) {

		double	x = plist[i].x - xmin,
				y = plist[i].y - ymin;

		DistributePixel( x, y, v[i], I, w, h );
	}
}

/* --------------------------------------------------------------- */
/* ValuesFromImageAndPoints -------------------------------------- */
/* --------------------------------------------------------------- */

// Takes a real-valued vector of pixel locations,
// returns a vector of pixel values.
// Uses bi-linear interpolation of adjacent pixels.
//
void ValuesFromImageAndPoints(
	vector<double>			&v,
	const uint8*			raster,
	int						w,
	const vector<Point>		&plist )
{
	int		nPts = plist.size();

	v.resize( nPts );

	for( int i = 0; i < nPts; ++i )
		v[i] = InterpolatePixel( plist[i].x, plist[i].y, raster, w );
}


// Takes a real-valued vector of pixel locations,
// returns a vector of pixel values.
// Uses bi-linear interpolation of adjacent pixels.
//
void ValuesFromImageAndPoints(
	vector<double>			&v,
	double					&dvx,
	double					&dvy,
	const vector<double>	&raster,
	int						w,
	const vector<Point>		&plist,
	const vector<double>	&spv )
{
	double	derx = 0.0, dery = 0.0; // local copies for derivatives
	int		nPts = plist.size();

	v.resize( nPts );

	for( int i = 0; i < nPts; ++i ) {

		// find the 4 points surrounding real valued coordinate

		if( plist[i].x <  0.0	||
			plist[i].x >= w - 1	||
			plist[i].y <  0.0	||
			plist[i].y >= 4095 ) {

			// anything outside the raster is 0.0
			v[i] = 0.0;
			continue;
		}

		int xl	= (int)plist[i].x;
		int yl	= (int)plist[i].y;
		int xr	= xl + 1;
		int yu	= yl + 1;

		double alpha	= plist[i].x - xl;
		double beta		= plist[i].y - yl;
		double ll		= raster[w*yl + xl];
		double ul		= raster[w*yu + xl];
		double ur		= raster[w*yu + xr];
		double lr		= raster[w*yl + xr];
		double dx		= lr - ll;
		double dy		= ul - ll;
		double dxy		= ll - lr - ul + ur;
		double d		= ll + alpha*dx + beta*dy + alpha*beta*dxy;

		derx += spv[i] * (dx + beta*dxy);
		dery += spv[i] * (dy + alpha*dxy);

		v[i] = d;
	}

	dvx = derx;
	dvy = dery;
	printf( "ValsFrmImg: Local derivatives %f %f\n", derx, dery );
}

/* --------------------------------------------------------------- */
/* IDistSqr ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int vertex::IDistSqr( const vertex& rhs ) const
{
	int	dx = x - rhs.x,
		dy = y - rhs.y;

	return dx*dx + dy*dy;
}

/* --------------------------------------------------------------- */
/* DistSqr ------------------------------------------------------- */
/* --------------------------------------------------------------- */

double vertex::DistSqr( const vertex& rhs ) const
{
	double	dx = x - rhs.x,
			dy = y - rhs.y;

	return dx*dx + dy*dy;
}

/* --------------------------------------------------------------- */
/* SegPointDist -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return closest approach of point 2 to any point
// on closed line segment from point 0 to point 1.
//
double SegPointDist(
	const vertex	&v0,
	const vertex	&v1,
	const vertex	&v2 )
{
	return SegPointDist( v0.x, v0.y, v1.x, v1.y, v2.x, v2.y );
}

/* --------------------------------------------------------------- */
/* LeftSide ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return true if point c is on the left side of vector a -> b.
//
bool LeftSide( const Point &a, const Point &b, const Point &c )
{
	return CROSS( a, b, c ) > 0;
}


// Return true if point c is on the left side of vector a -> b.
//
bool LeftSide( const vertex &a, const vertex &b, const vertex &c )
{
	return CROSS( a, b, c ) > 0;
}

/* --------------------------------------------------------------- */
/* LinesCross ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if the infinite lines passing through points {1,2}
// and through {3,4} intersect.
//
bool LinesCross(
	const vertex	&p1,
	const vertex	&p2,
	const vertex	&p3,
	const vertex	&p4 )
{
	int		x21	= p2.x - p1.x,
			x43	= p4.x - p3.x,
			y21	= p2.y - p1.y,
			y43	= p4.y - p3.y;

// cross if 1-2, 3-4 diff slopes

	if( x21*y43 != x43*y21 )
		return true;

// else if same slope, cross (colinear) if 1-3, 1-2 same slope

	int		x31	= p3.x - p1.x,
			y31	= p3.y - p1.y;

	return x21*y31 == x31*y21;
}

/* --------------------------------------------------------------- */
/* ClosedSegIsects ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Test if closed line segments [1,2] and [3,4] (endpoints
// included) intersect. Return {count; intersection pts.}.
// Possible results are {0}, {1;pi}, {2;pi,pj}. Two points
// are returned if [1,2] and [3,4] are colinear and their
// overlap is the non-zero-length segment [i,j].
//
// Set A = vector from 1 to 2 and B = vector from 3 to 4.
//
// Represent point Pt on 1-2 as: P1 + tA, t range [0,1].
// Represent point Ps on 3-4 as: P3 + sB, s range [0,1].
//
// At an intersection point, Ps lies on 1-2, so (Ps-P1)xA = 0.
// At an intersection point, Pt lies on 3-4, so (Pt-P3)xB = 0.
//
// s = (P1-P3)xA / BxA.
// t = (P1-P3)xB / BxA.
//
// Note that args p3, p4 are passed by value here (copies of
// caller data) because we might need to swap them.
//
int ClosedSegIsects(
	vertex			&pi,
	vertex			&pj,
	const vertex	&p1,
	const vertex	&p2,
	vertex			p3,
	vertex			p4 )
{
	int	x12	= p2.x - p1.x,
		y12	= p2.y - p1.y,
		x34	= p4.x - p3.x,
		y34	= p4.y - p3.y,
		BxA	= x34*y12 - x12*y34;

	if( BxA ) {	// nondegenerate case

		int	x31	= p1.x - p3.x,
			y31	= p1.y - p3.y;

		double	s = (x31*y12 - x12*y31) / (double)BxA;

		if( s < 0.0 || s > 1.0 )
			return 0;

		double	t = (x31*y34 - x34*y31) / (double)BxA;

		if( t < 0.0 || t > 1.0 )
			return 0;

		pi.x = int(p3.x + s * x34);
		pi.y = int(p3.y + s * y34);

		return 1;
	}
	else {	// degenerate

		// Since BxA == 0, we may have:
		// (1) |A| == 0 or |B| == 0,
		// (2) A parallel transported from B,
		// (3) A, B colinear, nonoverlapping,
		// (4) A, B colinear and overlapping.

		// Handle case (1).

		int	AdotA = x12*x12 + y12*y12,
			BdotB = x34*x34 + y34*y34;

		if( !AdotA && !BdotB ) {

			if( p1 == p3 ) {
				pi = p1;
				return 1;
			}

			return 0;
		}
		else if( !AdotA ) {

			double	d;

			d = SegPointDist( p3.x, p3.y, p4.x, p4.y, p1.x, p1.y );

			if( d < 0.001 ) {
				pi = p1;
				return 1;
			}

			return 0;
		}
		else if( !BdotB ) {

			double	d;

			d = SegPointDist( p1.x, p1.y, p2.x, p2.y, p3.x, p3.y );

			if( d < 0.001 ) {
				pi = p3;
				return 1;
			}

			return 0;
		}

		// Non-zero lengths!
		// To get intersection we first need colinearity.
		// We already know the slopes are the same. We next
		// need same intercepts. Intercept (b) for segment
		// 1-2 is y2 - y12*x2/x12. However, if x12 is zero
		// the vertical segments need same X to be colinear.

		if( x12 ) {	// non-vertical

			double	b = p2.y - (double)y12*p2.x/x12;

			if( b != p4.y - (double)y34*p4.x/x34 )
				return 0;
		}
		else {	// vertical

			if( p1.x != p3.x )
				return 0;
		}

		// Colinear!
		// Get [3,4] oriented same as [1,2]:
		// swap 3, 4 if dot product negative.

		if( x12*x34 + y12*y34 < 0 ) {

			pi = p3;
			p3 = p4;
			p4 = pi;
		}

		// Colinear and same orientation...but how are
		// they ordered? If [1,2] goes opposite [1,4]
		// there is no overlap, if zero, the endpoints
		// match, if positive then p4 is beyond p1 but
		// we don't know how far yet.

		int	x14	= p4.x - p1.x,
			y14	= p4.y - p1.y,
			dot;

		dot = x12*x14 + y12*y14;

		if( dot < 0 )
			return 0;

		if( !dot ) {
			pi = p1;
			return 1;
		}

		// Similarly compare [3,2] to [3,4].

		int	x32	= p2.x - p3.x,
			y32	= p2.y - p3.y;

		dot = x34*x32 + y34*y32;

		if( dot < 0 )
			return 0;

		if( !dot ) {
			pi = p2;
			return 1;
		}

		// We will return a non-zero-length segment!
		// Thus far, along the common direction [1,2],
		// we know that p4 > {p1,p3} and p2 > {p1,p3}.
		//
		// pi will either be p1 or p3, and pj will be
		// whichever of {p2,p4} is closer to pi.

		int	x13	= p3.x - p1.x,
			y13	= p3.y - p1.y;

		dot = x12*x13 + y12*y13;

		if( dot <= 0 ) {

			pi = p1;

			if( x12*x12 + y12*y12 >= x14*x14 + y14*y14 )
				pj = p4;
			else
				pj = p2;
		}
		else {

			pi = p3;

			if( x34*x34 + y34*y34 >= x32*x32 + y32*y32 )
				pj = p2;
			else
				pj = p4;
		}

		return 2;
	}
}

/* --------------------------------------------------------------- */
/* OpenSegsCross ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if open line segments (1,2) and (3,4) have
// an interior intersection point (excluding endpoints).
//
// Set A = vector from 1 to 2 and B = vector from 3 to 4.
//
// Represent point Pt on 1-2 as: P1 + tA, t range (0,1).
// Represent point Ps on 3-4 as: P3 + sB, s range (0,1).
//
// At an intersection point, Ps lies on 1-2, so (Ps-P1)xA = 0.
// At an intersection point, Pt lies on 3-4, so (Pt-P3)xB = 0.
//
// s = (P1-P3)xA / BxA.
// t = (P1-P3)xB / BxA.
//
bool OpenSegsCross(
	const vertex	&p1,
	const vertex	&p2,
	const vertex	&p3,
	const vertex	&p4 )
{
	int	x12	= p2.x - p1.x,
		y12	= p2.y - p1.y,
		x34	= p4.x - p3.x,
		y34	= p4.y - p3.y,
		BxA	= x34*y12 - x12*y34;

	if( BxA ) {	// nondegenerate case

		int	x31	= p1.x - p3.x,
			y31	= p1.y - p3.y;

		double	s = (x31*y12 - x12*y31) / (double)BxA;

		if( s <= 0.0 || s >= 1.0 )
			return false;

		double	t = (x31*y34 - x34*y31) / (double)BxA;

		return t > 0.0 && t < 1.0;
	}
	else {	// degenerate

		// Since BxA == 0, we may have:
		// (1) |A| == 0 or |B| == 0,
		// (2) A parallel transported from B,
		// (3) A, B colinear, nonoverlapping,
		// (4) A, B colinear and overlapping.

		// In case (1) a zero-length segment IS its endpoints,
		// and endpoint touching is not counted as intersecting
		// for open segments, so the result is always false.

		int	AdotA = x12*x12 + y12*y12,
			BdotB = x34*x34 + y34*y34;

		if( !AdotA || !BdotB )
			return false;

		// Non-zero lengths!
		// To get intersection we first need colinearity.
		// We already know the slopes are the same. We next
		// need same intercepts. Intercept (b) for segment
		// 1-2 is y2 - y12*x2/x12. However, if x12 is zero
		// the vertical segments need same X to be colinear.

		if( x12 ) {	// non-vertical

			double	b = p2.y - (double)y12*p2.x/x12;

			if( b != p4.y - (double)y34*p4.x/x34 )
				return false;
		}
		else {	// vertical

			if( p1.x != p3.x )
				return false;
		}

		// Colinear!
		// Now we have to take a segment A and see if one of
		// segment B's endpoints falls between. But there is
		// a tricky case. Suppose A is shorter and sits entirely
		// interior to B. Neither of B's ends is inside A so we
		// would wrongly conclude no overlap. Therefore we need
		// the longer segment.

		if( AdotA >= BdotB ) {

			// test B ends against seg A

			int	x13	= p3.x - p1.x,
				y13	= p3.y - p1.y,
				dot;

			dot = x12*x13 + y12*y13;

			if( dot > 0 && dot < AdotA )
				return true;

			int	x14	= p4.x - p1.x,
				y14	= p4.y - p1.y;

			dot = x12*x14 + y12*y14;

			if( dot > 0 && dot < AdotA )
				return true;
		}
		else {

			// test A ends against seg B

			int	x31	= p1.x - p3.x,
				y31	= p1.y - p3.y,
				dot;

			dot = x34*x31 + y34*y31;

			if( dot > 0 && dot < BdotB )
				return true;

			int	x32	= p2.x - p3.x,
				y32	= p2.y - p3.y;

			dot = x34*x32 + y34*y32;

			if( dot > 0 && dot < BdotB )
				return true;
		}

		// no interposed ends
		return false;
	}
}

/* --------------------------------------------------------------- */
/* AnyCrossing --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if segment from a to b crosses any other edges.
//
bool AnyCrossing(
	const vector<lineseg>	&s,
	const vertex			&a,
	const vertex			&b )
{
	int		ns = s.size();

	for( int i = 0; i < ns; ++i ) {

		if( OpenSegsCross( s[i].v[0], s[i].v[1], a, b ) )
			return true;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* CountCrossings ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return number of crossings of segment from a to b.
//
int CountCrossings(
	const vector<lineseg>	&s,
	const vertex			&a,
	const vertex			&b )
{
	int		ns = s.size(), N = 0;

	for( int i = 0; i < ns; ++i )
		N += OpenSegsCross( s[i].v[0], s[i].v[1], a, b );

	return N;
}

/* --------------------------------------------------------------- */
/* IsSubseg ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Is an edge 'e' contained within segment 'W'?
// It is if both ends have SegPointDist zero from W.
//
// If contained, push onto stack any non-zero-length
// sections of W that extend beyond e.
//
bool IsSubseg(
	stack<lineseg>	&stk,
	const lineseg	&e,
	const lineseg	&W )
{
	if( 0.001 < SegPointDist( W.v[0], W.v[1], e.v[0] ) ||
		0.001 < SegPointDist( W.v[0], W.v[1], e.v[1] ) ) {

		return false;
	}

/* -------------------------------------------- */
/* It's a subset; push leftover sections to stk */
/* -------------------------------------------- */

// Which e-end is closer to W.v[0]?
// The 'other' e-end pairs with W.v[1].

	int	other, d0 = W.v[0].IDistSqr( e.v[0] );

	if( !d0 ) {
		// nothing to push at this end
		other = 1;
	}
	else {

		int	d1 = W.v[0].IDistSqr( e.v[1] );

		if( d0 < d1 ) {
			stk.push( lineseg( W.v[0], e.v[0] ) );
			other = 1;
		}
		else {

			if( d1 )
				stk.push( lineseg( W.v[0], e.v[1] ) );

			other = 0;
		}
	}

// Now do other end

	if( W.v[1].IDistSqr( e.v[other] ) )
		stk.push( lineseg( e.v[other], W.v[1] ) );

	return true;
}

/* --------------------------------------------------------------- */
/* InTriangle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// True if point p on interior of ABC.
//
bool InTriangle(
	const vertex	&va,
	const vertex	&vb,
	const vertex	&vc,
	const vertex	&p )
{
	return	LeftSide( va, vb, p ) &&
			LeftSide( vb, vc, p ) &&
			LeftSide( vc, va, p );
}

/* --------------------------------------------------------------- */
/* AnyInside ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// For the proposed triangle ABC, return true if any other line
// segments (s) or region-interior points (ips) are on the triangle
// interior.
//
// Note: To test a line segment we test each of its ends, but there
// is a tricky case that this can miss wherein both ends are on the
// boundary but the segment cuts through the interior. We can catch
// this case by testing the segment midpoint. However, we do not want
// to do this midpoint test if both ends are vertices of ABC:
//
// (1) The test is unnecessary since this segment is a boundary.
// (2) The test may be foiled by integer round-off.
//
bool AnyInside(
	const vertex			&va,
	const vertex			&vb,
	const vertex			&vc,
	const vector<lineseg>	&s,
	const vector<vertex>	&ips )
{
	int		i, N;

/* ----------------------- */
/* Check the line segments */
/* ----------------------- */

	N = s.size();

	for( i = 0; i < N; ++i ) {

		const vertex	&v0 = s[i].v[0];
		const vertex	&v1 = s[i].v[1];

		if( InTriangle( va, vb, vc, v0 ) )
			return true;

		if( InTriangle( va, vb, vc, v1 ) )
			return true;

		// if either point is not a vertex, test midpoint

		if( !((v0 == va) || (v0 == vb) || (v0 == vc)) ||
			!((v1 == va) || (v1 == vb) || (v1 == vc)) ) {

			vertex	mid( (v0.x + v1.x)/2, (v0.y + v1.y)/2 );

			if( InTriangle( va, vb, vc, mid ) )
				return true;
		}
	}

/* ------------------------- */
/* Check the internal points */
/* ------------------------- */

	N = ips.size();

	for( i = 0; i < N; ++i ) {

		if( InTriangle( va, vb, vc, ips[i] ) )
			return true;
	}

	return false;	// none inside
}

/* --------------------------------------------------------------- */
/* Area ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

double triangle::Area( const vector<Point> &ctl ) const
{
	return  0.5 * abs( CROSS( ctl[v[0]], ctl[v[1]], ctl[v[2]] ) );
}

/* --------------------------------------------------------------- */
/* AreaOfTriangle ------------------------------------------------ */
/* --------------------------------------------------------------- */

double AreaOfTriangle(
	const Point		&v0,
	const Point		&v1,
	const Point		&v2 )
{
	return 0.5 * abs( CROSS( v0, v1, v2 ) );
}


double AreaOfTriangle(
	const vertex	&v0,
	const vertex	&v1,
	const vertex	&v2 )
{
	return 0.5 * abs( (double)CROSS( v0, v1, v2 ) );
}

/* --------------------------------------------------------------- */
/* BestTriangle -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Find the best triangle.
// If we are inside it's obvious; otherwise pick the closest.
//
int BestTriangle(
	const vector<triangle>	&T,
	const vector<vertex>	&C,
	const Point				&p )
{
	int		nT = T.size();
	vertex	v( int(p.x), int(p.y) );

	for( int i = 0; i < nT; ++i ) {

		const vertex	&v0 = C[T[i].v[0]],
						&v1 = C[T[i].v[1]],
						&v2 = C[T[i].v[2]];

		if( InTriangle( v0, v1, v2, v ) )
			return i;
	}

// Not inside any of them. Find the closest one...

	double	dbest	= BIG;
	int		ibest	= -1;

	for( int i = 0; i < nT; ++i ) {

		const vertex	&v0 = C[T[i].v[0]],
						&v1 = C[T[i].v[1]],
						&v2 = C[T[i].v[2]];
		double			d;

		d = SegPointDist( v0, v1, v );

		if( d < dbest ) {
			dbest = d;
			ibest = i;
		}

		d = SegPointDist( v1, v2, v );

		if( d < dbest ) {
			dbest = d;
			ibest = i;
		}

		d = SegPointDist( v2, v0, v );

		if( d < dbest ) {
			dbest = d;
			ibest = i;
		}

		if( dbest <= 0.1 )
			break;
	}

	return ibest;
}

/* --------------------------------------------------------------- */
/* AreaOfPolygon ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return area of polygon traversed in CCW direction.
//
// The area inside a polygon is just the sum of the areas under
// each of its directed line segments. A line segment can be viewed
// as one side of a tetrahedron, where the other three sides are the
// corresponding length of x-axis, and the two verticals, which are
// the bases in this case. The area "under" a line segment is then
// .5*(y1+y2)*dx. This quantity is positive for a segment directed
// from left to right.
//
// However, conventionally, figures are traced in the CCW direction.
// In that case, the upper-most segments are traced from right to
// left and should contribute positive area to the figure, while the
// lower left-to-right segments count exterior area between the axis
// and the figure which should be subtracted from the sum. We adopt
// the CCW usage convention, so that left-directed segments produce
// positive contributions.
//
double AreaOfPolygon( const vector<lineseg> &edges )
{
	long	A	= 0;
	int		ne	= edges.size();

	for( int i = 0; i < ne; ++i ) {

		const lineseg&	L = edges[i];

		A += (L.v[0].x - L.v[1].x) * (L.v[0].y + L.v[1].y);
	}

	return A / 2.0;
}

/* --------------------------------------------------------------- */
/* RemoveFromMap ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Remove from uint8 map all pixels within 'dist' of vertex.
//
void RemoveFromMap(
	vector<uint8>	&map,
	int				w,
	int				h,
	const vertex	&v,
	int				dist )
{
	int	d	= dist + 1,
		dsqr,
		xlo	= max( v.x - d, 0 ),
		ylo	= max( v.y - d, 0 ),
		xhi	= min( v.x + d, w - 1 ),
		yhi	= min( v.y + d, h - 1 );

	dsqr = d * d;

	for( int y = ylo; y <= yhi; ++y ) {

		int	dy		= y - v.y,
			xmaxsqr	= dsqr - dy*dy;

		for( int x = xlo; x <= xhi; ++x ) {

			int	dx = x - v.x;

			if( dx*dx <= xmaxsqr )
				map[x + w*y] = 0;
		}
	}
}


