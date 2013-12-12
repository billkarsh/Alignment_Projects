

/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* Culled from Lou code ------------------------------------------ */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// ----------------------------------------

class BBox {

public:
	double	xmin, xmax, ymin, ymax;

public:
	BBox()
		{Reset();};

	void Reset()
		{xmin = ymin = BIGD; xmax = ymax = -BIGD;};

	void AddPoint( const Point &p )
		{
			xmin = min( xmin, p.x );
			ymin = min( ymin, p.y );
			xmax = max( xmax, p.x );
			ymax = max( ymax, p.y );
		};
};

/* --------------------------------------------------------------- */
/* IncludeThis --------------------------------------------------- */
/* --------------------------------------------------------------- */

#if 0
// Return true if name in user's 'include only these' list.
//
static bool IncludeThis( const char *name )
{
	int	z	= ZFromName( name ),
		n	= gArgs.TileIDFromName( name ),
		ni	= gArgs.include_only.size();

	for( int i = 0; i < ni; ++i ) {

		if( z == gArgs.include_only[i].z &&
			n == gArgs.include_only[i].id ) {

			return true;
		}
	}

	return false;
}
#endif

/* --------------------------------------------------------------- */
/* IgnoreNames --------------------------------------------------- */
/* --------------------------------------------------------------- */

#if 0
// Return true if line with {name1, name2} should be ignored.
//
static bool IgnoreNames(
	const char				*name1,
	const char				*name2,
	const set<string>		&BadPairs,
	const set<string>		&ign )
{
// If in the bad pairs, return true

	set<string>::iterator	si;

	si = BadPairs.find( string(name1) + string(",") + string(name2) );

	if( si != BadPairs.end() ) {
		//printf( "Reject %s, %s\n", name1, name2 );
		return true;
	}

// If either name is in the ignore list, return true

	si = ign.find( string(name1) );
	if( si != ign.end() )
		return true;

	si = ign.find( string(name2) );
	if( si != ign.end() )
		return true;

// If there is an include_only list, AND if either name
// is not on it, return true.

	if( !gArgs.include_only.size() )
		return false;

	return !IncludeThis( name1 ) || !IncludeThis( name2 );
}
#endif

/* --------------------------------------------------------------- */
/* LookupByNumber ------------------------------------------------ */
/* --------------------------------------------------------------- */

#if 0
// Lookup by number.
// Very inefficient; could (should) keep tables both ways.
//
static string LookupByNumber( map<string,int> &mapping, int k )
{
	map<string,int>::iterator	it;

	for( it = mapping.begin(); it != mapping.end(); ++it ) {

		if( it->second == k )
			return it->first;
	}

	printf( "Looking for mapped string %d, not there?\n", k );
	exit( 42 );
}
#endif

/* --------------------------------------------------------------- */
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#if 0
static void Bounds(
	double			&xbnd,
	double			&ybnd,
	vector<double>	&X )
{
// For plotting, we'd like to know the extent of each patch.
// So for each transform, create a BBox for enclosed points.

	vector<BBox>	Boxes( gNTr );
	int				nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint	&C = vAllC[i];

		if( C.used ) {
			Boxes[vImg[C.n1].itr].AddPoint( C.p1 );
			Boxes[vImg[C.n2].itr].AddPoint( C.p2 );
		}
	}

// Calculate the global XY-bounds

	double	xmin = BIGD, xmax = -BIGD,
			ymin = BIGD, ymax = -BIGD;
	int		nimg = vImg.size();

	for( int i = 0; i < nimg; ++i ) {

		int	itr = vImg[i].itr;

		if( itr < 0 )
			continue;

		TAffine			T( &X[itr * 6] );
		vector<Point>	pts( 4 );

		pts[0] = Point( Boxes[itr].xmin, Boxes[itr].ymin );
		pts[1] = Point( Boxes[itr].xmax, Boxes[itr].ymin );
		pts[2] = Point( Boxes[itr].xmax, Boxes[itr].ymax );
		pts[3] = Point( Boxes[itr].xmin, Boxes[itr].ymax );

		T.Transform( pts );

		for( int k = 0; k < 4; ++k ) {
			xmin = fmin( xmin, pts[k].x );
			xmax = fmax( xmax, pts[k].x );
			ymin = fmin( ymin, pts[k].y );
			ymax = fmax( ymax, pts[k].y );
		}
	}

// Translate all transforms to put global origin at ~(0,0).
// An integer change makes layer-to-layer alignment easier.

	int	xfl = int(floor( xmin )),
		yfl = int(floor( ymin ));

	xmin -= xfl;
	xmax -= xfl;
	ymin -= yfl;
	ymax -= yfl;

	for( int i = 0; i < nimg; ++i ) {

		int	j = vImg[i].itr;

		if( j >= 0 ) {
			j		*= 6;
			X[j+2]	-= xfl;
			X[j+5]	-= yfl;
		}
	}

// Open GNUPLOT files for debugging

	FILE	*fEven		= FileOpenOrDie( "pf.even", "w" ),
			*fOdd		= FileOpenOrDie( "pf.odd", "w" ),
			*fLabEven	= FileOpenOrDie( "pf.labels.even", "w" ),
			*fLabOdd	= FileOpenOrDie( "pf.labels.odd", "w" );

// Write rects and labels

	for( int i = 0; i < nimg; ++i ) {

		int	itr = vImg[i].itr;

		if( itr < 0 )
			continue;

		TAffine			T( &X[itr * 6] );
		vector<Point>	pts( 4 );
		double			xmid = 0.0, ymid = 0.0;

		pts[0] = Point( Boxes[itr].xmin, Boxes[itr].ymin );
		pts[1] = Point( Boxes[itr].xmax, Boxes[itr].ymin );
		pts[2] = Point( Boxes[itr].xmax, Boxes[itr].ymax );
		pts[3] = Point( Boxes[itr].xmin, Boxes[itr].ymax );

		T.Transform( pts );

		for( int k = 0; k < 4; ++k ) {
			xmid += pts[k].x;
			ymid += pts[k].y;
		}

		xmid /= 4.0;
		ymid /= 4.0;

		// select even or odd reporting

		FILE	*f, *flab;
		int		color;

		if( vImg[i].z & 1 ) {
			f		= fOdd;
			flab	= fLabOdd;
			color	= 1;
		}
		else {
			f		= fEven;
			flab	= fLabEven;
			color	= 2;
		}

		// transformed rect

		for( int k = 0; k < 5; ++k )
			fprintf( f, "%f %f\n", pts[k%4].x, pts[k%4].y );

		fprintf( f, "\n" );

		// label

		const char *rgn = strrchr( vImg[i].name.c_str(), ':' ) + 1;

		fprintf( flab, "set label \"%d:%d.%s \" at %f,%f tc lt %d\n",
		vImg[i].z, vImg[i].id, rgn, xmid, ymid, color );
	}

// Close files

	fclose( fLabOdd );
	fclose( fLabEven );
	fclose( fOdd );
	fclose( fEven );

// Report

	printf( "Bounds of global image are x=[%f %f] y=[%f %f].\n",
	xmin, xmax, ymin, ymax );

	fprintf( FOUT, "BBOX %f %f %f %f\n", xmin, ymin, xmax, ymax );

	xbnd = xmax;
	ybnd = ymax;
}
#endif

/* --------------------------------------------------------------- */
/* PrintMagnitude ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void PrintMagnitude( const vector<double> &X, int nvars )
{
	int	k = X.size() - nvars;

	if( k >= 0 ) {

		double	mag	= sqrt( X[k]*X[k] + X[k+1]*X[k+1] );

		printf( "Final magnitude is %g\n", mag );
	}

	fflush( stdout );
}

/* --------------------------------------------------------------- */
/* SolveWithMontageSqr ------------------------------------------- */
/* --------------------------------------------------------------- */

// Add some constraints so the left edge of the array
// needs to be the same length as the right edge, and
// repeat for top and bottom. Of course, this assumes
// a symmetric montage.
//
#if 1
static void SolveWithMontageSqr(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
	const int	MAXINT = 0x7FFFFFFF;

	int	nr = vRgn.size();

/* ------------ */
/* Which layer? */
/* ------------ */

	if( gArgs.ref_layer < 0 ) {

		gArgs.ref_layer = MAXINT;

		for( int i = 0; i < nr; ++i ) {

			if( vRgn[i].itr >= 0 ) {

				gArgs.ref_layer =
					min( gArgs.ref_layer, vRgn[i].z );
			}
		}

		printf(
		"\nNo reference layer specified,"
		" using lowest layer %d\n", gArgs.ref_layer );
	}

/* ---------------------------------------- */
/* Search for greatest outward translations */
/* ---------------------------------------- */

	double	vNW = -MAXINT,
			vNE = -MAXINT,
			vSE = -MAXINT,
			vSW = -MAXINT;
	int		iNW, iNE, iSE, iSW,
			jNW, jNE, jSE, jSW;

	for( int i = 0; i < nr; ++i ) {

		if( vRgn[i].z != gArgs.ref_layer )
			continue;

		if( vRgn[i].itr < 0 )
			continue;

		double	v;
		int		j = vRgn[i].itr * 6;

		if( (v =  X[j+2] + X[j+5]) > vNE )
			{iNE = i; jNE = j; vNE = v;}

		if( (v =  X[j+2] - X[j+5]) > vSE )
			{iSE = i; jSE = j; vSE = v;}

		if( (v = -X[j+2] + X[j+5]) > vNW )
			{iNW = i; jNW = j; vNW = v;}

		if( (v = -X[j+2] - X[j+5]) > vSW )
			{iSW = i; jSW = j; vSW = v;}
	}

	printf(
	"Corner tiles are:"
	" se %d (%f %f),"
	" ne %d (%f %f),"
	" nw %d (%f %f),"
	" sw %d (%f %f).\n",
	vRgn[iSE].id, X[jSE+2], X[jSE+5],
	vRgn[iNE].id, X[jNE+2], X[jNE+5],
	vRgn[iNW].id, X[jNW+2], X[jNW+5],
	vRgn[iSW].id, X[jSW+2], X[jSW+5] );

/* ------------------------------------------- */
/* Use these corner tiles to impose squareness */
/* ------------------------------------------- */

	double	stiff = gArgs.square_strength;

// Top = bottom (DX = DX)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+2,  jSW+2,  jNE+2, jNW+2};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

// Left = right (DY = DY)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+5,  jSW+5,  jNE+5, jNW+5};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

/* --------------- */
/* Update solution */
/* --------------- */

	WriteSolveRead( X, LHS, RHS, "A-Square", 1, false );
	PrintMagnitude( X, 6 );
}
#endif

#if 0
// Experiment to simply hardcode which tiles to use as corners.
// In practice, though, looks like this constraint causes montage
// to buckle if there really is a natural warp like a banana shape,
// so not recommended.
//
static void SolveWithMontageSqr(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
	int	nr = vRgn.size();

/* ------------------------ */
/* Assign hand-picked tiles */
/* ------------------------ */

	int	jNW, jNE, jSE, jSW, nass = 0;

	for( int i = 0; i < nr && nass < 4; ++i ) {

		if( vRgn[i].itr < 0 )
			continue;

		int id = vRgn[i].id;

		if( id == 19001000 ) {
			jNW = i;
			++nass;
		}
		else if( id == 19069000 ) {
			jNE = i;
			++nass;
		}
		else if( id == 19001149 ) {
			jSW = i;
			++nass;
		}
		else if( id == 19069149 ) {
			jSE = i;
			++nass;
		}
	}

	if( nass != 4 ) {
		printf( "***   ***   *** Missing squareness corner.\n" );
		return;
	}

/* ------------------------------------------- */
/* Use these corner tiles to impose squareness */
/* ------------------------------------------- */

	double	stiff = gArgs.square_strength;

// Top = bottom (DX = DX)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+2,  jSW+2,  jNE+2, jNW+2};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

// Left = right (DY = DY)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+5,  jSW+5,  jNE+5, jNW+5};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

/* --------------- */
/* Update solution */
/* --------------- */

	WriteSolveRead( X, LHS, RHS, "A-Square", 1, false );
	PrintMagnitude( X, 6 );
}
#endif


/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* Bill's code --------------------------------------------------- */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* Bill's experimental constraints ------------------------------- */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* --------------------------------------------------------------- */
/* AffineFromTrans ----------------------------------------------- */
/* --------------------------------------------------------------- */

// This method for solving montages omits calling SetPointPairs
// to assert A1(p1) - A2(p2). Rather, we only set approximations
// A(pi) = T(pj). This isn't nearly as accurate as the current
// favorite AffineFromTransWt, which combines SetPointPairs with
// downweighted scaffold constraints A(pi) = T(pj).
//
// Rather than give up on this method, I also tried to improve
// results by a series of scaffold refinements. That is, first
// get T's, use those to get A's, use those to get better A's
// do that again. The result of that slowly converges toward
// AffineFromTransWt but a great computaion cost. These iterative
// experiment snippets follow.
//
void MAffine::AffineFromTrans( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;

	M.SetModelParams( gW, gH, -1, -1, -1, -1,
		nproc, -1, NULL, NULL, zs );
	M.SolveSystem( T, nTr );

// SetPointPairs: A(pi) = T(pj)

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * scaf_strength / sc,
					y1 = C.p1.y * scaf_strength / sc,
					x2 = (C.p2.x + T[k  ]) * scaf_strength / sc,
					y2 = (C.p2.y + T[k+1]) * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * scaf_strength / sc,
					y1 = C.p2.y * scaf_strength / sc,
					x2 = (C.p1.x + T[k  ]) * scaf_strength / sc,
					y2 = (C.p1.y + T[k+1]) * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-FrmT", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* AffineFromAffine ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Use AffineFromTrans as starting point to get first-pass A's
// now get better A's from that scaffold.
//
void MAffine::AffineFromAffine( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Affines A

	vector<double>	A;
	AffineFromTrans( A, nTr );

// SetPointPairs: A(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x * fz / sc;
			y2 = g2.y * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x * fz / sc;
			y2 = g2.y * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-FrmA", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* AffineFromAffine2 --------------------------------------------- */
/* --------------------------------------------------------------- */

// Like AffineFromAffine, try another pass of scaffold refinement.
//
void MAffine::AffineFromAffine2( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Affines A

	vector<double>	A;
	AffineFromAffine( A, nTr );

// SetPointPairs: A(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x * fz / sc;
			y2 = g2.y * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x * fz / sc;
			y2 = g2.y * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-FrmA2", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyFromTrans ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to use approximations H(pi) = T(pj). In practice
// translations lack the accuracy of affines and are always an
// inferior starting point for homographies.
//
void MHmgphy::HmgphyFromTrans( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;

	M.SetModelParams( gW, gH, -1, -1, -1, -1,
		nproc, -1, NULL, NULL, zs );
	M.SolveSystem( T, nTr );

// SetPointPairs: H(pi) = T(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2 = (C.p2.x + T[k  ]) / sc,
					y2 = (C.p2.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2 = (C.p1.x + T[k  ]) / sc,
					y2 = (C.p1.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Solve

	WriteSolveRead( X, LHS, RHS, "H-FrmT", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* AveHTerms ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::AveHTerms(
	double					g[4],
	double					h[4],
	const vector<double>	&X )
{
	int		nr		= vRgn.size(),
			nt[4]	= {0,0,0,0};

	for( int i = 0; i < 4; ++i ) {
		g[i] = 0;
		h[i] = 0;
	}

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		const Til2Img	*m;
		RGN::GetMeta( &m, NULL, vRgn[i], vRgn[i] );

		g[m->cam] += X[itr+6];
		h[m->cam] += X[itr+7];
		++nt[p->cam];
	}

	for( int i = 0; i < 4; ++i ) {

		if( nt[i] ) {
			g[i] /= nt[i];
			h[i] /= nt[i];
		}

		printf( "Cam,G,H: %d\t%.12g\t%.12g\n", i, g[i], h[i] );
	}

	IDBTil2ImgClear();
}

/* --------------------------------------------------------------- */
/* MedHTerms ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::MedHTerms(
	double					g[4],
	double					h[4],
	const vector<double>	&X )
{
	int						nr = vRgn.size();
	vector<vector<double> >	G( 4 ), H( 4 );

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		const Til2Img	*m;
		RGN::GetMeta( &m, NULL, vRgn[i], vRgn[i] );

		G[m->cam].push_back( X[itr+6] );
		H[m->cam].push_back( X[itr+7] );
	}

	for( int i = 0; i < 4; ++i ) {

		g[i] = MedianVal( G[i] );
		h[i] = MedianVal( H[i] );

		printf( "Cam,G,H: %d\t%.12g\t%.12g\n", i, g[i], h[i] );
	}

	IDBTil2ImgClear();
}

/* --------------------------------------------------------------- */
/* ForceHTerms --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::ForceHTerms(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	const double	g[4],
	const double	h[4] )
{
	double	wt	= 0.1;
	int		nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		int	j = vRgn[i].itr;

		if( j < 0 )
			continue;

		j *= NX;

		const Til2Img	*m;
		RGN::GetMeta( &m, NULL, vRgn[i], vRgn[i] );

		double	v[1]	= { wt };
		int		i1[1]	= { j+6 },
				i2[1]	= { j+7 };

		AddConstraint( LHS, RHS, 1, i1, v, wt*g[m->cam] );
		AddConstraint( LHS, RHS, 1, i2, v, wt*h[m->cam] );
	}

	IDBTil2ImgClear();
}

/* --------------------------------------------------------------- */
/* HmgphyFromHmgphy_ForceGH -------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to get homographies and either average or median
// their (g,h) components (per camera). Then add constraints to
// push all g,h to these values.
//
// Note median and average yield nearly same results on full
// size montages (not so many outliers). Don't see improvement
// in error distribution from this, slightly worse, but spread
// in g,h params is only one tenth of standard spread.
//
void MHmgphy::HmgphyFromHmgphy_ForceGH( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Homographies A

	vector<double>	A;
	double			g[4], h[4];

	HmgphyFromAffine( A, nTr );

	AveHTerms( g, h, A );
	//MedHTerms( g, h, A );
	ForceHTerms( LHS, RHS, g, h );

// SetPointPairs: H(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = A(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Solve

	WriteSolveRead( X, LHS, RHS, "H-FrmA", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}


