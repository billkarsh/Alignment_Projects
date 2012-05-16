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

		TForm			T( &X[itr * 6] );
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

		TForm			T( &X[itr * 6] );
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


