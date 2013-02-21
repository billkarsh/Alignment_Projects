/* --------------------------------------------------------------- */
/* AontoBOverlap ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return approximated fraction of a on b area overlap.
//
static double AontoBOverlap( const Picture &a, const Picture &b )
{
	TAffine			T = b.inv * a.tr;	// map a->b
	vector<Point>	corners( 4 );
	int				w = a.w, h = a.h;

	corners[0] = Point( 0.0, 0.0 );
	corners[1] = Point( w-1, 0.0 );
	corners[2] = Point( 0.0, h-1 );
	corners[3] = Point( w-1, h-1 );

// bounding box.

	double xmin = 1E9, xmax = -1E9;
	double ymin = 1E9, ymax = -1E9;

	for( int k = 0; k < 4; ++k ) {

		T.Transform( corners[k] );

		xmin = fmin( xmin, corners[k].x );
		ymin = fmin( ymin, corners[k].y );
		xmax = fmax( xmax, corners[k].x );
		ymax = fmax( ymax, corners[k].y );
	}

// any overlap possibility?

	if( xmin > w-1 || xmax < 0 || ymin > h-1 || ymax < 0 )
		return 0.0;

// approximate area using sampling of random b-points.

	const int count = 4000;

	double	wf	= double(w-1) / RAND_MAX;
	double	hf	= double(h-1) / RAND_MAX;
	int		in	= 0;

	for( int i = 0; i < count; ++i ) {

		Point p( wf*rand(), hf*rand() );

		T.Transform( p );

		if( p.x >= 0 && p.x < w && p.y >= 0.0 && p.y < h )
			++in;
	}

	fprintf( flog, "----AontoBOverlap fraction %f\n",
		double(in)/count );

	return double(in)/count;
}


