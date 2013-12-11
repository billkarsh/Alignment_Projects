/* --------------------------------------------------------------- */
/* AontoBOverlap ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return approximated fraction of a on b area overlap.
//
static double AontoBOverlap( const Picture &a, const Picture &b )
{
	TAffine			T = b.inv * a.tr;	// map a->b
	vector<Point>	cnr;
	int				w = a.w, h = a.h;

	Set4Corners( cnr, w, h );
	T.Transform( cnr );

// bounding box.

	double xmin = 1E9, xmax = -1E9;
	double ymin = 1E9, ymax = -1E9;

	for( int k = 0; k < 4; ++k ) {
		xmin = fmin( xmin, cnr[k].x );
		ymin = fmin( ymin, cnr[k].y );
		xmax = fmax( xmax, cnr[k].x );
		ymax = fmax( ymax, cnr[k].y );
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


