/* --------------------------------------------------------------- */
/* PeakHunt ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Repeated parabola fits for peak on a narrowing angle range.
//
// *** Retired in favor of current bracket search.
//
static double PeakHunt(
	CorRec	&best,
	double	hlfwid,
	ThmRec	&thm,
	FILE*	flog )
{
	CorRec	B0	= best;
	double	L	= best.A - hlfwid,
			R	= best.A + hlfwid;
	int		k	= 0;

	clock_t	t0 = StartTiming();

	while( R - L > 0.0001 ) {

		CorRec	C;
		double	x1, y0, y2;

		++k;

		// left
		RFromAngle( C, L, thm, flog );
		y0 = C.R;

		// right
		RFromAngle( C, R, thm, flog );
		y2 = C.R;

		// middle
		RFromAngle( C, x1 = (L+R)/2.0, thm, flog );

		// estimate peak position
		C.A = NewXFromParabola( x1, (R-L)/2.0, y0, C.R, y2 );

		// sanity check
		if( C.A <= L || C.A >= R )
			break;

		// compete estimate against best so far
		RFromAngle( C, C.A, thm, flog );

		if( C.R > best.R )
			best = C;

		//fprintf( flog,
		//"*** [%f,%f] [%f,%f] [%f,%f] <newx=%f>.\n",
		//L,y0, best.A,best.R, R,y2, x1 );

		x1	= fmin( best.A - L, R - best.A ) / 4.0;
		L	= best.A - x1;
		R	= best.A + x1;
	}

	if( B0.R > best.R )
		best = B0;

	fprintf( flog,
	"PeakHunt: Best: K=%d, R=%.3f, A=%.3f, X=%.3f, Y=%.3f\n",
	k, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );

	return best.R;
}


