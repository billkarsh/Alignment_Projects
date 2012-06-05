

#include	"CThmScan.h"
#include	"Correlation.h"
#include	"Maths.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CThmScan*	ME;






/* --------------------------------------------------------------- */
/* BigEnough ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool BigEnough( int sx, int sy, void *a )
{
	return	sx >= ME->olap1D &&
			sy >= ME->olap1D &&
			(long)sx * sy > (long)a;
}

/* --------------------------------------------------------------- */
/* EnoughPoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool EnoughPoints( int count1, int count2, void *a )
{
	return count1 > (long)a && count2 > (long)a;
}

/* --------------------------------------------------------------- */
/* NewXFromParabola ---------------------------------------------- */
/* --------------------------------------------------------------- */

static double NewXFromParabola(
	double	x1,
	double	d,
	double	y0,
	double	y1,
	double	y2 )
{
	double D = d * (y0 - y2) / (2 * (y0 + y2 - y1 - y1));

	if( fabs( D ) <= 2 * d )
		x1 += D;

	return x1;
}

/* --------------------------------------------------------------- */
/* RotatePoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::RotatePoints(
	vector<Point>	&pts,
	TForm			&T,
	double			rads )
{
	TForm	A, temp;

	A.NUSetRot( rads );
	MultiplyTrans( temp, Tuser, Tptwk );
	MultiplyTrans( T, A, temp );
	T.Apply_R_Part( pts );
}

/* --------------------------------------------------------------- */
/* PTWApply1 ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Create product of Ttry with Tptwk and return new corr.
// Make change to Tptwk permanent if keep = true.
//
double CThmScan::PTWApply1(
	const TForm	&Ttry,
	double		deg,
	ThmRec		&thm,
	bool		keep )
{
	CorRec	C;
	TForm	Tback = Tptwk;

	MultiplyTrans( Tptwk, Ttry, Tback );
	RFromAngle( C, deg, thm );

	if( !keep )
		Tptwk = Tback;

	return C.R;
}

/* --------------------------------------------------------------- */
/* PTWSweep ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// For near unity transform chosen by sel, perform a magnitude
// sweep about unity in +/- nsteps of size astep. Return the
// best NU-magnitude and its corresponding R.
//
double CThmScan::PTWSweep(
	double	&rbest,
	int		sel,
	int		nstep,
	double	astep,
	double	deg,
	ThmRec	&thm )
{
	double	abase = 0.0;
	int		ibest = 0;

	rbest = -1.0;

	fprintf( flog, "PTWSweep %2d:", sel );

	if( sel >= tfnuScl && sel <= tfnuYScl )
		abase = 1.0;

	for( int i = -nstep; i <= nstep; ++i ) {

		double	R;
		TForm	T;

		T.NUSelect( sel, abase + i * astep );
		R = PTWApply1( T, deg, thm, false );
		fprintf( flog, " %5.3f", R );

		if( R > rbest ) {
			rbest = R;
			ibest = i;
		}
	}

	fprintf( flog, "\n" );

	return abase + ibest * astep;
}

/* --------------------------------------------------------------- */
/* PTWInterp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Use NewXFromParabola to interpolate a better NU-magnitude
// where the center guess is (x1,y1) from the sweep and the
// side points are +/1 d from that.
//
// Return xnew.
//
double CThmScan::PTWInterp(
	double	&ynew,
	int		sel,
	double	x1,
	double	y1,
	double	d,
	double	deg,
	ThmRec	&thm )
{
	double	y0, y2, xnew;
	TForm	T;

	T.NUSelect( sel, x1 - d );
	y0 = PTWApply1( T, deg, thm, false );

	T.NUSelect( sel, x1 + d );
	y2 = PTWApply1( T, deg, thm, false );

	xnew = NewXFromParabola( x1, d, y0, y1, y2 );

	T.NUSelect( sel, xnew );
	ynew = PTWApply1( T, deg, thm, false );

	if( ynew < y1 ) {
		xnew = x1;
		ynew = y1;
	}

	fprintf( flog, "PTWIntrp %2d: %5.3f @ %.3f\n", sel, ynew, xnew );

	return xnew;
}

/* --------------------------------------------------------------- */
/* PeakHunt ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Bracket search for peak.
//
double CThmScan::PeakHunt( CorRec &best, double hlfwid, ThmRec &thm )
{
	CorRec	C,
			B0	= best;
	double	L	= best.A - hlfwid,
			R	= best.A + hlfwid,
			M, lr, rr;
	int		k	= 1;

	clock_t	t0 = StartTiming();

	RFromAngle( C, L, thm );
	lr = C.R;

	RFromAngle( C, R, thm );
	rr = C.R;

	RFromAngle( best, M = (L+R)/2.0, thm );

	while( R - L > 0.0001 ) {

		double	a;

		++k;

		// move left up
		RFromAngle( C, a = (L+M)/2.0, thm );

		if( C.R >= best.R ) {
			rr		= best.R;
			R		= M;
			best	= C;
			M		= a;
		}
		else {
			lr		= C.R;
			L		= a;
		}

		// move right back
		RFromAngle( C, a = (M+R)/2.0, thm );

		if( C.R >= best.R ) {
			rr		= best.R;
			L		= M;
			best	= C;
			M		= a;
		}
		else {
			rr		= C.R;
			R		= a;
		}
	}

	if( B0.R > best.R )
		best = B0;

	fprintf( flog,
	"PeakHunt: Best: K=%d, R=%.3f, A=%.3f, X=%.3f, Y=%.3f\n",
	k, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );

	return best.R;
}

/* --------------------------------------------------------------- */
/* CThmScan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

CThmScan::CThmScan()
{
	ME			= this;
	flog		= stdout;
	rthresh		= 0.30;
	nbmaxht		= 0.75;
	err			= errOK;
	swpSimple	= false;
	swpPretweak	= true;
	useCorrR	= false;
	olap1D		= 0;
	Ox			= 0;
	Oy			= 0;
	Or			= -1;
}

/* --------------------------------------------------------------- */
/* CThmScan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CThmScan::Initialize( FILE* flog, CorRec &best )
{
	SetFlog( flog );

	fprintf( flog, "\n---- Thumbnail matching ----\n" );

	best.X	= -999.0;
	best.Y	= -999.0;
	best.A	= -999.0;
	best.R	= -999.0;
}

/* --------------------------------------------------------------- */
/* SetTuser ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CThmScan::SetTuser(
	double	scl,
	double	xscl,
	double	yscl,
	double	xskw,
	double	yskw )
{
	TForm	X, Y, S, temp;

	X.NUSetXSkw( xskw );
	Y.NUSetYSkw( yskw );

	S.t[0] = xscl*scl; S.t[1] = 0;        S.t[2] = 0;
	S.t[3] = 0;        S.t[4] = yscl*scl; S.t[5] = 0;

	MultiplyTrans( temp, Y, X );
	MultiplyTrans( Tuser, S, temp );
}

/* --------------------------------------------------------------- */
/* Pretweaks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Several near-unity transforms are applied to A to boost initial
// correlation at the given center angle.
//
// Return true if any changes made.
//
bool CThmScan::Pretweaks( double bestR, double deg, ThmRec &thm )
{
	fprintf( flog, "Pretweaks start, best R=%.3f\n", bestR );

	clock_t	t0 = StartTiming();
	bool	anychange = false;

	Tptwk.NUSetOne();

// We will compose a Tptwk by multiplying NU transforms. We examine
// 5 transform types vsel = {Scl, XScl, YScl, XSkw, YSkw}. Each type
// will enter the product at most once, which is tracked with vused.

	vector<int>	vsel( 5 );
	vector<int>	vused( 5, 0 );

	for( int i = tfnuScl; i <= tfnuYSkw; ++i )
		vsel[i] = i;

// To decide which transform type to use, we will do a magnitude
// sweep with each and see which type got the highest peak R. We
// then use interpolation to improve the winner. If using that
// type gives a better R than before, it goes into the product.
// Repeat with other unused types.

	for( int itype = 0; itype < 5; ++itype ) {

		// Do the sweeps

		vector<double>	vrbest( 5, 0.0 );
		vector<double>	vabest( 5 );

		for( int i = 0; i < 5; ++i ) {

			if( vused[i] )
				continue;

			vabest[i] = PTWSweep( vrbest[i], vsel[i],
							5, .02, deg, thm );
		}

		// Find the best sweep

		double	rbest	= -1;
		int		selbest = -1;

		for( int i = 0; i < 5; ++i ) {

			if( vused[i] )
				continue;

			if( vrbest[i] > rbest ) {
				rbest	= vrbest[i];
				selbest	= i;
			}
		}

		if( selbest == -1 )
			break;

		// Improve candidate with interpolator

		double	a = PTWInterp( rbest, vsel[selbest],
						vabest[selbest], vrbest[selbest],
						.02, deg, thm );

		// Is it better than before?

		if( rbest > bestR ) {

			TForm	T;

			fprintf( flog, "PTWUsing %2d\n", vsel[selbest] );
			T.NUSelect( vsel[selbest], a );
			MultiplyTrans( Tptwk, T, TForm( Tptwk ) );
			vused[selbest]	= 1;
			bestR			= rbest;
			anychange		= true;
		}
		else
			break;
	}

	StopTiming( stdout, "Pretweaks", t0 );

	fprintf( flog, "Approx: Pretweak " );
	Tptwk.PrintTransform( flog );

	return anychange;
}

/* --------------------------------------------------------------- */
/* RFromAngle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::RFromAngle( CorRec &C, double deg, ThmRec &thm )
{
	vector<Point>	pts = thm.ap;

	C.A = deg;

	RotatePoints( pts, C.T, deg * PI/180.0 );

	C.R = (useCorrR ? CorrImagesR : CorrImages)(
		flog, false, C.X, C.Y,
		pts, thm.av, thm.bp, thm.bv,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		0.0, nbmaxht, Ox, Oy, Or, thm.ftc );
}

/* --------------------------------------------------------------- */
/* DebugAngs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::DebugAngs(
	int		alyr,
	int		atil,
	int		blyr,
	int		btil,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	char	file[256];

	sprintf( file, "angs_%d_%d_@_%d_%d.log", alyr, atil, blyr, btil );
	FILE	*f = fopen( file, "w" );

	fprintf( f, "Deg\tR\tX\tY\n" );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		RFromAngle( C, a, thm );

		fprintf( f, "%.3f\t%.4f\t%.3f\t%.3f\n",
			a, C.R, C.X, C.Y );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* RecordAngle --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RecordAngle(
	FILE			*f,
	const char		*label,
	const CorRec	&C )
{
	fprintf( f,
	"%s: R=%6.3f, A=%8.3f, X=%8.2f, Y=%8.2f\n",
	label, C.R, C.A, C.X, C.Y );
}

/* --------------------------------------------------------------- */
/* AngleScanMaxR ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Original version selecting best angle based on maximal R.
//
double CThmScan::AngleScanMaxR(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	fprintf( flog,
	"AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f\n",
	center, hlfwid, step );

	clock_t	t0 = StartTiming();

	RFromAngle( best, center, thm );
	RecordAngle( flog, "Center", best );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		RFromAngle( C, a, thm );
		RecordAngle( flog, "  Scan", C );

		if( C.R > best.R )
			best = C;
	}

	RecordAngle( flog, "  Best", best );

	StopTiming( stdout, "AngleScan", t0 );

	return best.R;
}

/* --------------------------------------------------------------- */
/* AngleScanConstXY ---------------------------------------------- */
/* --------------------------------------------------------------- */

// New version selecting best angle based on maximal R, but only if
// coordinates XY vary smoothly in a neighborhood about that angle.
// Smoothness measured as linear corr. coeff. from fitting X vs A.

static const vector<CorRec>	*_vC;


static bool Sort_vC_dec( int a, int b )
{
	return (*_vC)[a].R > (*_vC)[b].R;
}


double CThmScan::AngleScanConstXY(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	fprintf( flog,
	"AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f\n",
	center, hlfwid, step );

//RFromAngle( best, center, thm );
//return;

	clock_t	t0 = StartTiming();

	best.X	= 0.0;
	best.Y	= 0.0;
	best.A	= 0.0;
	best.R	= 0.0;

// Sweep and collect

	vector<CorRec>	vC;

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		RFromAngle( C, a, thm );
		RecordAngle( flog, "Scan", C );
		vC.push_back( C );
	}

// Make indices sorted by decreasing R

	int	nC = vC.size();

	vector<int>	order( nC );

	for( int i = 0; i < nC; ++i )
		order[i] = i;

	_vC = &vC;

	sort( order.begin(), order.end(), Sort_vC_dec );

// Evaluate sweep
// Search through decreasing R.
// Take first result for which there are also +- m data points, and
// the X and Y coords vary smoothly over that range (r > rthresh).
// Also note that coords are rounded to nearest pixel, and that
// lincor values that are NAN or INF are generally ok: they arise
// from a zero variance in the coords, meaning that the point is
// stationary, and that's a good thing because it's not noise.

	const double	rthresh = 0.7;
	const			int m = 3;
	const			int M = 2*m + 1;

	for( int i = 0; i < nC; ++i ) {

		int	ic = order[i];

		if( ic < m )
			continue;

		if( nC-1 - ic < m )
			continue;

		vector<double>	A( M ), X( M ), Y( M );
		double			lincor;
		int				n = 0;

		for( int j = ic - m; j <= ic + m; ++j ) {

			A[n] = vC[j].A;
			X[n] = ROUND( vC[j].X );
			Y[n] = ROUND( vC[j].Y );
			++n;
		}

		LineFit( NULL, NULL, &lincor, &A[0], &X[0], 0, n );
		fprintf( flog, "LCOR: A=%8.3f, RX=%6.3f\n", A[m], lincor );

		if( isfinite( lincor ) && fabs( lincor ) < rthresh )
			continue;

		LineFit( NULL, NULL, &lincor, &A[0], &Y[0], 0, n );
		fprintf( flog, "LCOR: A=%8.3f, RY=%6.3f\n", A[m], lincor );

		if( isfinite( lincor ) && fabs( lincor ) < rthresh )
			continue;

		best = vC[ic];
		break;
	}

// Report

	RecordAngle( flog, "Best", best );

	StopTiming( stdout, "AngleScan", t0 );

	return best.R;
}

/* --------------------------------------------------------------- */
/* AngleScanSel -------------------------------------------------- */
/* --------------------------------------------------------------- */

double CThmScan::AngleScanSel(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	if( swpSimple )
		return AngleScanMaxR( best, center, hlfwid, step, thm );
	else
		return AngleScanConstXY( best, center, hlfwid, step, thm );
}

/* --------------------------------------------------------------- */
/* AngleScanWithTweaks ------------------------------------------- */
/* --------------------------------------------------------------- */

double CThmScan::AngleScanWithTweaks(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	if( AngleScanSel( best, center, hlfwid, step, thm ) < rthresh ) {

		if( swpPretweak ) {

			if( Pretweaks( best.R,
					(best.R > 0.0 ? best.A : center), thm ) ) {

				AngleScanSel( best, center, hlfwid, step, thm );
			}
		}
	}

	return best.R;
}

/* --------------------------------------------------------------- */
/* UsePriorAngles ------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CThmScan::UsePriorAngles(
	CorRec	&best,
	int		nprior,
	double	ang0,
	double	hfangpr,
	ThmRec	&thm )
{
	fprintf( flog, "Approx: Using prior angles n=%d, med=%f\n",
	nprior, ang0 );

	if( AngleScanWithTweaks( best, ang0, hfangpr, 0.1, thm )
			< rthresh ||
		PeakHunt( best, 0.3, thm )
			< rthresh ) {

		fprintf( flog,
		"FAIL: Approx: Prior angles R=%g below thresh=%g\n",
		best.R, rthresh );

		err = errLowRPrior;
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* DenovoBestAngle ----------------------------------------------- */
/* --------------------------------------------------------------- */

bool CThmScan::DenovoBestAngle(
	CorRec	&best,
	double	ang0,
	double	hfangdn,
	ThmRec	&thm )
{
	if( AngleScanWithTweaks( best, ang0, hfangdn, 0.5, thm )
			< rthresh ||
		AngleScanSel( best, best.A, 1.0, 0.1, thm )
			< rthresh ||
		PeakHunt( best, 0.3, thm )
			< rthresh ) {

		fprintf( flog,
		"FAIL: Approx: Denovo R=%g below thresh=%g\n",
		best.R, rthresh );

		err = errLowRDenov;
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* PostTweaks ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::PostTweaks( CorRec &best, ThmRec &thm )
{
	SetUseCorrR( true );
	SetDisc( (int)best.X, (int)best.Y, 40 );
	Pretweaks( best.R, best.A, thm );
	RFromAngle( best, best.A, thm );
}

/* --------------------------------------------------------------- */
/* FinishAtFullRes ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Get final offsets at full resolution.
//
void CThmScan::FinishAtFullRes( CorRec &best, ThmRec &thm )
{
	CorRec	b0 = best;
	int		ok;

	fprintf( flog,
	"Approx: LowRes  R=%.3f, X=%.3f, Y=%.3f\n",
	best.R, best.X, best.Y );

	clock_t	t0 = StartTiming();

	SetUseCorrR( true );
	SetDisc( (int)best.X, (int)best.Y, 40 );
	RFromAngle( best, best.A, thm );

	ok = (fabs( best.X - b0.X ) <= 20)
	  && (fabs( best.Y - b0.Y ) <= 20);

	fprintf( flog,
	"Approx: FullRes R=%.3f, X=%.3f, Y=%.3f, use=%c.\n",
	best.R, best.X, best.Y, (ok ? 'Y' : 'N') );

// Always report the thumb-R instead of the fullres value because
// fullres values are scaled differently. We prefer all values in
// the file to be comparable.

	if( ok )
		best.R = b0.R;
	else
		best = b0;

	StopTiming( stdout, "FinishAtFullRes", t0 );
}


