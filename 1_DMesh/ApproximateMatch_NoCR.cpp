

#include	"CGBL_dmesh.h"
#include	"ApproximateMatch.h"

#include	"Disk.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Timer.h"
#include	"Debug.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum thmerrs {
	errOK			= 0,
	errLowRDenov	= 1,
	errLowRPrior	= 2
};

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	vector<double>	av, bv;
	vector<Point>	ap, bp;
	Point			aO, bO;			// subimage ref. points
	int				ow, oh;			// common sub-img dims
} OlapRec;

typedef struct {
	vector<double>	av, bv;
	vector<Point>	ap, bp;
	vector<CD>		ftc;			// fourier transform cache
	long			reqArea;
	int				scl;
} ThmRec;

typedef struct {
	TForm			T;
	double			X, Y,
					A, R;
} CorRec;

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static int		gErr	= errOK;
static double	ang0	= 0.0;
static TForm	Tptwk;






/* --------------------------------------------------------------- */
/* Orig ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Orig( vector<TForm> &G )
{
	TForm	T;

	AToBTrans( T, GBL.A.t2i.T, GBL.B.t2i.T );

	G.push_back( T );

	return true;
}

/* --------------------------------------------------------------- */
/* FromLog ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool FromLog( vector<TForm> &G, FILE* flog )
{
	ThmPair	tpr;
	int		ok;

	ok = ReadThmPair( tpr, GBL.A.layer, GBL.A.tile, 1,
			GBL.B.layer, GBL.B.tile, 1, flog )
		&& !tpr.err;

	if( ok )
		G.push_back( tpr.T );

	return ok;
}

/* --------------------------------------------------------------- */
/* SetStartingAngle ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Estimate starting angle (ang0) and use it to adjust OLAP2D_XL.
//
// Return {0=use denovo method, n=count of prior angles}.
//
static int SetStartingAngle( FILE* flog )
{
	vector<ThmPair>	tpr;
	int				ntpr, nprior = 0;

// Handle user override

	if( GBL.arg.CTR != 999.0 ) {

		ang0	= GBL.arg.CTR;
		nprior	= 1;
		goto adjust_olap;
	}

// Try to get prior angles

	if( ReadAllThmPair( tpr, GBL.A.layer, GBL.B.layer, flog )
		&& (ntpr = tpr.size()) ) {

		vector<double>	A;

		for( int i = 0; i < ntpr; ++i ) {

			if( !tpr[i].err )
				A.push_back( tpr[i].A );
		}

		if( (nprior = A.size()) >= 4 )
			ang0 = MedianVal( A );
		else
			nprior = 0;
	}

// Otherwise estimate from initial tforms

	if( !nprior ) {

		TForm	atob;

		AToBTrans( atob, GBL.A.t2i.T, GBL.B.t2i.T );

		ang0 = 180.0/PI * RadiansFromAffine( atob );
	}

// Force tiny ang0 to zero

	if( fabs( ang0 ) < 0.001 )
		ang0 = 0.0;

// Adjust OLAP2D_XL

adjust_olap:
	if( GBL.A.layer != GBL.B.layer ) {

		double	a = ang0 * PI/180.0,
				c = cos( a ),
				s = sin( a );

		GBL.ctx.OLAP2D =
			int(GBL.ctx.OLAP2D / fmax( c*c, s*s ));
	}

	return nprior;
}

/* --------------------------------------------------------------- */
/* WholeImage ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set fields of olp appropriate for whole imaging searching.
//
// We do that for all cross-layer work, or if t2i data suggest
// overlaps are unset or inaccurate.
//
// Always return true from this case.
//
static bool WholeImage(
	OlapRec			&olp,
	const PixPair	&px,
	FILE*			flog )
{
	int	w = px.ws,
		h = px.hs;

	fprintf( flog,
	"Subimage: Using whole images, pix=%d.\n", w * h );

// Values

	olp.av = *px.avs_aln;
	olp.bv = *px.bvs_aln;

// Points

	MakeZeroBasedPoints( olp.ap, w, h );
	olp.bp = olp.ap;

// Dimensions

	olp.aO	= Point( 0, 0 );
	olp.bO	= olp.aO;
	olp.ow	= w;
	olp.oh	= h;

	return true;
}

/* --------------------------------------------------------------- */
/* SelectSubimage ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Describe the conservative a/b intersection region based upon
// data from TileToImage.txt.
//
static bool SelectSubimage(
	OlapRec			&olp,
	const PixPair	&px,
	FILE*			flog )
{
// Use whole images if cross layer

	if( GBL.A.layer != GBL.B.layer )
		return WholeImage( olp, px, flog );

// Get displacements of a in b-system

	int	dx, dy;

	{
		TForm	atob;
		Point	delta;

		AToBTrans( atob, GBL.A.t2i.T, GBL.B.t2i.T );
		atob.Transform( delta );

		dx = int(delta.x) / px.scl;
		dy = int(delta.y) / px.scl;

		if( GBL.mch.SLOPPY_SL ) {
			dx /= 2;
			dy /= 2;
		}
	}

// Use offsets and image size to determine:
// - {ow, oh} common overlap dimesnions.
// - {ax, ay} intersection corner in a-coordinates.
// - {bx, by} ...and in b-coordinates.

	int	w = px.ws,
		h = px.hs,
		np,
		ax, ay,
		bx, by, min1d;

	if( dx < 0 ) {
		olp.ow = w + dx;
		ax = -dx;
		bx = 0;
	}
	else {
		olp.ow = w - dx;
		ax = 0;
		bx = dx;
	}

	if( dy < 0 ) {
		olp.oh = h + dy;
		ay = -dy;
		by = 0;
	}
	else {
		olp.oh = h - dy;
		ay = 0;
		by = dy;
	}

// Double-check that there was sufficient overlap

	min1d = max( GBL.ctx.OLAP1D, 8 );

	if( olp.ow < min1d || olp.oh < min1d ) {
		fprintf( flog, "Subimage: Overlap looks small.\n" );
		return WholeImage( olp, px, flog );
	}

	np = olp.ow * olp.oh;

	fprintf( flog,
	"Subimage: Using intersection, w=%d, h=%d, pix=%d.\n",
	olp.ow, olp.oh, np );

// Corners and a-centroid

	olp.aO	= Point( ax, ay );
	olp.bO	= Point( bx, by );

// Set points

	MakeZeroBasedPoints( olp.ap, olp.ow, olp.oh );
	olp.bp = olp.ap;

// Set values

	const vector<double>&	av = *px.avs_aln;
	const vector<double>&	bv = *px.bvs_aln;

	olp.av.resize( np );
	olp.bv.resize( np );

	for( int i = 0; i < np; ++i ) {

		int	y = i / olp.ow;
		int	x = i - olp.ow * y;

		olp.av[i] = av[ax + x + w*(ay + y)];
		olp.bv[i] = bv[bx + x + w*(by + y)];
	}

	return true;
}

/* --------------------------------------------------------------- */
/* MakeThumbs ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeThumbs(
	ThmRec			&thm,
	const OlapRec	&olp,
	int				decfactor,
	FILE*			flog )
{
	thm.av		= olp.av;
	thm.bv		= olp.bv;
	thm.ap		= olp.ap;
	thm.bp		= olp.bp;
	thm.ftc.clear();
	thm.reqArea	= GBL.ctx.OLAP2D;
	thm.scl		= decfactor;

	if( decfactor > 1 ) {

		DecimateVector( thm.ap, thm.av, olp.ow, olp.oh, decfactor );
		DecimateVector( thm.bp, thm.bv, olp.ow, olp.oh, decfactor );

		thm.reqArea	/= decfactor * decfactor;

		fprintf( flog,
		"Thumbs: After decimation %d pts, reqArea %d, thmscl %d.\n",
		thm.ap.size(), thm.reqArea, thm.scl );
	}

	Normalize( thm.av );
	Normalize( thm.bv );
}

/* --------------------------------------------------------------- */
/* RotatePoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RotatePoints(
	vector<Point>	&pts,
	TForm			&T,
	const TForm		&Tskew,
	double			theta )
{
	double	c	= cos( theta ) * GBL.ctx.SCALE,
			s	= sin( theta ) * GBL.ctx.SCALE;
	TForm	ao( GBL.ctx.XSCALE * c, -GBL.ctx.YSCALE * s, 0.0,
				GBL.ctx.XSCALE * s,  GBL.ctx.YSCALE * c, 0.0 );
	TForm	T0;

	MultiplyTrans( T0, Tskew, Tptwk );
	MultiplyTrans( T, ao, T0 );
	T.Apply_R_Part( pts );
}

/* --------------------------------------------------------------- */
/* BigEnough ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool BigEnough( int sx, int sy, void *a )
{
	return	sx >= GBL.ctx.OLAP1D &&
			sy >= GBL.ctx.OLAP1D &&
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
/* RFromAngle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RFromAngle(
	CorRec	&C,
	double	a,
	ThmRec	&thm,
	FILE*	flog )
{
	TForm			Tskew( 1.0, 0.0, 0.0, GBL.ctx.SKEW, 1.0, 0.0 );
	vector<Point>	ps = thm.ap;

	C.A = a;

	RotatePoints( ps, C.T, Tskew, a * PI/180.0 );

	C.R = CorrImages(
		flog, false, C.X, C.Y,
		ps, thm.av, thm.bp, thm.bv,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		0.0, GBL.ctx.NBMXHT, thm.ftc );
}

/* --------------------------------------------------------------- */
/* DebugAngs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DebugAngs(
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm )
{
	char	file[256];

	sprintf( file, "angs_%d_%d_@_%d_%d.log",
		GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile );
	FILE	*f = fopen( file, "w" );

	fprintf( f, "Deg\tR\tX\tY\n" );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		RFromAngle( C, a, thm, stdout );

		fprintf( f, "%.3f\t%.4f\t%.3f\t%.3f\n",
			a, C.R, C.X, C.Y );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* AngleScan ----------------------------------------------------- */
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


#if 0

// Original version selecting best angle based on maximal R.

static double AngleScan(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm,
	FILE*	flog )
{
	fprintf( flog,
	"AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f.\n",
	center, hlfwid, step );

	clock_t	t0 = StartTiming();

	RFromAngle( best, center, thm, flog );
	RecordAngle( flog, "Center", best );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		RFromAngle( C, a, thm, flog );
		RecordAngle( flog, "  Scan", C );

		if( C.R > best.R )
			best = C;
	}

	RecordAngle( flog, "  Best", best );

	StopTiming( stdout, "AngleScan", t0 );

	return best.R;
}

#else

// New version selecting best angle based on maximal R, but only if
// coordinates XY vary smoothly in a neighborhood about that angle.
// Smoothness measured as linear corr. coeff. from fitting X vs A.

static const vector<CorRec>	*_vC;


static bool Sort_vC_dec( int a, int b )
{
	return (*_vC)[a].R > (*_vC)[b].R;
}


static double AngleScan(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm,
	FILE*	flog )
{
	fprintf( flog,
	"AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f.\n",
	center, hlfwid, step );

//RFromAngle( best, center, thm, flog );
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

		RFromAngle( C, a, thm, flog );
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

#endif

/* --------------------------------------------------------------- */
/* PTWApply1 ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Create product of Ttry with Tptwk and return new corr.
// Make change to Tptwk permanent if keep = true.
//
static double PTWApply1(
	const TForm	&Ttry,
	double		center,
	ThmRec		&thm,
	FILE*		flog,
	bool		keep )
{
	CorRec	C;
	TForm	Tback = Tptwk;

	MultiplyTrans( Tptwk, Ttry, Tback );
	RFromAngle( C, center, thm, flog );

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
static double PTWSweep(
	double	&rbest,
	int		sel,
	int		nstep,
	double	astep,
	double	center,
	ThmRec	&thm,
	FILE*	flog )
{
	double	abase = 0.0;
	int		ibest;

	rbest = -1.0;

	fprintf( flog, "PTWSweep %2d:", sel );

	if( sel >= tfnuScl && sel <= tfnuYScl )
		abase = 1.0;

	for( int i = -nstep; i <= nstep; ++i ) {

		double	R;
		TForm	T;

		T.NUSelect( sel, abase + i * astep );
		R = PTWApply1( T, center, thm, flog, false );
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
/* PTWInterp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Use NewXFromParabola to interpolate a better NU-magnitude
// where the center guess is (x1,y1) from the sweep and the
// side points are +/1 d from that.
//
static double PTWInterp(
	double	&ynew,
	int		sel,
	double	x1,
	double	y1,
	double	d,
	double	center,
	ThmRec	&thm,
	FILE*	flog )
{
	double	y0, y2, xnew;
	TForm	T;

	T.NUSelect( sel, x1 - d );
	y0 = PTWApply1( T, center, thm, flog, false );

	T.NUSelect( sel, x1 + d );
	y2 = PTWApply1( T, center, thm, flog, false );

	xnew = NewXFromParabola( x1, d, y0, y1, y2 );

	T.NUSelect( sel, xnew );
	ynew = PTWApply1( T, center, thm, flog, false );

	if( ynew < y1 ) {
		xnew = x1;
		ynew = y1;
	}

	fprintf( flog, "PTWIntrp %2d: %5.3f @ %.3f\n", sel, ynew, xnew );

	return xnew;
}

/* --------------------------------------------------------------- */
/* Pretweaks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Several near-unity transforms are applied to A to boost initial
// correlation at the given center angle.
//
// Return true if any changes made.
//
static bool Pretweaks(
	double	bestR,
	double	center,
	ThmRec	&thm,
	FILE*	flog )
{
	fprintf( flog, "Pretweaks start, best R=%.3f.\n", bestR );

	clock_t	t0 = StartTiming();
	bool	anychange = false;

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
							5, .02, center, thm, flog );
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
						.02, center, thm, flog );

		// Is it better than before?

		if( rbest > bestR ) {

			TForm	T, Tback = Tptwk;

			fprintf( flog, "PTWUsing %2d\n", vsel[selbest] );
			T.NUSelect( vsel[selbest], a );
			MultiplyTrans( Tptwk, T, Tback );
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
/* AngleScanWithTweaks ------------------------------------------- */
/* --------------------------------------------------------------- */

static double AngleScanWithTweaks(
	CorRec	&best,
	double	center,
	double	hlfwid,
	double	step,
	ThmRec	&thm,
	FILE*	flog )
{
	if( AngleScan( best, center, hlfwid, step, thm, flog )
		< GBL.ctx.RTRSH ) {

		if( GBL.mch.PRETWEAK ) {

			if( Pretweaks( best.R,
				(best.R > 0.0 ? best.A : center), thm, flog ) ) {

				AngleScan( best, center, hlfwid, step, thm, flog );
			}
		}
	}

	return best.R;
}

/* --------------------------------------------------------------- */
/* PeakHunt ------------------------------------------------------ */
/* --------------------------------------------------------------- */

#if 0
// Repeated parabola fits for peak on a narrowing angle range.
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
	"PeakHunt: Best: K=%d, R=%.3f, A=%.3f, X=%.3f, Y=%.3f.\n",
	k, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );

	return best.R;
}
#endif

// Bracket search for peak.
//
static double PeakHunt(
	CorRec	&best,
	double	hlfwid,
	ThmRec	&thm,
	FILE*	flog )
{
	CorRec	C,
			B0	= best;
	double	L	= best.A - hlfwid,
			R	= best.A + hlfwid,
			M, lr, rr;
	int		k	= 1;

	clock_t	t0 = StartTiming();

	RFromAngle( C, L, thm, flog );
	lr = C.R;

	RFromAngle( C, R, thm, flog );
	rr = C.R;

	RFromAngle( best, M = (L+R)/2.0, thm, flog );

	while( R - L > 0.0001 ) {

		double	a;

		++k;

		// move left up
		RFromAngle( C, a = (L+M)/2.0, thm, flog );

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
		RFromAngle( C, a = (M+R)/2.0, thm, flog );

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
	"PeakHunt: Best: K=%d, R=%.3f, A=%.3f, X=%.3f, Y=%.3f.\n",
	k, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );

	return best.R;
}

/* --------------------------------------------------------------- */
/* UsePriorAngles ------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool UsePriorAngles(
	CorRec	&best,
	int		nprior,
	ThmRec	&thm,
	FILE*	flog )
{
	fprintf( flog, "Approx: Using prior angles n=%d, med=%f\n",
	nprior, ang0 );

	if( AngleScanWithTweaks(
			best, ang0, GBL.ctx.HFANGPR, 0.1, thm, flog )
		< GBL.ctx.RTRSH ||
		PeakHunt( best, 0.3, thm, flog )
		< GBL.ctx.RTRSH ) {

		fprintf( flog,
		"FAIL: Approx: Prior angles R=%g below thresh=%g.\n",
		best.R, GBL.ctx.RTRSH );

		gErr = errLowRPrior;
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* DenovoBestAngle ----------------------------------------------- */
/* --------------------------------------------------------------- */

static bool DenovoBestAngle( CorRec &best, ThmRec &thm, FILE* flog )
{
	if( AngleScanWithTweaks(
			best, ang0, GBL.ctx.HFANGDN, 0.5, thm, flog )
		< GBL.ctx.RTRSH ||
		AngleScan( best, best.A, 1.0, 0.1, thm, flog )
		< GBL.ctx.RTRSH ||
		PeakHunt( best, 0.3, thm, flog )
		< GBL.ctx.RTRSH ) {

		fprintf( flog,
		"FAIL: Approx: Denovo R=%g below thresh=%g.\n",
		best.R, GBL.ctx.RTRSH );

		gErr = errLowRDenov;
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* TryTweaks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// The best transform is repeatedly multiplied with each of several
// near-unity tweak transforms to see if we get a better result.
//
static void TryTweaks( CorRec &best, ThmRec &thm, FILE* flog )
{
	vector<TForm>	twk( 12 );

	twk[0].NUSetScl( 1.005 );
	twk[1].NUSetScl( 0.995 );
	twk[2].NUSetXScl( 1.005 );
	twk[3].NUSetXScl( 0.995 );
	twk[4].NUSetYScl( 1.005 );
	twk[5].NUSetYScl( 0.995 );
	twk[6].NUSetXSkw( 0.005 );
	twk[7].NUSetXSkw( -.005 );
	twk[8].NUSetYSkw( 0.005 );
	twk[9].NUSetYSkw( -.005 );
	twk[10].NUSetRot( 0.005 );
	twk[11].NUSetRot( -.005 );

	fprintf( flog, "Tweaks start, best R=%.3f.\n", best.R );

	clock_t	t0 = StartTiming();

	for( int changed = true; changed; ) {

		changed = false;

		for( int i = 0; i < 12; ++i ) {

			CorRec			C;
			vector<Point>	ps = thm.ap;

			C.A = best.A;
			MultiplyTrans( C.T, best.T, twk[i] );

			C.T.Apply_R_Part( ps );

			C.R = CorrImages(
				flog, false, C.X, C.Y,
				ps, thm.av, thm.bp, thm.bv,
				BigEnough, (void*)thm.reqArea,
				EnoughPoints, (void*)thm.reqArea,
				0.0, GBL.ctx.NBMXHT, thm.ftc );

			fprintf( flog, "Tweak %d R=%.3f", i, C.R );

			if( C.R > best.R ) {

				fprintf( flog, "  *\n" );
				changed	= true;
				best	= C;
			}
			else
				fprintf( flog, "\n" );
		}
	}

	StopTiming( stdout, "TryTweaks", t0 );
}

/* --------------------------------------------------------------- */
/* FinishAtFullRes ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Get final offsets at full resolution.
//
static void FinishAtFullRes( CorRec &best, ThmRec &thm, FILE* flog )
{
	CorRec	b0 = best;
	int		ok;

	fprintf( flog,
	"Approx: LowRes  R=%.3f, X=%.3f, Y=%.3f.\n",
	best.R, best.X, best.Y );

	clock_t	t0 = StartTiming();

	vector<Point>	ps = thm.ap;

	best.T.Apply_R_Part( ps );

	best.R = CorrImages(
		flog, false, best.X, best.Y,
		ps, thm.av, thm.bp, thm.bv,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		0.0, GBL.ctx.NBMXHT, thm.ftc );

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

/* --------------------------------------------------------------- */
/* OffsetXY ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Let T(A) -> B be resolved into parts: R(A) + V -> B.
// Thus far we have found R,V for points (a') relative to
// the intersection origins, which we have labeled aO, bO
// in their respective coord systems. So we have found:
//
//		R( a' = a - aO ) + V -> b' = b - bO.
//
// We can rearrange this to:
//
//		R( a ) + ( V + bO - R(aO) ) -> b,
//
// which shows directly how to translate V to the image corner.
//
static void OffsetXY( CorRec &best, OlapRec &olp )
{
	best.T.Apply_R_Part( olp.aO );

	best.X += olp.bO.x - olp.aO.x;
	best.Y += olp.bO.y - olp.aO.y;
}

/* --------------------------------------------------------------- */
/* RecordSumSqDif ------------------------------------------------ */
/* --------------------------------------------------------------- */

// ******************************
// ALL CALCS USING SCALED PX DATA
// ******************************

static void SQD(
	double			&sqd,
	double			&prd,
	int				&N,
	const PixPair	&px,
	const TForm		&T )
{
	sqd	= 0.0;
	prd	= 0.0;
	N	= 0;

	const vector<double>&	av = *px.avs_vfy;
	const vector<double>&	bv = *px.bvs_vfy;

	vector<Point>	ap, bp;
	int				w  = px.ws,
					h  = px.hs,
					np = w * h;

// fill points

	MakeZeroBasedPoints( ap, w, h );
	bp = ap;

// sums

	for( int i = 0; i < np; ++i ) {

		Point	p = ap[i];

		T.Transform( p );

		if( p.x >= 0.0 && p.x < w-1 &&
			p.y >= 0.0 && p.y < h-1 ) {

			double	d = InterpolatePixel( p.x, p.y, bv, w );

			++N;
			prd += av[i] * d;
			d   -= av[i];
			sqd += d * d;
		}
	}
}


static void RecordSumSqDif(
	const PixPair	&px,
	const TForm		&T )
{
	CMutex	M;
	char	name[256];

	sprintf( name, "sqd_%d_%d", GBL.A.layer, GBL.B.layer );

	if( M.Get( name ) ) {

		sprintf( name, "SmSqDf_%d_@_%d.log",
			GBL.A.layer, GBL.B.layer );

		FILE *f;
		int  is;

		f = fopen( name, "r" );
		if( is = (f != NULL) )
			fclose( f );

		f = fopen( name, "a" );

		if( f ) {

			if( !is )
				fprintf( f, "TileA\tTileB\tSQ\tR\tN\tSQ/N\tR/N\n" );

			double	sqd, prd;
			int		N;

			SQD( sqd, prd, N, px, T );

			fprintf( f, "%d\t%d\t%f\t%f\t%d\t%f\t%f\n",
				GBL.A.tile, GBL.B.tile,
				sqd, prd, N, sqd/N, prd/N );
			fflush( f );
			fclose( f );
		}
	}

	M.Release();
}

/* --------------------------------------------------------------- */
/* RecordResult -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RecordResult( const CorRec &best )
{
	ThmPair	tpr;

	tpr.T	= best.T;
	tpr.A	= best.A;
	tpr.R	= best.R;
	tpr.err	= gErr;

	WriteThmPair( tpr, GBL.A.layer, GBL.A.tile, 1,
		GBL.B.layer, GBL.B.tile, 1 );
}

/* --------------------------------------------------------------- */
/* Failure ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Failure( CorRec &best )
{
	RecordResult( best );
	return false;
}

/* --------------------------------------------------------------- */
/* ApproximateMatch_NoCR ----------------------------------------- */
/* --------------------------------------------------------------- */

// Seek rough transform (mapping) from A to B.
//
// Return true and add one entry to guesses if successful.
//
// Also, output an entry in the ThmPair_lyrA_@_lyrB.txt file.
//
// Overview
// --------
// Matching in-plane gives much stronger signals since the
// features are essentially the same, but imaged twice. Cross-
// plane the features are a mix of nearly same and not same.
// Job parameters (required overlap area, thresholds) are cased
// out accordingly.
//
// We work on thumbnails, typically reduced by x8 on a side. The
// basic scheme is to calculate the correlation for angles in
// range +/-halfAng about estimated center angle. Best value is
// pretty distinct (simple maximum). That's best (rough) angle.
//
// Next the rough angle is refined by a narrowing bracket search
// for improved R until bracket width is 0.001 degree.
//
// The angles within a layer should really be the same. To improve
// angle consistency and to skip the wide angle sweep, we also read
// the file to get best angles from prior matches for this layer.
// Currently the test for skipping the sweep is just whether there
// are 4 or more entries. We could be smarter by looking at quality
// and/or clustering or distribution, but simple median works well
// on test data.
//
// After getting best angle on thumbnails we get final XY using
// full resolution images.
//
bool ApproximateMatch_NoCR(
	vector<TForm>	&guesses,
	const PixPair	&px,
	FILE*			flog )
{
	fprintf( flog, "\n---- Thumbnail matching ----\n" );

	CorRec	best;

	best.X	= -999.0;
	best.Y	= -999.0;
	best.A	= -999.0;
	best.R	= -999.0;

/* --------------------------------------------- */
/* Optionally reiterate result from existing log */
/* --------------------------------------------- */

//	return Orig( guesses );

//	return FromLog( guesses, flog );

/* ----------------------- */
/* Estimate starting angle */
/* ----------------------- */

	int	nPriorAngles = SetStartingAngle( flog );

/* ----------------------- */
/* Create image thumbnails */
/* ----------------------- */

	OlapRec	olp;
	ThmRec	thm;

	SelectSubimage( olp, px, flog );
	MakeThumbs( thm, olp, 8, flog );

/* --------------------- */
/* Debug the angle sweep */
/* --------------------- */

// -----------------------------------------------------------
//	DebugAngs( ang0, 45, .1, thm );
//	DebugAngs( ang0, 1, .01, thm );
//	exit( 42 );
// -----------------------------------------------------------

	if( dbgCor ) {
		RFromAngle( best, ang0, thm, flog );
		return false;
	}

/* ------------------------------------------- */
/* Search range of angles for best correlation */
/* ------------------------------------------- */

	if( nPriorAngles ) {

		if( !UsePriorAngles( best, nPriorAngles, thm, flog ) )
			return Failure( best );
	}
	else if( !DenovoBestAngle( best, thm, flog ) )
		return Failure( best );

/* ----------------------- */
/* Apply distortion tweaks */
/* ----------------------- */

// These tweaks can be enabled for EM. They don't make sense for
// optical, nor does RecordSumSqDif() show improvement in optical.

	if( GBL.mch.TWEAKS )
		TryTweaks( best, thm, flog );

/* --------------------------------- */
/* Always do this scaling adjustment */
/* --------------------------------- */

	best.X *= thm.scl;
	best.Y *= thm.scl;

/* -------------------------------- */
/* Final full resolution adjustment */
/* -------------------------------- */

// Using RecordSumSqDif() the observation is that this can improve
// starting overlap correlation as much as 15% for in-plane cases.
// For cross-plane the improvement is usually < 0.5% but it doesn't
// hurt, so we do it always.

	MakeThumbs( thm, olp, 1, flog );
	FinishAtFullRes( best, thm, flog );

/* ------------------------------------------ */
/* Translate from intersection to full coords */
/* ------------------------------------------ */

	OffsetXY( best, olp );

	best.T.SetXY( best.X, best.Y );

/* ------------- */
/* Report to log */
/* ------------- */

	//RecordSumSqDif( px, best.T );

	best.T.MulXY( px.scl );

	fprintf( flog, "Approx: Returning A=%f, R=%f, X=%f, Y=%f.\n",
	best.A, best.R, best.T.t[2], best.T.t[5] );

	fprintf( flog, "Approx: Best transform " );
	best.T.PrintTransform( flog );

/* --------------------- */
/* Constrain translation */
/* --------------------- */

// Stay close to original transform, assuming some preliminary
// alignment was done.

	if( GBL.mch.INPALN ) {

		TForm	T, Tinv, I;

		AToBTrans( T, GBL.A.t2i.T, GBL.B.t2i.T );
		fprintf( flog, "Approx: Orig transform " );
		T.PrintTransform( flog );

		InvertTrans( Tinv, T );
		MultiplyTrans( I, Tinv, best.T );
		fprintf( flog, "Approx: Idnt transform " );
		I.PrintTransform( flog );

		double	err = sqrt( I.t[2]*I.t[2] + I.t[5]*I.t[5] );

		fprintf( flog, "Approx: err = %g, max = %d.\n",
			err, GBL.mch.DINPUT );

		if( err > GBL.mch.DINPUT ) {

			fprintf( flog,
			"FAIL: Approx: Too different from Tinput"
			" err=%g, max=%d.\n",
			err, GBL.mch.DINPUT );

			return false;
		}
	}

/* --------------- */
/* Report to world */
/* --------------- */

	RecordResult( best );

	guesses.push_back( best.T );

	return true;
}


