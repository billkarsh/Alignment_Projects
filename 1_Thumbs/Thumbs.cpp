

#include	"CGBL_Thumbs.h"
#include	"Thumbs.h"

#include	"Disk.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum thmerrs {
	errOK			= 0,
	errLowQ			= 1,
	errLowRDenov	= 2,
	errLowRPrior	= 3
};

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	USE_Q		0
#define	QTOL_SWEEP	0.20
#define	QTOL_ZERO	0.00

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	vector<double>	v;
	vector<Point>	p;
	Point			O;			// subimage ref. point
	int				w, h;		// sub-img dims
} SubI;

typedef struct {
	SubI			a, b;
	int				nnegx, nnegy,	// ful-res search range
					nposx, nposy;
} OlapRec;

typedef struct {
	vector<double>	av, bv;
	vector<Point>	ap, bp;
	vector<CD>		ftc;			// fourier transform cache
	long			reqArea;
	int				scl;
	int				nnegx, nnegy,	// scaled search range
					nposx, nposy;
} ThmRec;

typedef struct {
	TForm			T;
	double			X, Y,
					A, Q, R;
} CorRec;

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static int		gErr		= errOK;
static long		min_2D_olap;
static double	ang0		= 0.0;
static double	halfAngDN;
static double	halfAngPR;
static double	scale		= 1.0;
static double	xscale		= 1.0;
static double	yscale		= 1.0;
static double	skew		= 0.0;
static double	qthresh;
static double	rthresh;






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

		ang0 = 180.0/PI * atan2( atob.t[3], atob.t[0] );
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

		GBL.thm.OLAP2D_XL = int(GBL.thm.OLAP2D_XL / fmax( c*c, s*s ));
	}

	return nprior;
}

/* --------------------------------------------------------------- */
/* SetLayerType -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Differentially set parameters according to job type.
//
static void SetLayerType()
{
	if( GBL.A.layer == GBL.B.layer ) {

		min_2D_olap	= GBL.thm.OLAP2D_SL;
		halfAngDN	= GBL.thm.HFANGDN_SL;
		halfAngPR	= GBL.thm.HFANGPR_SL;
		qthresh		= GBL.thm.QTRSH_SL;
		rthresh		= GBL.thm.RTRSH_SL;
	}
	else {

		min_2D_olap	= GBL.thm.OLAP2D_XL;
		halfAngDN	= GBL.thm.HFANGDN_XL;
		halfAngPR	= GBL.thm.HFANGPR_XL;
		scale		= GBL.thm.SCALE;
		xscale		= GBL.thm.XSCALE;
		yscale		= GBL.thm.YSCALE;
		skew		= GBL.thm.SKEW;
		qthresh		= GBL.thm.QTRSH_XL;
		rthresh		= GBL.thm.RTRSH_XL;
	}
}

/* --------------------------------------------------------------- */
/* WholeImageSubI ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WholeImageSubI(
	SubI					&S,
	const vector<double>	&v,
	const vector<Point>		&p,
	int						w )
{
	IBox	B;
	int		np = p.size();

	BBoxFromPoints( B, p );

	S.v.resize( np );
	S.p.resize( np );

	for( int i = 0; i < np; ++i ) {

		const Point&	P = p[i];

		S.v[i]		= v[int(P.x) + w * int(P.y)];
		S.p[i].x	= P.x - B.L;
		S.p[i].y	= P.y - B.B;
	}

	S.O	= Point( B.L, B.B );
	S.w = B.R - B.L + 1;
	S.h = B.T - B.B + 1;
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
	OlapRec				&olp,
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog )
{
	fprintf( flog,
	"Subimage: Using whole images, apix=%d, bpix=%d.\n",
	acr.pts.size(), bcr.pts.size() );

	WholeImageSubI( olp.a, *px.avs_aln, acr.pts, px.ws );
	WholeImageSubI( olp.b, *px.bvs_aln, bcr.pts, px.ws );

	olp.nnegx	= olp.a.w - GBL.thm.OLAP1D;
	olp.nposx	= olp.b.w - GBL.thm.OLAP1D;
	olp.nnegy	= olp.a.h - GBL.thm.OLAP1D;
	olp.nposy	= olp.b.h - GBL.thm.OLAP1D;

	return true;
}

/* --------------------------------------------------------------- */
/* SubimageSubI -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SubimageSubI(
	SubI					&S,
	const IBox				&Bolap,
	const vector<double>	&v,
	const vector<Point>		&p,
	int						w )
{
// 1st push_back data within intersection

	int	np = p.size();

	for( int i = 0; i < np; ++i ) {

		const Point&	P = p[i];

		if( P.x >= Bolap.L && P.x <= Bolap.R &&
			P.y >= Bolap.B && P.y <= Bolap.T ) {

			S.v.push_back( v[int(P.x) + w * int(P.y)] );
			S.p.push_back( P );
		}
	}

	np = S.p.size();

	S.v.resize( np );
	S.p.resize( np );

// Get true origin and dimension for selected points

	IBox	B;

	BBoxFromPoints( B, S.p );

	S.O	= Point( B.L, B.B );
	S.w = B.R - B.L + 1;
	S.h = B.T - B.B + 1;

// Finally make points zero-based

	for( int i = 0; i < np; ++i ) {

		S.p[i].x -= B.L;
		S.p[i].y -= B.B;
	}
}

/* --------------------------------------------------------------- */
/* SelectSubimage ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Describe the conservative a/b intersection region based upon
// data from TileToImage.txt.
//
static bool SelectSubimage(
	OlapRec				&olp,
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog )
{
// Use whole images if cross layer

	if( GBL.A.layer != GBL.B.layer )
		return WholeImage( olp, px, acr, bcr, flog );

// Get displacements of a in b-system

	int	dx, dy;

	{
		TForm	atob;
		Point	delta;

		AToBTrans( atob, GBL.A.t2i.T, GBL.B.t2i.T );
		atob.Transform( delta );

		dx = int(delta.x) / px.scl;
		dy = int(delta.y) / px.scl;

		if( GBL.thm.SLOPPY_SL ) {
			dx /= 2;
			dy /= 2;
		}
	}

// Use offsets and image size to determine:
// - {Ba, Bb} overlap boxes.

	IBox	Ba, Bb;
	int		w = px.ws,
			h = px.hs,
			ow, oh;

	BoxesFromShifts( Ba, Bb, w, h, w, h, dx, dy );

	ow = Ba.R - Ba.L + 1;
	oh = Ba.T - Ba.B + 1;

// Double-check that there was sufficient overlap

	if( ow < GBL.thm.OLAP1D || oh < GBL.thm.OLAP1D ) {
		fprintf( flog, "Subimage: Overlap looks small.\n" );
		return WholeImage( olp, px, acr, bcr, flog );
	}

	fprintf( flog,
	"Subimage: Using intersection, w=%d, h=%d, pix=%d.\n",
	ow, oh, ow * oh );

// Fill lists with ConnRegion data within the intersection

	SubimageSubI( olp.a, Ba, *px.avs_aln, acr.pts, w );
	SubimageSubI( olp.b, Bb, *px.bvs_aln, bcr.pts, w );

// Finish

	if( GBL.thm.SLOPPY_SL ) {
		olp.nnegx = olp.nposx = ow;
		olp.nnegy = olp.nposy = oh;
	}
	else {
		olp.nnegx = olp.nposx = ow / 2;
		olp.nnegy = olp.nposy = oh / 2;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* MakeThumbs ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool MakeThumbs(
	ThmRec			&thm,
	const OlapRec	&olp,
	int				decfactor,
	FILE*			flog )
{
	thm.av		= olp.a.v;
	thm.bv		= olp.b.v;
	thm.ap		= olp.a.p;
	thm.bp		= olp.b.p;
	thm.ftc.clear();
	thm.reqArea	= min_2D_olap;
	thm.scl		= decfactor;
	thm.nnegx	= olp.nnegx;
	thm.nnegy	= olp.nnegy;
	thm.nposx	= olp.nposx;
	thm.nposy	= olp.nposy;

	if( decfactor > 1 ) {

		DecimateVector( thm.ap, thm.av, olp.a.w, olp.a.h, decfactor );
		DecimateVector( thm.bp, thm.bv, olp.b.w, olp.b.h, decfactor );

		thm.reqArea	/= decfactor * decfactor;
		thm.nnegx	/= decfactor;
		thm.nnegy	/= decfactor;
		thm.nposx	/= decfactor;
		thm.nposy	/= decfactor;

		fprintf( flog,
		"Thumbs: After decimation %d pts, reqArea %d, thmscl %d.\n",
		thm.ap.size(), thm.reqArea, thm.scl );
	}

	Normalize( thm.av );
	Normalize( thm.bv );

	thm.nnegx	-= 1;
	thm.nnegy	-= 1;

	if( thm.nnegx <= 0 || thm.nposx <= 0 ||
		thm.nnegy <= 0 || thm.nposy <= 0 ) {

		fprintf( flog,
		"Thumbs: Search range below [+1,+1) x:[%d %d) y:[%d %d).\n",
		thm.nnegx, thm.nposx, thm.nnegy, thm.nposy );

		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* RotatePoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RotatePoints(
	vector<Point>	&pts,
	TForm			&T,
	const TForm		&T0,
	double			theta )
{
	double	c	= cos( theta ) * scale,
			s	= sin( theta ) * scale;
	TForm	ao( xscale * c, -yscale * s, 0.0,
				xscale * s,  yscale * c, 0.0 );

	MultiplyTrans( T, ao, T0 );
	T.Apply_R_Part( pts );
}

/* --------------------------------------------------------------- */
/* BigEnough ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool BigEnough( int sx, int sy, void *a )
{
	return (long)sx * sy > (long)a;
}

/* --------------------------------------------------------------- */
/* EnoughPoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool EnoughPoints( int count1, int count2, void *a )
{
	return count1 > (long)a && count2 > (long)a;
}

/* --------------------------------------------------------------- */
/* QFromAngle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void QFromAngle(
	CorRec	&C,
	double	a,
	double	qtol,
	ThmRec	&thm,
	FILE*	flog )
{
	TForm			Tskew( 1.0, 0.0, 0.0, skew, 1.0, 0.0 );
	vector<Point>	ps = thm.ap;

	C.A = a;

	RotatePoints( ps, C.T, Tskew, a * PI/180.0 );

#if USE_Q == 1
	C.R = CorrPatchesMaxQ(
		flog, false, C.Q, C.X, C.Y,
		ps, thm.av, thm.bp, thm.bv,
		thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea, thm.ftc );
#else
	C.Q = C.R = CorrPatchesMaxR(
		flog, false, C.X, C.Y,
		ps, thm.av, thm.bp, thm.bv,
		thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		qtol, thm.ftc );
#endif
}

/* --------------------------------------------------------------- */
/* DebugAngs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DebugAngs(
	double	center,
	double	hlfwid,
	double	step,
	double	qtol,
	ThmRec &thm )
{
	char	file[256];

	sprintf( file, "angs_%d_%d_@_%d_%d.log",
		GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile );
	FILE	*f = fopen( file, "w" );

	fprintf( f, "Deg\tR\tQ\tX\tY\n" );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		QFromAngle( C, a, qtol, thm, stdout );

		fprintf( f, "%.3f\t%.4f\t%.4f\t%.3f\t%.3f\n",
			a, C.R, C.Q, C.X, C.Y );
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
	"%s: Q=%6.3f, R=%6.3f, A=%8.3f, X=%8.2f, Y=%8.2f\n",
	label, C.Q, C.R, C.A, C.X, C.Y );
}


static void AngleScan(
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

	QFromAngle( best, center, QTOL_SWEEP, thm, flog );
	RecordAngle( flog, "Center", best );

	for( double a = center-hlfwid; a <= center+hlfwid; a += step ) {

		CorRec	C;

		QFromAngle( C, a, QTOL_SWEEP, thm, flog );
		RecordAngle( flog, "  Scan", C );

		if( C.Q > best.Q )
			best = C;
	}

	RecordAngle( flog, "  Best", best );

	StopTiming( stdout, "AngleScan", t0 );
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
	return x1 + d * (y0 - y2) / (2 * (y0 + y2 - y1 - y1));
}

/* --------------------------------------------------------------- */
/* PeakHunt ------------------------------------------------------ */
/* --------------------------------------------------------------- */

#if 0
// Repeated parabola fits for peak on a narrowing angle range.
//
static void PeakHunt(
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
		QFromAngle( C, L, QTOL_ZERO, thm, flog );
		y0 = C.Q;

		// right
		QFromAngle( C, R, QTOL_ZERO, thm, flog );
		y2 = C.Q;

		// middle
		QFromAngle( C, x1 = (L+R)/2.0, QTOL_ZERO, thm, flog );

		// estimate peak position
		C.A = NewXFromParabola( x1, (R-L)/2.0, y0, C.Q, y2 );

		// sanity check
		if( C.A <= L || C.A >= R )
			break;

		// compete estimate against best so far
		QFromAngle( C, C.A, QTOL_ZERO, thm, flog );

		if( C.Q > best.Q )
			best = C;

		//fprintf( flog,
		//"*** [%f,%f] [%f,%f] [%f,%f] <newx=%f>.\n",
		//L,y0, best.A,best.Q, R,y2, x1 );

		x1	= fmin( best.A - L, R - best.A ) / 4.0;
		L	= best.A - x1;
		R	= best.A + x1;
	}

	if( B0.R > best.R && fabs( B0.Q - best.Q ) < 0.005 )
		best = B0;

	fprintf( flog,
	"PeakHunt: Best: K=%d, Q=%.3f, R=%.3f, A=%.3f, X=%.3f, Y=%.3f.\n",
	k, best.Q, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );
}
#endif

// Bracket search for peak.
//
static void PeakHunt(
	CorRec	&best,
	double	hlfwid,
	ThmRec	&thm,
	FILE*	flog )
{
	CorRec	C,
			B0	= best;
	double	L	= best.A - hlfwid,
			R	= best.A + hlfwid,
			M, lq, rq;
	int		k	= 1;

	clock_t	t0 = StartTiming();

	QFromAngle( C, L, QTOL_ZERO, thm, flog );
	lq = C.Q;

	QFromAngle( C, R, QTOL_ZERO, thm, flog );
	rq = C.Q;

	QFromAngle( best, M = (L+R)/2.0, QTOL_ZERO, thm, flog );

	while( R - L > 0.0001 ) {

		double	a;

		++k;

		// move left up
		QFromAngle( C, a = (L+M)/2.0, QTOL_ZERO, thm, flog );

		if( C.Q >= best.Q ) {
			rq		= best.Q;
			R		= M;
			best	= C;
			M		= a;
		}
		else {
			lq		= C.Q;
			L		= a;
		}

		// move right back
		QFromAngle( C, a = (M+R)/2.0, QTOL_ZERO, thm, flog );

		if( C.Q >= best.Q ) {
			lq		= best.Q;
			L		= M;
			best	= C;
			M		= a;
		}
		else {
			rq		= C.Q;
			R		= a;
		}
	}

	if( B0.Q > best.Q )
		best = B0;

	fprintf( flog,
	"PeakHunt: Best: K=%d, Q=%.3f, R=%.3f, A=%.3f, X=%.3f, Y=%.3f.\n",
	k, best.Q, best.R, best.A, best.X, best.Y );

	StopTiming( stdout, "PeakHunt", t0 );
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
	fprintf( flog, "Thumbs: Using prior angles n=%d, med=%f\n",
	nprior, ang0 );

	AngleScan( best, ang0, halfAngPR, 0.1, thm, flog );

	if( best.Q < 0.001 ) {
		gErr = errLowQ;
		return false;
	}

	PeakHunt( best, 0.15, thm, flog );

	if( best.R < rthresh ) {

		fprintf( flog,
		"Thumbs: R=%g below thresh=%g.\n", best.R, rthresh );

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
	AngleScan( best, ang0, halfAngDN, 0.5, thm, flog );

	if( best.Q < qthresh ) {
		gErr = errLowQ;
		return false;
	}

	AngleScan( best, best.A, 1.0, 0.1, thm, flog );

	PeakHunt( best, 0.15, thm, flog );

	if( best.R < rthresh ) {

		fprintf( flog,
		"Thumbs: R=%g below thresh=%g.\n", best.R, rthresh );

		gErr = errLowRDenov;
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* TryTweaks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// The best transform is multiplied with each of eight near-unity
// tweak transforms to see if we get a better result.
//
static bool TryTweaks( CorRec &best, ThmRec &thm, FILE* flog )
{
	vector<TForm>	tweaks(8);
	TForm			tcur;

	tcur = best.T;

	tweaks[0] = TForm( 0.995,  0.0,   0.0,  0.0,   0.995, 0.0 ); // slightly smaller
	tweaks[1] = TForm( 1.005,  0.0,   0.0,  0.0,   1.005, 0.0 ); // slightly bigger
	tweaks[2] = TForm( 0.995,  0.0,   0.0,  0.0,   1.005, 0.0 ); // X smaller, Y bigger
	tweaks[3] = TForm( 1.005,  0.0,   0.0,  0.0,   0.995, 0.0 ); // slightly bigger
	tweaks[4] = TForm( 1.000,  0.005, 0.0,  0.005, 1.001, 0.0 ); // toe in
	tweaks[5] = TForm( 1.000, -0.005, 0.0, -0.005, 1.000, 0.0 ); // toe out
	tweaks[6] = TForm( 1.000,  0.005, 0.0, -0.005, 1.001, 0.0 ); // Small rotate
	tweaks[7] = TForm( 1.000, -0.005, 0.0,  0.005, 1.000, 0.0 ); // other way

	fprintf( flog,
	"Thumbs: Tweaks start, best Q=%.3f, R=%.3f.\n",
	best.Q, best.R );

	for( int i = 0; i < 8; ++i ) {

		CorRec			C;
		vector<Point>	ps = thm.ap;

		C.A = best.A;
		MultiplyTrans( C.T, tcur, tweaks[i] );

		C.T.Apply_R_Part( ps );

#if USE_Q == 1
		C.R = CorrPatchesMaxQ(
			flog, false, C.Q, C.X, C.Y,
			ps, thm.av, thm.bp, thm.bv,
			thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
			BigEnough, (void*)thm.reqArea,
			EnoughPoints, (void*)thm.reqArea, thm.ftc );
#else
		C.Q = C.R = CorrPatchesMaxR(
			flog, false, C.X, C.Y,
			ps, thm.av, thm.bp, thm.bv,
			thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
			BigEnough, (void*)thm.reqArea,
			EnoughPoints, (void*)thm.reqArea,
			0.0, thm.ftc );
#endif

		fprintf( flog, "Tweak %d Q=%.3f, R=%.3f.\n", i, C.Q, C.R );

		if( C.R > best.R ) {

			fprintf( flog,
			"Thumbs: Using tweak[%d] Q=%.3f, R=%.3f;"
			" was Q=%.3f, R=%.3f.\n",
			i, C.Q, C.R, best.Q, best.R );

			best = C;
		}
	}
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
	"Thumbs: LowRes  Q=%.3f, R=%.3f, X=%.3f, Y=%.3f.\n",
	best.Q, best.R, best.X, best.Y );

	clock_t	t0 = StartTiming();

	vector<Point>	ps = thm.ap;

	best.T.Apply_R_Part( ps );

#if USE_Q == 1
	best.R = CorrPatchesMaxQ(
		flog, false, best.Q, best.X, best.Y,
		ps, thm.av, thm.bp, thm.bv,
		thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea, thm.ftc );
#else
	best.Q = best.R = CorrPatchesMaxR(
		flog, false, best.X, best.Y,
		ps, thm.av, thm.bp, thm.bv,
		thm.nnegx, thm.nposx, thm.nnegy, thm.nposy,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		0.0, thm.ftc );
#endif

	ok = (fabs( best.X - b0.X ) <= 20)
	  && (fabs( best.Y - b0.Y ) <= 20);

	fprintf( flog,
	"Thumbs: FullRes Q=%.3f, R=%.3f, X=%.3f, Y=%.3f, use=%c.\n",
	best.Q, best.R, best.X, best.Y, (ok ? 'Y' : 'N') );

// Always report the thumb-Q instead of the fullres value because
// fullres values are scaled differently. We prefer all values in
// the file to be comparable.

	if( ok )
		best.Q = b0.Q;
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
	best.T.Apply_R_Part( olp.a.O );

	best.X += olp.b.O.x - olp.a.O.x;
	best.Y += olp.b.O.y - olp.a.O.y;
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

static void RecordResult( const CorRec &best, int aid, int bid )
{
	ThmPair	tpr;

	tpr.T	= best.T;
	tpr.A	= best.A;
	tpr.Q	= best.Q;
	tpr.R	= best.R;
	tpr.err	= gErr;

	WriteThmPair( tpr, GBL.A.layer, GBL.A.tile, aid,
		GBL.B.layer, GBL.B.tile, bid );
}

/* --------------------------------------------------------------- */
/* Failure ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Failure( CorRec &best, int aid, int bid )
{
	RecordResult( best, aid, bid );
	return false;
}

/* --------------------------------------------------------------- */
/* Thumbs -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Seek rough transform (mapping) from A to B.
//
// The output is an entry in the ThmPair_lyrA_@_lyrB.txt file.
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
// basic scheme is to calculate the correlation {Q,R} for angles
// in range +/-halfAng about estimated center angle. Best Q is
// pretty distinct (simple maximum). That's best (rough) angle.
//
// Next the rough angle is refined by a narrowing bracket search
// for improved Q until bracket width is 0.001 degree.
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
bool Thumbs(
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog )
{
	fprintf( flog, "\n---- Thumbs ----\n" );

	CorRec	best;

	best.X	= -999.0;
	best.Y	= -999.0;
	best.A	= -999.0;
	best.Q	= -999.0;
	best.R	= -999.0;

/* ----------------------- */
/* Estimate starting angle */
/* ----------------------- */

	int	nPriorAngles = SetStartingAngle( flog );

/* ----------------------------------------------------- */
/* Set parameters according to same- or cross-layer type */
/* ----------------------------------------------------- */

	SetLayerType();

/* ----------------------- */
/* Create image thumbnails */
/* ----------------------- */

	OlapRec	olp;
	ThmRec	thm;

	SelectSubimage( olp, px, acr, bcr, flog );

	if( !MakeThumbs( thm, olp, 8, flog ) )
		return false;

/* --------------------- */
/* Debug the angle sweep */
/* --------------------- */

// -----------------------------------------------------------
//	DebugAngs( ang0, halfAngPR, .1, QTOL_SWEEP, thm );
//	DebugAngs( ang0, 1, .01, 0, thm );
//	exit( 42 );
// -----------------------------------------------------------

#ifdef	CORR_DEBUG
	QFromAngle( best, ang0, QTOL_ZERO, thm, flog );
	return false;
#endif

/* ------------------------------------------- */
/* Search range of angles for best correlation */
/* ------------------------------------------- */

	if( nPriorAngles ) {

		if( !UsePriorAngles( best, nPriorAngles, thm, flog ) )
			return Failure( best, acr.id, bcr.id );
	}
	else if( !DenovoBestAngle( best, thm, flog ) )
		return Failure( best, acr.id, bcr.id );

/* ----------------------- */
/* Apply distortion tweaks */
/* ----------------------- */

// These tweaks can be enabled for EM. They don't make sense for
// optical, nor does RecordSumSqDif() show improvement in optical.

	if( GBL.thm.TWEAKS )
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

	fprintf( flog, "Thumbs: Returning A=%f, R=%f, X=%f, Y=%f.\n",
	best.A, best.R, best.T.t[2], best.T.t[5] );

	fprintf( flog, "Thumbs: Best transform " );
	best.T.PrintTransform( flog );

/* --------------------- */
/* Constrain translation */
/* --------------------- */

// Stay close to original transform, assuming some preliminary
// alignment was done.

	if( GBL.thm.INPALN ) {

		TForm	T, Tinv, I;

		AToBTrans( T, GBL.A.t2i.T, GBL.B.t2i.T );
		fprintf( flog, "Thumbs: Orig transform " );
		T.PrintTransform( flog );

		InvertTrans( Tinv, T );
		MultiplyTrans( I, Tinv, best.T );
		fprintf( flog, "Thumbs: Idnt transform " );
		I.PrintTransform( flog );

		double	err = sqrt( I.t[2]*I.t[2] + I.t[5]*I.t[5] );

		fprintf( flog, "Thumbs: err = %f, max = %d.\n",
			err, GBL.thm.DINPUT );

		if( err > GBL.thm.DINPUT )
			return false;
	}

/* --------------- */
/* Report to world */
/* --------------- */

	RecordResult( best, acr.id, bcr.id );

	return true;
}


