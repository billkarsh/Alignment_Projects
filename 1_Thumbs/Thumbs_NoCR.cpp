

#include	"CGBL_Thumbs.h"
#include	"Thumbs.h"

#include	"Disk.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Timer.h"


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
		bx, by;

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

	if( olp.ow < GBL.mch.OLAP1D || olp.oh < GBL.mch.OLAP1D ) {
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
	const TForm		&T0,
	double			theta )
{
	double	c	= cos( theta ) * GBL.ctx.SCALE,
			s	= sin( theta ) * GBL.ctx.SCALE;
	TForm	ao( GBL.ctx.XSCALE * c, -GBL.ctx.YSCALE * s, 0.0,
				GBL.ctx.XSCALE * s,  GBL.ctx.YSCALE * c, 0.0 );

	MultiplyTrans( T, ao, T0 );
	T.Apply_R_Part( pts );
}

/* --------------------------------------------------------------- */
/* BigEnough ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool BigEnough( int sx, int sy, void *a )
{
	return	sx >= GBL.mch.OLAP1D &&
			sy >= GBL.mch.OLAP1D &&
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

	if( AngleScan( best, ang0, GBL.ctx.HFANGPR, 0.1, thm, flog )
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
	if( AngleScan( best, ang0, GBL.ctx.HFANGDN, 0.5, thm, flog )
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

// The best transform is repeatedly multiplied with each of eight
// near-unity tweak transforms to see if we get a better result.
//
static bool TryTweaks( CorRec &best, ThmRec &thm, FILE* flog )
{
	vector<TForm>	tweaks( 8 );

	tweaks[0] = TForm( 0.995,  0.0,   0.0,  0.0,   0.995, 0.0 ); // slightly smaller
	tweaks[1] = TForm( 1.005,  0.0,   0.0,  0.0,   1.005, 0.0 ); // slightly bigger
	tweaks[2] = TForm( 0.995,  0.0,   0.0,  0.0,   1.005, 0.0 ); // X smaller, Y bigger
	tweaks[3] = TForm( 1.005,  0.0,   0.0,  0.0,   0.995, 0.0 ); // slightly bigger
	tweaks[4] = TForm( 1.000,  0.005, 0.0,  0.005, 1.001, 0.0 ); // toe in
	tweaks[5] = TForm( 1.000, -0.005, 0.0, -0.005, 1.000, 0.0 ); // toe out
	tweaks[6] = TForm( 1.000,  0.005, 0.0, -0.005, 1.001, 0.0 ); // Small rotate
	tweaks[7] = TForm( 1.000, -0.005, 0.0,  0.005, 1.000, 0.0 ); // other way

	fprintf( flog, "Tweaks start, best R=%.3f.\n", best.R );

	clock_t	t0 = StartTiming();

	for( int changed = true; changed; ) {

		changed = false;

		for( int i = 0; i < 8; ++i ) {

			CorRec			C;
			vector<Point>	ps = thm.ap;

			C.A = best.A;
			MultiplyTrans( C.T, best.T, tweaks[i] );

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
/* Thumbs_NoCR --------------------------------------------------- */
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
bool Thumbs_NoCR( const PixPair &px, FILE* flog )
{
	fprintf( flog, "\n---- Thumbnail matching ----\n" );

	CorRec	best;

	best.X	= -999.0;
	best.Y	= -999.0;
	best.A	= -999.0;
	best.R	= -999.0;

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

#ifdef	CORR_DEBUG
	RFromAngle( best, ang0, thm, flog );
	return false;
#endif

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

		fprintf( flog, "Approx: err = %f, max = %d.\n",
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

	return true;
}


