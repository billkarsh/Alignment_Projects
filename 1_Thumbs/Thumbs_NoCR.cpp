

#include	"CGBL_Thumbs.h"
#include	"Thumbs.h"

#include	"CThmScan.h"
#include	"Disk.h"
#include	"Maths.h"
#include	"Debug.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Zero angle hypothesis: That within a layer the tile-tile
// angles are all zero, unless sample actually distorts during
// imaging, which may be true in EM. However, zero may still be
// good enough angle estimator.
//
#define	ZEROANG	1

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	vector<double>	av, bv;
	vector<Point>	ap, bp;
	Point			aO, bO;			// subimage ref. points
	int				ow, oh;			// common sub-img dims
} OlapRec;

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

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

// In-layer zero angle hypothesis

#if ZEROANG

	if( GBL.A.layer == GBL.B.layer ) {

		ang0	= 0.0;
		return 1;
	}

#endif

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

	if( !nprior )
		ang0 = 180.0/PI * RadiansFromAffine( GBL.Tab[0] );

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
	"Subimage: Using whole images, pix=%d\n", w * h );

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
		Point	delta;

		GBL.Tab[0].Transform( delta );

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
		fprintf( flog, "Subimage: 1D overlap too small.\n" );
		return WholeImage( olp, px, flog );
	}

	np = olp.ow * olp.oh;

	fprintf( flog,
	"Subimage: Using intersection, w=%d, h=%d, pix=%d\n",
	olp.ow, olp.oh, np );

// Corners

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

static bool MakeThumbs(
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
	thm.olap1D	= GBL.ctx.OLAP1D;
	thm.scl		= decfactor;

	if( decfactor > 1 ) {

		DecimateVector( thm.ap, thm.av, olp.ow, olp.oh, decfactor );
		DecimateVector( thm.bp, thm.bv, olp.ow, olp.oh, decfactor );

		thm.reqArea	/= decfactor * decfactor;
		thm.olap1D  /= decfactor;

		fprintf( flog,
		"Thumbs: After decimation %d pts, reqArea %d, thmscl %d\n",
		thm.ap.size(), thm.reqArea, thm.scl );
	}

	double	sd = Normalize( thm.av );

	if( !sd || !isfinite( sd ) ) {

		fprintf( flog,
		"Thumbs: Image A intersection region has stdev: %f\n", sd );

		return false;
	}

	sd = Normalize( thm.bv );

	if( !sd || !isfinite( sd ) ) {

		fprintf( flog,
		"Thumbs: Image B intersection region has stdev: %f\n", sd );

		return false;
	}

	return true;
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

static void RecordResult( const CorRec &best, int err )
{
	ThmPair	tpr;

	tpr.T	= best.T;
	tpr.A	= best.A;
	tpr.R	= best.R;
	tpr.err	= err;

	WriteThmPair( tpr, GBL.A.layer, GBL.A.tile, 1,
		GBL.B.layer, GBL.B.tile, 1 );
}

/* --------------------------------------------------------------- */
/* Failure ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Failure( CorRec &best, int err )
{
	RecordResult( best, err );
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
	CThmScan	S;
	CorRec		best;

	S.Initialize( flog, best );
	S.SetTuser(
		GBL.ctx.SCALE, GBL.ctx.XSCALE, GBL.ctx.YSCALE,
		0, GBL.ctx.SKEW );
	S.SetRThresh( GBL.ctx.RTRSH );
	S.SetNbMaxHt( GBL.ctx.NBMXHT );
	S.SetSweepConstXY( true );
	S.SetSweepPretweak( GBL.mch.PRETWEAK );
	S.SetUseCorrR( false );
	S.SetDisc( 0, 0, -1, -1 );

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

	if( !MakeThumbs( thm, olp, GBL.mch.THMDEC, flog ) )
		return false;

/* --------------------- */
/* Debug the angle sweep */
/* --------------------- */

// -----------------------------------------------------------
//	S.DebugAngs( GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile,
//		ang0, 45, .1, thm );
//	S.DebugAngs( GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile,
//		ang0, 1, .01, thm );
//	exit( 42 );
// -----------------------------------------------------------

#if ZEROANG

	if( dbgCor && GBL.A.layer != GBL.B.layer ) {
		S.RFromAngle( best, ang0, thm );
		return false;
	}

#else

	if( dbgCor ) {
		S.RFromAngle( best, ang0, thm );
		return false;
	}

#endif

/* --------------------------------------------- */
/* Zero angle search constrained by stage coords */
/* --------------------------------------------- */

#if ZEROANG

	if( GBL.A.layer == GBL.B.layer ) {

		Point	delta, TaO = olp.aO;
		int		Ox, Oy, Rx;

		GBL.Tab[0].Transform( delta );
		GBL.Tab[0].Apply_R_Part( TaO );

		Ox = ROUND((delta.x / px.scl - olp.bO.x + TaO.x) / thm.scl);
		Oy = ROUND((delta.y / px.scl - olp.bO.y + TaO.y) / thm.scl);
		Rx = GBL.ctx.DINPUT / (thm.scl * px.scl);

		fprintf( flog, "SetDisc( %d, %d, %d, %d )\n", Ox, Oy, Rx, Rx );

		S.SetUseCorrR( true );
		S.SetDisc( Ox, Oy, Rx, Rx );
		S.RFromAngle( best, 0, thm );

		fprintf( flog,
		"Initial: R=%6.3f, A=%8.3f, X=%8.2f, Y=%8.2f\n",
		best.R, best.A, best.X, best.Y );

		if( dbgCor )
			return false;

		if( best.R < GBL.ctx.RTRSH ) {

			if( GBL.mch.PRETWEAK ) {

				S.Pretweaks( best.R, 0, thm );
				S.RFromAngle( best, 0, thm );

				fprintf( flog,
				"Tweaked: R=%6.3f, A=%8.3f, X=%8.2f, Y=%8.2f\n",
				best.R, best.A, best.X, best.Y );
			}

			if( best.R < GBL.ctx.RTRSH ) {

				fprintf( flog,
				"FAIL: Approx: Zeroang R=%g below thresh=%g\n",
				best.R, GBL.ctx.RTRSH );

				return Failure( best, 2 );
			}
		}

		Point	dS(
				(best.X - Ox) * (thm.scl * px.scl),
				(best.Y - Oy) * (thm.scl * px.scl) );

		fprintf( flog, "Peak-Disc: dR %d dX %d dY %d\n",
		int(sqrt( dS.RSqr() )), int(dS.x), int(dS.y) );
	}
	else

#endif

/* ------------------------------------------- */
/* Search range of angles for best correlation */
/* ------------------------------------------- */

	if( nPriorAngles ) {

		if( !S.UsePriorAngles( best, nPriorAngles, ang0,
				GBL.ctx.HFANGPR, thm ) ) {

			return Failure( best, S.GetErr() );
		}
	}
	else if( !S.DenovoBestAngle( best,
				ang0, GBL.ctx.HFANGDN, 0.5, thm ) ) {

		return Failure( best, S.GetErr() );
	}

/* ----------------------- */
/* Apply distortion tweaks */
/* ----------------------- */

// These tweaks can be enabled for EM. They don't make sense for
// optical, nor does RecordSumSqDif() show improvement in optical.

	if( GBL.mch.TWEAKS )
		S.PostTweaks( best, thm );

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

	if( !MakeThumbs( thm, olp, 1, flog ) )
		return false;

	S.FinishAtFullRes( best, thm );

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

	fprintf( flog, "Approx: Returning A=%f, R=%f, X=%f, Y=%f\n",
	best.A, best.R, best.T.t[2], best.T.t[5] );

	fprintf( flog, "Approx: Best transform " );
	best.T.PrintTransform( flog );

/* --------------------- */
/* Constrain translation */
/* --------------------- */

// Stay close to original transform, assuming some preliminary
// alignment was done.

#if ZEROANG
// Soln already known close enough to input
	if( GBL.A.layer == GBL.B.layer )
		GBL.ctx.INPALN = false;
#endif

	{
		TForm	Tinv, I;

		fprintf( flog, "Approx: Orig transform " );
		GBL.Tab[0].PrintTransform( flog );

		InvertTrans( Tinv, GBL.Tab[0] );
		MultiplyTrans( I, Tinv, best.T );
		fprintf( flog, "Approx: Idnt transform " );
		I.PrintTransform( flog );

		double	err = sqrt( I.t[2]*I.t[2] + I.t[5]*I.t[5] );

		fprintf( flog, "Approx: err = %g, max = %d\n",
			err, GBL.ctx.DINPUT );

		if( GBL.ctx.INPALN && err > GBL.ctx.DINPUT ) {

			fprintf( flog,
			"FAIL: Approx: Too different from Tinput"
			" err=%g, max=%d\n",
			err, GBL.ctx.DINPUT );

			return false;
		}
	}

/* --------------- */
/* Report to world */
/* --------------- */

	RecordResult( best, S.GetErr() );

	return true;
}


