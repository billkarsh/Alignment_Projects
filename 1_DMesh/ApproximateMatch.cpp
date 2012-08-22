

#include	"CGBL_dmesh.h"
#include	"ApproximateMatch.h"

#include	"CThmScan.h"
#include	"Disk.h"
#include	"Maths.h"
#include	"Geometry.h"
#include	"Debug.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

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
} OlapRec;

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static double	ang0	= 0.0;






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

static bool FromLog( vector<TForm> &G, FILE* flog, int aid, int bid )
{
	ThmPair	tpr;
	int		ok;

	ok = ReadThmPair( tpr, GBL.A.layer, GBL.A.tile, aid,
			GBL.B.layer, GBL.B.tile, bid, flog )
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
	"Subimage: Using whole images, apix=%d, bpix=%d\n",
	acr.pts.size(), bcr.pts.size() );

	WholeImageSubI( olp.a, *px.avs_aln, acr.pts, px.ws );
	WholeImageSubI( olp.b, *px.bvs_aln, bcr.pts, px.ws );

	return true;
}

/* --------------------------------------------------------------- */
/* SubimageSubI -------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool SubimageSubI(
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

	if( (np = S.p.size()) <= GBL.ctx.OLAP2D )
		return false;

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

	return true;
}

/* --------------------------------------------------------------- */
/* SelectSubimage ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Describe the conservative a/b intersection region based upon
// data from TileToImage.txt.
//
// Return true if non-empty olap.
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

		if( GBL.mch.SLOPPY_SL ) {
			dx /= 2;
			dy /= 2;
		}
	}

// Use offsets and image size to determine:
// - {Ba, Bb} overlap boxes.

	IBox	Ba, Bb;
	int		w = px.ws,
			h = px.hs,
			ow, oh, min1d;

	BoxesFromShifts( Ba, Bb, w, h, w, h, dx, dy );

	ow = Ba.R - Ba.L + 1;
	oh = Ba.T - Ba.B + 1;

// Double-check that there was sufficient overlap

	min1d = max( GBL.ctx.OLAP1D, 8 );

	if( ow < min1d || oh < min1d ) {
		fprintf( flog, "Subimage: 1D overlap too small.\n" );
		return WholeImage( olp, px, acr, bcr, flog );
	}

	fprintf( flog,
	"Subimage: Using intersection, w=%d, h=%d, pix=%d\n",
	ow, oh, ow * oh );

// Fill lists with ConnRegion data within the intersection

	if( !SubimageSubI( olp.a, Ba, *px.avs_aln, acr.pts, w ) ||
		!SubimageSubI( olp.b, Bb, *px.bvs_aln, bcr.pts, w ) ) {

		fprintf( flog, "Subimage: 2D overlap too small.\n" );
		return false;
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
	thm.reqArea	= GBL.ctx.OLAP2D;
	thm.olap1D	= GBL.ctx.OLAP1D;
	thm.scl		= decfactor;

	if( decfactor > 1 ) {

		DecimateVector( thm.ap, thm.av, olp.a.w, olp.a.h, decfactor );
		DecimateVector( thm.bp, thm.bv, olp.b.w, olp.b.h, decfactor );

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

static void RecordResult(
	const CorRec	&best,
	int				aid,
	int				bid,
	int				err )
{
	ThmPair	tpr;

	tpr.T	= best.T;
	tpr.A	= best.A;
	tpr.R	= best.R;
	tpr.err	= err;

	WriteThmPair( tpr, GBL.A.layer, GBL.A.tile, aid,
		GBL.B.layer, GBL.B.tile, bid );
}

/* --------------------------------------------------------------- */
/* Failure ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Failure( CorRec &best, int aid, int bid, int err )
{
	RecordResult( best, aid, bid, err );
	return false;
}

/* --------------------------------------------------------------- */
/* ApproximateMatch ---------------------------------------------- */
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
bool ApproximateMatch(
	vector<TForm>		&guesses,
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog )
{
	CThmScan	S;
	CorRec		best;

	S.Initialize( flog, best );
	S.SetTuser(
		GBL.ctx.SCALE, GBL.ctx.XSCALE, GBL.ctx.YSCALE,
		0, GBL.ctx.SKEW );
	S.SetRThresh( GBL.ctx.RTRSH );
	S.SetNbMaxHt( GBL.ctx.NBMXHT );
	S.SetSweepType( false, GBL.mch.PRETWEAK );
	S.SetUseCorrR( false );
	S.SetDisc( 0, 0, -1, -1 );

/* --------------------------------------------- */
/* Optionally reiterate result from existing log */
/* --------------------------------------------- */

//	return Orig( guesses );

//	return FromLog( guesses, flog, acr.id, bcr.id );

/* ----------------------- */
/* Estimate starting angle */
/* ----------------------- */

	int	nPriorAngles = SetStartingAngle( flog );

/* ----------------------- */
/* Create image thumbnails */
/* ----------------------- */

	OlapRec	olp;
	ThmRec	thm;

	if( !SelectSubimage( olp, px, acr, bcr, flog ) )
		return false;

	if( !MakeThumbs( thm, olp, GBL.mch.THMDEC, flog ) )
		return false;

/* --------------------- */
/* Debug the angle sweep */
/* --------------------- */

// -----------------------------------------------------------
//	S.DebugAngs( GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile,
//		ang0, GBL.ctx.HFANGPR, .1, thm );
//	S.DebugAngs( GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile
//		ang0, 1, .01, thm );
//	exit( 42 );
// -----------------------------------------------------------

	if( dbgCor ) {
		S.RFromAngle( best, ang0, thm );
		return false;
	}

/* ------------------------------------------- */
/* Search range of angles for best correlation */
/* ------------------------------------------- */

	if( nPriorAngles ) {

		if( !S.UsePriorAngles( best, nPriorAngles, ang0,
				GBL.ctx.HFANGPR, thm ) ) {

			return Failure( best, acr.id, bcr.id, S.GetErr() );
		}
	}
	else if( !S.DenovoBestAngle( best,
				ang0, GBL.ctx.HFANGDN, 0.5, thm ) ) {

		return Failure( best, acr.id, bcr.id, S.GetErr() );
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

	{
		TForm	T, Tinv, I;

		AToBTrans( T, GBL.A.t2i.T, GBL.B.t2i.T );
		fprintf( flog, "Approx: Orig transform " );
		T.PrintTransform( flog );

		InvertTrans( Tinv, T );
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

	RecordResult( best, acr.id, bcr.id, S.GetErr() );

	guesses.push_back( best.T );

	return true;
}


