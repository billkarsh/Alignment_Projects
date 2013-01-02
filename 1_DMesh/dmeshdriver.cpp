

#include	"CGBL_dmesh.h"
#include	"FoldMask.h"
#include	"dmesh.h"
#include	"InSectionOverlap.h"

#include	"ImageIO.h"
#include	"Inspect.h"
#include	"Timer.h"
#include	"Memory.h"


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* CalcTransforms ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Calculate/allocate a set of transforms tfs (and inverses ifs),
// each mapping a connected region of points in image A to image B.
//
// Also allocate uint16 image rmap of A's connected regions
// whose pixel values on exit are:
//
//    < 10:	no mapping
//   >= 10: mapped by tfs[pix - 10].
//
// Caller must release memory using:
//
//	delete [] tfs;
//	delete [] ifs;
//	free( rmap );
//
static void CalcTransforms(
	uint16*			&rmap,
	int				&Ntrans,
	TForm*			&tfs,
	TForm*			&ifs,
	const PixPair	&px )
{
	char		sfile[256];
	FILE*		f = NULL;
	uint8		*fold_mask_a,
				*fold_mask_b;
	double*		tr_array = NULL;
	clock_t		t0;
	int			wf, hf;

/* --------------- */
/* Initialize data */
/* --------------- */

	printf( "\n---- Foldmaps ----\n" );

	wf		= px.wf;
	hf		= px.hf;
	Ntrans	= 0;
	rmap	= (uint16*)malloc( wf * hf * sizeof(uint16) );

	// Note that the foldmasks are always at full resolution.

	fold_mask_a = GetFoldMask(
					GBL.idb, GBL.A.layer, GBL.A.tile,
					GBL.arg.fma, wf, hf, (GBL.ctx.FLD == 'N'),
					GBL.arg.Transpose, GBL.arg.SingleFold );

	fold_mask_b = GetFoldMask(
					GBL.idb, GBL.B.layer, GBL.B.tile,
					GBL.arg.fmb, wf, hf, (GBL.ctx.FLD == 'N'),
					GBL.arg.Transpose, GBL.arg.SingleFold );

/* ------------- */
/* Call Pipeline */
/* ------------- */

	t0 = StartTiming();

	PipelineDeformableMap(
		Ntrans, tr_array, rmap,
		px, fold_mask_a, fold_mask_b, stdout );

	StopTiming( stdout, "Alignment", t0 );

	RasterFree( fold_mask_a );
	RasterFree( fold_mask_b );

/* -------------- */
/* Report results */
/* -------------- */

	if( GBL.mch.WMT ) {

		sprintf( sfile, "../%d/%d.%d.map.tif",
		GBL.A.tile, GBL.B.layer, GBL.B.tile );

		Raster16ToTif8( sfile, rmap, wf, hf );
	}

	printf( "\nmain: Got %d mapping regions.\n", Ntrans );

/* ------------------------------- */
/* Convert tform arrays to objects */
/* ------------------------------- */

	tfs = new TForm[Ntrans];
	ifs = new TForm[Ntrans];

	if( GBL.mch.WTT ) {

		sprintf( sfile, "../%d/%d.%d.tf.txt",
		GBL.A.tile, GBL.B.layer, GBL.B.tile );

		f = fopen( sfile, "w" );
	}

	for( int i = 0; i < Ntrans; ++i ) {

		// copy-in Matlab-style values
		for( int j = 0; j < 6; ++j )
			tfs[i].t[j] = tr_array[i*6+j];

		// print in Matlab format
		if( f ) {
			fprintf( f, "%9.6f %9.6f %9.6f %9.6f %10.2f %10.2f\n",
			tfs[i].t[0], tfs[i].t[1], tfs[i].t[2],
			tfs[i].t[3], tfs[i].t[4], tfs[i].t[5] );
		}

		// now convert to our format
		tfs[i].FromMatlab();
		InvertTrans( ifs[i], tfs[i] );

		printf( "main: Transform %3d:", i );
		tfs[i].PrintTransform();
	}

	printf( "\n" );

	if( f )
		fclose( f );

	if( tr_array )
		free( tr_array );
}

/* --------------------------------------------------------------- */
/* Decomp -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Decomp( const TForm &T, const char *label )
{
	printf( "main: %s:", label );
	T.PrintTransform();

	double	r = RadiansFromAffine( T );
	TForm	R, D;

	R.NUSetRot( -r );
	MultiplyTrans( D, R, T );

	printf( "main: Degrees: %g\n", r*180/PI );
	printf( "main: Residue:" );
	D.PrintTransform();

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* ReportAveTForm ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ReportAveTForm( int Ntrans, const TForm* tfs )
{
	if( !Ntrans )
		return;

	TForm	I, T = tfs[0];

	if( Ntrans > 1 ) {

		for( int i = 1; i < Ntrans; ++i ) {

			for( int j = 0; j < 6; ++j )
				T.t[j] += tfs[i].t[j];
		}

		for( int j = 0; j < 6; ++j )
			T.t[j] /= Ntrans;
	}

	T.t[2] = T.t[5] = 0.0;

	Decomp( T, "Average" );
	InvertTrans( I, T );
	Decomp( I, "Inverse" );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	clock_t	t0 = StartTiming();

/* ------------------ */
/* Parse command line */
/* ------------------ */

	if( !GBL.SetCmdLine( argc, argv ) )
		return 42;

/* ---------- */
/* Get images */
/* ---------- */

	PixPair	px;
	uint16*	rmap	= NULL;
	TForm*	tfs		= NULL;
	TForm*	ifs		= NULL;
	int		Ntrans	= 0;

	if( !px.Load(
			GBL.A.t2i.path.c_str(),
			GBL.B.t2i.path.c_str(),
			(GBL.mch.PXLENS ? &GBL.idb : NULL),
			GBL.mch.PXBRO, GBL.mch.PXDOG,
			GBL.mch.PXDOG_R1, GBL.mch.PXDOG_R2,
			stdout, GBL.arg.Transpose ) ) {

		goto exit;
	}

/* ------------------- */
/* Scaling adjustments */
/* ------------------- */

	GBL.ctx.OLAP1D	/=  px.scl;
	GBL.ctx.OLAP2D	/= (px.scl * px.scl);
	GBL.mch.MNL		/=  px.scl;
	GBL.mch.MTA		/= (px.scl * px.scl);
	GBL.mch.MMA		/= (px.scl * px.scl);

/* ----------------------- */
/* Just test overlap code? */
/* ----------------------- */

	if( GBL.arg.WithinSection ) {

		double	*apts = NULL;
		double	*bpts = NULL;
		int		Npts;

		clock_t	t1 = StartTiming();

		InSectionOverlap( Npts, apts, bpts, px, stdout );

		StopTiming( stdout, "InSectionOverlap", t1 );

		printf( "main: InSectionOverlap returned %d points.\n", Npts );

		if( apts )
			free( apts );

		if( bpts )
			free( bpts );

		goto exit;
	}

/* ------------- */
/* Call Pipeline */
/* ------------- */

	CalcTransforms( rmap, Ntrans, tfs, ifs, px );

	ReportAveTForm( Ntrans, tfs );

/* ----------- */
/* Diagnostics */
/* ----------- */

	if( GBL.arg.Verbose ) {

		ABOverlay( px, rmap, Ntrans, tfs, ifs );

		RunCorrView( px, rmap, tfs, GBL.arg.Heatmap );
	}

/* ------- */
/* Cleanup */
/* ------- */

	printf( "main: Normal completion for dmesh run.\n" );

	if( ifs )
		delete [] ifs;

	if( tfs )
		delete [] tfs;

	if( rmap )
		free( rmap );

exit:
	StopTiming( stdout, "Total", t0 );
	VMStats( stdout );

	return 0;
}


