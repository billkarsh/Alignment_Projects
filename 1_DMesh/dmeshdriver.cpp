

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
	bool		nofile, tspose, force1;

/* --------------- */
/* Initialize data */
/* --------------- */

	printf( "\n---- Foldmaps ----\n" );

	wf		= px.wf;
	hf		= px.hf;
	Ntrans	= 0;
	rmap	= (uint16*)malloc( wf * hf * sizeof(uint16) );

	nofile	= GBL.arg.NoFolds || !GBL.msh.FLD;
	tspose	= GBL.arg.Transpose;
	force1	= GBL.arg.SingleFold;

	// Note that the foldmasks are always at full resolution.

	fold_mask_a = GetFoldMask(
					GBL.A.layer, GBL.A.tile,
					wf, hf, nofile, tspose, force1 );

	fold_mask_b = GetFoldMask(
					GBL.B.layer, GBL.B.tile,
					wf, hf, nofile, tspose, force1 );

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

	if( GBL.msh.WMT ) {

		sprintf( sfile, "../%d/%d/%d.%d.map.tif",
		GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile );

		Raster16ToTif8( sfile, rmap, wf, hf );
	}

	printf( "main: Got %d mapping regions.\n", Ntrans );

/* ------------------------------- */
/* Convert tform arrays to objects */
/* ------------------------------- */

	tfs = new TForm[Ntrans];
	ifs = new TForm[Ntrans];

	if( GBL.msh.WTT ) {

		sprintf( sfile, "../%d/%d/%d.%d.tf.txt",
		GBL.A.layer, GBL.A.tile, GBL.B.layer, GBL.B.tile );

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
		printf( "\n" );
	}

	if( f )
		fclose( f );

	if( tr_array )
		free( tr_array );
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

	if( !px.Load( GBL.A.t2i.path, GBL.B.t2i.path,
			GBL.msh.pxp.PXBRO, GBL.msh.pxp.PXDOG,
			GBL.msh.pxp.PXDOG_R1, GBL.msh.pxp.PXDOG_R2,
			stdout ) ) {

		goto exit;
	}

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

/* ----------- */
/* Diagnostics */
/* ----------- */

	if( GBL.arg.Verbose ) {

		ABOverlay( px, rmap, Ntrans, tfs, ifs );

		if( GBL.arg.Heatmap )
			RunCorrView( px, rmap, tfs, stdout );
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


