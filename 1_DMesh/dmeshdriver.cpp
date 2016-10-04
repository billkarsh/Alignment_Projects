

#include	"CGBL_dmesh.h"
#include	"FoldMask.h"
#include	"dmesh.h"
#include	"InSectionOverlap.h"

#include	"ImageIO.h"
#include	"Inspect.h"
#include	"Timer.h"
#include	"Memory.h"

#include	<stdlib.h>


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
    TAffine*		&tfs,
    TAffine*		&ifs,
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

    fprintf( stderr, "\n---- Foldmaps ----\n" );

    wf		= px.wf;
    hf		= px.hf;
    Ntrans	= 0;
    rmap	= (uint16*)malloc( wf * hf * sizeof(uint16) );

    // Note that the foldmasks are always at full resolution.

    {
        CCropMask	CM, *pCM = &CM;

        if( !CM.ReadIDB( GBL.idb, stderr ) )
            pCM = NULL;

        fold_mask_a = GetFoldMask(
                        GBL.idb, GBL.A, GBL.arg.fma,
                        px.resmska, pCM,
                        wf, hf, (GBL.ctx.FLD == 'N'),
                        GBL.arg.Transpose, GBL.arg.SingleFold,
                        stderr );

        fold_mask_b = GetFoldMask(
                        GBL.idb, GBL.B, GBL.arg.fmb,
                        px.resmskb, pCM,
                        wf, hf, (GBL.ctx.FLD == 'N'),
                        GBL.arg.Transpose, GBL.arg.SingleFold,
                        stderr );
    }

/* ------------- */
/* Call Pipeline */
/* ------------- */

    t0 = StartTiming();

    PipelineDeformableMap(
        Ntrans, tr_array, rmap,
        px, fold_mask_a, fold_mask_b, stderr );

    StopTiming( stderr, "Alignment", t0 );

    RasterFree( fold_mask_a );
    RasterFree( fold_mask_b );

/* -------------- */
/* Report results */
/* -------------- */

    if( GBL.mch.WMT ) {

        sprintf( sfile, "../%d/%d.%d.map.tif",
        GBL.A.id, GBL.B.z, GBL.B.id );

        Raster16ToTif8( sfile, rmap, wf, hf );
    }

    fprintf( stderr, "\nmain: Got %d mapping regions.\n", Ntrans );

/* ------------------------------- */
/* Convert tform arrays to objects */
/* ------------------------------- */

    tfs = new TAffine[Ntrans];
    ifs = new TAffine[Ntrans];

    if( GBL.mch.WTT ) {

        sprintf( sfile, "../%d/%d.%d.tf.txt",
        GBL.A.id, GBL.B.z, GBL.B.id );

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
        ifs[i].InverseOf( tfs[i] );

        fprintf( stderr, "main: Transform %3d: ", i );
        tfs[i].TPrint( stderr );
    }

    fprintf( stderr, "\n" );

    if( f )
        fclose( f );

    if( tr_array )
        free( tr_array );
}

/* --------------------------------------------------------------- */
/* Decomp -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Decomp( const TAffine &T, const char *label )
{
    fprintf( stderr, "main: %s: ", label );
    T.TPrint( stderr );

    double	r = T.GetRadians();
    TAffine	R, D;

    R.NUSetRot( -r );
    D = R * T;

    fprintf(  stderr, "main: Degrees: %g\n", r*180/PI );
    D.TPrint( stderr, "main: Residue: " );

    fprintf( stderr, "\n" );
}

/* --------------------------------------------------------------- */
/* ReportAveTAffine ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReportAveTAffine( int Ntrans, const TAffine* tfs )
{
    if( !Ntrans )
        return;

    TAffine	I, T = tfs[0];

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
    I.InverseOf( T );
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

    PixPair		px;
    uint16*		rmap	= NULL;
    TAffine*	tfs		= NULL;
    TAffine*	ifs		= NULL;
    int			Ntrans	= 0;

    if( !px.Load(
            GBL.A, GBL.B, GBL.idb,
            GBL.mch.PXLENS, GBL.mch.PXRESMSK, GBL.mch.PXBRO,
            GBL.mch.PXDOG, GBL.mch.PXDOG_R1, GBL.mch.PXDOG_R2,
            stderr, GBL.arg.Transpose ) ) {

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

        InSectionOverlap( Npts, apts, bpts, px, stderr );

        StopTiming( stderr, "InSectionOverlap", t1 );

        fprintf( stderr,
        "main: InSectionOverlap returned %d points.\n", Npts );

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

    ReportAveTAffine( Ntrans, tfs );

/* ----------- */
/* Diagnostics */
/* ----------- */

    if( Ntrans && GBL.arg.Verbose ) {

        ABOverlay( px, rmap, Ntrans, tfs, ifs,
            GBL.arg.comp_png, stderr );

        RunCorrView( px, rmap, tfs, GBL.arg.Heatmap,
            GBL.arg.registered_png, stderr );
    }

/* ------- */
/* Cleanup */
/* ------- */

    fprintf( stderr, "main: Normal completion for dmesh run.\n" );

    if( ifs )
        delete [] ifs;

    if( tfs )
        delete [] tfs;

    if( rmap )
        free( rmap );

exit:
    StopTiming( stderr, "Total", t0 );
    VMStats( stderr );

    return 0;
}


