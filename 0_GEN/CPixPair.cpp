

#include	"CPixPair.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"CAffineLens.h"
#include	"Correlation.h"
#include	"Timer.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	MAX1DPIX	2048






/* --------------------------------------------------------------- */
/* MakeDoGKernel ------------------------------------------------- */
/* --------------------------------------------------------------- */

// This is an edge enhancement filter. Each Gaussian is a blurring
// filter that passes frequencies longer than its (inverse) standard
// dev (here, specified as a radius). The difference, then, passes
// frequencies between the two (inverse) radii: r1 < r2. The length
// scale of interest is the number of pixels over which there is
// rapid intensity variation, usually, near object edges. This is
// usually a few pixels, so choose {r1, r2} around {3, 6}.
//
// Note: Radius of kernel: 3 sigma = 3 x r2.
//
static int MakeDoGKernel(
    vector<double>	&DoG,
    int				r1,
    int				r2,
    FILE*			flog )
{
    int	d = 3 * r2,
        D = 2 * d + 1;

// Normalization sums should ideally equal 2pi. They will differ
// slightly due to quantization and the finite tail length.

    double	sum1 = 0.0, sum2 = 0.0;

    for( int y = -d; y <= d; ++y ) {

        for( int x = -d; x <= d; ++x ) {

            double	R = x*x + y*y;

            sum1 += exp( -( R/(2*r1*r1) ) );
            sum2 += exp( -( R/(2*r2*r2) ) );
        }
    }

    fprintf( flog,
    "MakeDoGKernel: Check %f equals %f equals 2pi.\n",
    sum1/(r1*r1), sum2/(r2*r2) );

// Store values

    DoG.resize( D * D );

    for( int y = -d; y <= d; ++y ) {

        for( int x = -d; x <= d; ++x ) {

            double	R	= x*x + y*y;
            double	v1	= exp( -( R/(2*r1*r1) ) );
            double	v2	= exp( -( R/(2*r2*r2) ) );

            DoG[x+d + D*(y+d)] = v1/sum1 - v2/sum2;
        }
    }

    return D;
}

/* --------------------------------------------------------------- */
/* PixPair::Downsample ------------------------------------------- */
/* --------------------------------------------------------------- */

void PixPair::Downsample(
    vector<double>			&dst,
    const vector<double>	&src )
{
    int	n = scl * scl;

    dst.resize( ws * hs );

    for( int iy = 0; iy < hs; ++iy ) {

        for( int ix = 0; ix < ws; ++ix ) {

            double	sum = 0.0;

            for( int dy = 0; dy < scl; ++dy ) {

                for( int dx = 0; dx < scl; ++dx )
                    sum += src[ix*scl+dx + wf*(iy*scl+dy)];
            }

            dst[ix+ws*iy] = sum / n;
        }
    }
}

/* --------------------------------------------------------------- */
/* Trim ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Trim Nathan images and adjust w, h.
//
static void Trim( uint8* ras, uint32 &w, uint32 &h )
{
// Trim at least top two rows and make row count even

//	int	htrim = (h & 1 ? 3 : 2);	// Nathan
    int	htrim = (h & 1 ? 1 : 0);	// Davi (just make even)

    if( htrim )
        memmove( ras, ras + htrim*w, w*(h-=htrim) );

// Trim rightmost col to make width even

    if( w & 1 ) {

        uint8	*dst, *src;

        dst = ras + w - 1;
        src = ras + w;

        for( int i = 1; i < h; ++i, dst += w - 1, src += w )
            memmove( dst, src, w - 1 );

        w -= 1;
    }
}

/* --------------------------------------------------------------- */
/* HasTissue ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to kick out mostly agar images in Nathan data.
//
static bool HasTissue( const char *path, FILE *flog )
{
    uint32	w, h;
    uint16	*ras = Raster16FromTif16( path, w, h, flog );
    int		N = w*h, nlow = 0;

    for( int i = 0; i < N; ++i ) {
        if( ras[i] < 4000 )
            ++nlow;
    }

    RasterFree( ras );

    return nlow < 0.20 * N;
}

/* --------------------------------------------------------------- */
/* Lens ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Lens(
    vector<double>	&vout,
    CAffineLens		&LN,
    const uint8		*ras,
    int				w,
    int				h,
    int				order,
    int				cam )
{
// Flatten and release raster

    vector<double>	vflat;
    LegPolyFlatten( vflat, ras, w, h, order );

// Transform into vout

    TAffine	T = LN.GetTf( cam );
    int		np = w * h;

    vout.resize( np, 0.0 );

    for( int i = 0; i < np; ++i ) {

        int		y = i / w,
                x = i - w*y;
        Point	p( x, y );

        T.Transform( p );
        DistributePixel( p.x, p.y, vflat[i], vout, w, h );
    }
}

/* --------------------------------------------------------------- */
/* PixPair::Load ------------------------------------------------- */
/* --------------------------------------------------------------- */

bool PixPair::Load(
    const PicSpec	&A,
    const PicSpec	&B,
    const string	&idb,
    bool			lens,
    bool			resmsk,
    int				order,
    int				bDoG,
    int				r1,
    int				r2,
    FILE*			flog,
    bool			transpose )
{
    fprintf( flog, "\n---- Image loading ----\n" );

/* ----------------------------- */
/* Load and sanity check rasters */
/* ----------------------------- */

    clock_t	t0 = StartTiming();
    uint8	*aras, *bras;
    uint32	wa, ha, wb, hb;
    int		ok = false;

    aras = Raster8FromAny( A.t2i.path.c_str(),
            wa, ha, flog, transpose );

    bras = Raster8FromAny( B.t2i.path.c_str(),
            wb, hb, flog, transpose );

    if( !aras || !bras ) {
        fprintf( flog,
        "FAIL: PixPair: Picture load failure.\n" );
        goto exit;
    }

    if( wa != wb || ha != hb ) {
        fprintf( flog,
        "FAIL: PixPair: Nonmatching picture dimensions.\n" );
        goto exit;
    }

//-----------------------------------------------------------
// Experiment to filter out mostly-agar Nathan images (not a good
// way to do this because uses absolute intensity assumptions).
// Also we trim off some rows and columns to fix aperture and
// odd sizing issues (if needed).
//
    //if( !HasTissue( A.t2i.path.c_str(), flog ) ||
    //	!HasTissue( B.t2i.path.c_str(), flog ) ) {

    //	goto exit;
    //}
    //Trim( aras, wa, ha );
    //Trim( bras, wb, hb );
//-----------------------------------------------------------

//-----------------------------------------------------------
// Quick fix if y-dim not multiple of two.
//
    //if( ha & 1 ) --ha, --hb;
//-----------------------------------------------------------

    wf		= wa;
    hf		= ha;
    ws		= wa;
    hs		= ha;
    scl		= 1;

/* ------- */
/* Flatten */
/* ------- */

    if( lens ) {

        CAffineLens	LN;

        if( !LN.ReadIDB( idb, flog ) )
            goto exit;

        Lens( _avf, LN, aras, wf, hf, order, A.t2i.cam );
        Lens( _bvf, LN, bras, wf, hf, order, B.t2i.cam );
        //VectorDblToTif8( "LensA.tif", _avf, wf, hf, flog );
        //VectorDblToTif8( "LensB.tif", _bvf, wf, hf, flog );
    }
    else {

        LegPolyFlatten( _avf, aras, wf, hf, order );
        LegPolyFlatten( _bvf, bras, wf, hf, order );
    }

/* ------------- */
/* Resin masking */
/* ------------- */

    if( resmsk ) {

        double	tisfraca, tisfracb;
        int		n = wf * hf, suma = 0, sumb = 0;

        // first make smoothest masks
        ResinMask8( resmska, aras, wf, hf, false );
        ResinMask8( resmskb, bras, wf, hf, false );

        // reject if no tissue
        for( int i = 0; i < n; ++i ) {
            suma += resmska[i];
            sumb += resmskb[i];
        }

        tisfraca = (double)suma / n;
        tisfracb = (double)sumb / n;

        fprintf( flog, "Tissue frac: A %.3f B %.3f\n",
        tisfraca, tisfracb );

        if( tisfraca < 0.15 && tisfracb < 0.15 ) {
            fprintf( flog,
            "FAIL: PixPair: Low tissue fraction [%.3f %.3f]\n",
            tisfraca, tisfracb );
            goto exit;
        }

        // remake masks if same layer
        if( A.z == B.z ) {
            ResinMask8( resmska, aras, wf, hf, true );
            ResinMask8( resmskb, bras, wf, hf, true );
        }

        //Raster8ToTif8( "resinA.tif", &resmska[0], wf, hf, flog );
        //Raster8ToTif8( "resinB.tif", &resmskb[0], wf, hf, flog );

#if 0
// This section moves resin pixel values toward zero
// but should'nt be necessary if removing resin from
// point lists via fold mask machinery.

        for( int i = 0; i < n; ++i ) {

            if( !resmska[i] )
                _avf[i] = 0.0;

            if( !resmskb[i] )
                _bvf[i] = 0.0;
        }

        Normalize( _avf );
        Normalize( _bvf );
#endif
    }

/* ------------- */
/* Init pointers */
/* ------------- */

    avs_vfy	= avs_aln = avf_vfy	= avf_aln = &_avf;
    bvs_vfy	= bvs_aln = bvf_vfy	= bvf_aln = &_bvf;

/* ------------- */
/* Apply filters */
/* ------------- */

    if( bDoG ) {

        vector<double>	DoG;
        vector<CD>		kfft;
        int				dim = MakeDoGKernel( DoG, r1, r2, flog );

        Convolve( _avfflt, _avf, wf, hf,
            &DoG[0], dim, dim, true, true, kfft, flog );
        Normalize( _avfflt );

        Convolve( _bvfflt, _bvf, wf, hf,
            &DoG[0], dim, dim, true, true, kfft, flog );
        Normalize( _bvfflt );

        avs_aln = avf_aln = &_avfflt;
        bvs_aln = bvf_aln = &_bvfflt;
    }

//-----------------------------------------------------------
// Experiment to apply spatial averaging to Nathan data...
// Not needed if Nathan averages several (10) layers in z.
//
    //VectorDblToTif8( "PRE.tif", *avs_aln, ws, hs, flog );
    //{
    //	double		K[] = {
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    //	vector<CD>	kfft;
    //	int			dm = 7;

    //	Convolve( _avfflt, _avf, wf, hf, K, dm, dm, true, true, kfft, flog );
    //	Normalize( _avfflt );

    //	Convolve( _bvfflt, _bvf, wf, hf, K, dm, dm, true, true, kfft, flog );
    //	Normalize( _bvfflt );

    //	avs_aln = avf_aln = &_avfflt;
    //	bvs_aln = bvf_aln = &_bvfflt;

    //	bDoG = 1;
    //}
    //VectorDblToTif8( "POST.tif", *avs_aln, ws, hs, flog );
//-----------------------------------------------------------

/* --------------------- */
/* Downsample all images */
/* --------------------- */

    if( ws > MAX1DPIX || hs > MAX1DPIX ) {

        do {
            ws		/= 2;
            hs		/= 2;
            scl		*= 2;
        } while( ws > MAX1DPIX || hs > MAX1DPIX );

        fprintf( flog, "PixPair: Scaling by %d\n", scl );

        if( ws * scl != wf || hs * scl != hf ) {

            fprintf( flog,
            "FAIL: PixPair: Dimensions not multiple of scale!\n" );
            goto exit;
        }

        Downsample( _avs, _avf );
        Downsample( _bvs, _bvf );

        avs_vfy = avs_aln = &_avs;
        bvs_vfy = bvs_aln = &_bvs;

        if( bDoG ) {

            if( _avfflt.size() ) {
                Downsample( _avsflt, _avfflt );
                avs_aln = &_avsflt;
            }

            if( _bvfflt.size() ) {
                Downsample( _bvsflt, _bvfflt );
                bvs_aln = &_bvsflt;
            }
        }
    }
    else
        fprintf( flog, "PixPair: Using image scale=1.\n" );

/* ------------------------------ */
/* Write DoG images for debugging */
/* ------------------------------ */

#if 0
    if( bDoG ) {
        VectorDblToTif8( "DoGa.tif", *avs_aln, ws, hs, flog );
        VectorDblToTif8( "DoGb.tif", *bvs_aln, ws, hs, flog );
    }
#endif

    ok = true;

/* -------- */
/* Clean up */
/* -------- */

exit:
    if( aras )
        RasterFree( aras );

    if( bras )
        RasterFree( bras );

    StopTiming( flog, "Image conditioning", t0 );

    return ok;
}


