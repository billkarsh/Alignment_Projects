

#include	"ImageIO.h"
#include	"File.h"
#include	"mrc.h"
#include	"Maths.h"

#include	"png.h"
#include	"tiffio.h"

#include	<limits.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	USE_TIF_DEFLATE		0
#define	GENEMEYERSTIFF		1

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// write out converted images as 1.tif, 2.tif, etc.
static int num = 1;






/* --------------------------------------------------------------- */
/* RasterAlloc --------------------------------------------------- */
/* --------------------------------------------------------------- */

void* RasterAlloc( unsigned long bytes )
{
    return _TIFFmalloc( bytes );
}

/* --------------------------------------------------------------- */
/* RasterFree ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void _RasterFree( void** praster )
{
    if( *praster ) {
        _TIFFfree( *praster );
        *praster = NULL;
    }
}

/* --------------------------------------------------------------- */
/* Raster8FromAny ------------------------------------------------ */
/* --------------------------------------------------------------- */

uint8* Raster8FromAny(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog,
    bool		transpose )
{
    const char*	p;

    p = strstr( name, ".tif" );
    if( p && !p[4] )
        return Raster8FromTif( name, w, h, flog, transpose );

    p = strstr( name, ".mrc" );
    if( p && !p[4] )
        return ReadAnMRCFile( name, w, h, flog, transpose );

    p = strstr( name, ".png" );
    if( p && !p[4] )
        return Raster8FromPng( name, w, h, flog, transpose );

// unknown format
    fprintf( flog,
    "Raster8FromAny: Unknown file suffix for [%s].\n", name );
    exit( 42 );
}

/* --------------------------------------------------------------- */
/* Raster8FromTifFloat ------------------------------------------- */
/* --------------------------------------------------------------- */

// Some MRC images are written with 32-bit floating point samples.
// Convert them to 8-bit integer samples.
//
// Remap to {u=127.5, sd=25}.
//
// (This is inefficient since we later compute floating values again,
// but it makes the code easier.)
//
static void Raster8FromTifFloat(
    TIFF	*tif,
    int		w,
    int		h,
    uint8*	raster,
    FILE*	flog,
    bool	writeDebug = false )
{
    vector<double>	v( w*h );
    tdata_t			raw;
    float			*b2;
    uint32			row;
    int				sls, j;

// float -> doubles vector
    sls = TIFFScanlineSize( tif );
    raw = malloc( sls );

    fprintf( flog, "TIF(f): Image length is %d\n", h );
    fprintf( flog, "TIF(f): Scan line size is %d bytes.\n", sls );

    for( row = 0; row < h; ++row ) {

        TIFFReadScanline( tif, raw, row );
        b2 = (float*)raw;

        for( j = 0; j < w; ++j )
            v[j + w * row] = b2[j];
    }

    free( raw );

// Normalize vector
    Normalize( v );

// Copy vector -> raster.
// Change mean to 127, std dev to 25

    int	np = w * h;

    for( j = 0; j < np; ++j ) {

        int pix	= int(127.5 + v[j]*25.0);

        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;

        raster[j] = pix;
    }

// Write in GIMP format?
    if( writeDebug ) {

        char	name[32];

        sprintf( name, "%d.tif", num );

        fprintf( flog,
        "TIF(f): ----> writing as %d.tif, w=%d h=%d\n",
        num, w, h );

        Raster8ToTif8( name, raster, w, h, flog );
    }

// Always advance name-index whether write or not
    ++num;
}

/* --------------------------------------------------------------- */
/* Raster8FromTif16Bit ------------------------------------------- */
/* --------------------------------------------------------------- */

// Some MRC images are written as 16-bit TIFFs.
// Convert them to 8-bit integer samples.
//
// Remap to {u=127.5, sd=25}.
//
static void Raster8FromTif16Bit(
    TIFF	*tif,
    int		w,
    int		h,
    uint8*	raster,
    FILE*	flog,
    bool	writeDebug = false )
{
    int				npixels = w * h;
    uint16			*raw;
    vector<double>	v( npixels );
    int				NB, nb, j;
    uint16			fmt;

// Check for signed data if 16-bit
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLEFORMAT, &fmt );

// Read encoded strips
    raw = (uint16*)malloc( npixels * sizeof(uint16) );

    for( j = 0, NB = 0; NB < npixels * 2; NB += nb, ++j ) {

        nb = TIFFReadEncodedStrip( tif, j,
                raw + NB / 2, (tsize_t)-1 );
    }

    fprintf( flog,
    "TIF(16): Last Encoded strip had %d bytes.\n", nb );

// Pixels -> doubles vector

    if( fmt == SAMPLEFORMAT_INT ) {

        for( j = 0; j < npixels; ++j )
            v[j] = ((int16*)raw)[j] - (int)SHRT_MIN;
    }
    else {

        for( j = 0; j < npixels; ++j )
            v[j] = raw[j];
    }

    free( raw );

// Normalize vector
    Normalize( v );

// Copy vector -> raster.
// Apply optional background subtraction
// Change mean to 127, std dev to 25 (or other as specified)

    const char	*p	= getenv( "Convert16BitBkgSub" );
    double		bkg	= (p == NULL ? 0.0 : atof( p ));
                 p  = getenv( "Convert16BitStdDev" );
    double		std	= (p == NULL ? 25.0 : atof( p )),
                mn	= 127.5 - bkg;

    fprintf( flog,
    "TIF(16): Converting intensity using bkg %f, stddev %f\n",
    bkg, std );

    for( j = 0; j < npixels; ++j ) {

        int pix	= int(mn + v[j]*std);

        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;

        raster[j] = pix;
    }

// Write in GIMP format?
    if( writeDebug ) {

        char	name[32];

        sprintf( name, "%d.tif", num );

        fprintf( flog,
        "TIF(16):  ----> writing as %d.tif, w=%d h=%d\n",
        num, w, h );

        Raster8ToTif8( name, raster, w, h, flog );
    }

// Always advance name-index whether write or not
    ++num;
}

/* --------------------------------------------------------------- */
/* Raster8FromTif8Bit -------------------------------------------- */
/* --------------------------------------------------------------- */

// Read 8-bit TIFFs.
//
//
// No intensity remapping.
//
static void Raster8FromTif8Bit(
    TIFF	*tif,
    int		w,
    int		h,
    uint8*	raster,
    FILE*	flog )
{
    int		npixels = w * h;
    int		NB, nb, j;

// read encoded strips
    for( j = 0, NB = 0; NB < npixels; NB += nb, ++j ) {

        nb = TIFFReadEncodedStrip( tif, j,
                raster + NB, (tsize_t)-1 );
    }

    fprintf( flog,
    "TIF(8): Last Encoded strip had %d bytes.\n", nb );
}

/* --------------------------------------------------------------- */
/* Raster8FromTifRGBA -------------------------------------------- */
/* --------------------------------------------------------------- */

// Read RGBA TIFFs.
//
// Return only R-component.
//
// No intensity remapping.
//
static void Raster8FromTifRGBA(
    TIFF	*tif,
    int		w,
    int		h,
    uint8*	raster )
{
    int		npixels = w * h;
    uint32	*raw = (uint32*)malloc( npixels * sizeof(uint32) );

    TIFFReadRGBAImage( tif, w, h, raw, 0 );

    for( int i = 0; i < npixels; ++i )
        raster[i] = raw[i] & 0xFF;

    free( raw );
}

/* --------------------------------------------------------------- */
/* Transpose ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Transpose( uint8* raster, uint32 &w, uint32 &h )
{
// swap pixels

    int	np = w * h;

    for( int i = 0; i < np; ++i ) {

        int	t;
        int	y = i / w,
            x = i - w * y,
            k = y + h * x;

        t			= raster[k];
        raster[k]	= raster[i];
        raster[i]	= t;
    }

// swap w and h

    uint32	t;

    t = w;
    w = h;
    h = t;
}

/* --------------------------------------------------------------- */
/* FlipY --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void FlipY( uint8* raster, uint32 w, uint32 h )
{
// swap rows

    vector<uint8>	t( w );

    for( int y = 0; y < h/2; ++y ) {

        memcpy( &t[0], &raster[w*y], w );
        memcpy( &raster[w*y], &raster[w*(h-1-y)], w );
        memcpy( &raster[w*(h-1-y)], &t[0], w );
    }
}

/* --------------------------------------------------------------- */
/* Raster8FromTif ------------------------------------------------ */
/* --------------------------------------------------------------- */

uint8* Raster8FromTif(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog,
    bool		transpose )
{
    TIFF*	tif;
    uint8*	raster;
    size_t	npixels;
    uint16	bps, spp;

    fprintf( flog, "Raster8FromTif: Opening [%s].\n", name );

#if GENEMEYERSTIFF
// suppress exotic field warnings
TIFFErrorHandler	oldEH = TIFFSetWarningHandler( NULL );
#endif

    tif = TIFFOpen( name, "r" );

#if GENEMEYERSTIFF
TIFFSetWarningHandler( oldEH );
#endif

    if( !tif ) {
        fprintf( flog,
        "Raster8FromTif: Cannot open [%s] for read.\n", name );
        exit( 42 );
    }

    if( !TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &w ) ) {
        fprintf( flog, "Raster8FromTif: Missing image width tag in [%s].\n", name );
        exit( 42 );
    }

    if( !TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &h ) ) {
        fprintf( flog, "Raster8FromTif: Missing image length tag in [%s].\n", name );
        exit( 42 );
    }

    TIFFGetFieldDefaulted( tif, TIFFTAG_BITSPERSAMPLE, &bps );
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLESPERPIXEL, &spp );
    npixels = w * h;

    fprintf( flog,
    "Raster8FromTif: Bits per sample %d, samples per pixel %d\n",
    bps, spp );

    fprintf( flog,
    "Raster8FromTif: Picture is %d by %d, %ld pixels total.\n",
    w, h, npixels );

    raster = (uint8*)_TIFFmalloc( npixels * sizeof(uint8) );

    if( raster ) {

        if( spp == 1 || spp == 0 ) {

            if( bps == 32 )
                Raster8FromTifFloat( tif, w, h, raster, flog );
            else if( bps == 16 )
                Raster8FromTif16Bit( tif, w, h, raster, flog );
            else
                Raster8FromTif8Bit( tif, w, h, raster, flog );
        }
        else if( bps == 8 && spp == 4 )
            Raster8FromTifRGBA( tif, w, h, raster );
        else {
            fprintf( flog, "Raster8FromTif: Unknown format.\n" );
            exit( 42 );
        }
    }
    else {
        fprintf( flog, "Raster8FromTif: Malloc failed.\n" );
        exit( 42 );
    }

    TIFFClose( tif );

// Enable following line to flip image vertically
//	FlipY( raster, w, h );

// optional transposition
    if( transpose )
        Transpose( raster, w, h );

    return raster;
}

/* --------------------------------------------------------------- */
/* Raster16FromTif16 --------------------------------------------- */
/* --------------------------------------------------------------- */

uint16* Raster16FromTif16(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog )
{
    TIFF*	tif;
    uint16*	raster;
    size_t	npixels;
    uint16	bps, spp, fmt;

#if GENEMEYERSTIFF
// suppress exotic field warnings
TIFFErrorHandler	oldEH = TIFFSetWarningHandler( NULL );
#endif

    tif = TIFFOpen( name, "r" );

#if GENEMEYERSTIFF
TIFFSetWarningHandler( oldEH );
#endif

    if( !tif ) {
        fprintf( flog,
        "Raster16FromTif16: Cannot open [%s] for read.\n", name );
        exit( 42 );
    }

    if( !TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &w ) ) {
        fprintf( flog, "Raster16FromTif16: Missing image width tag in [%s].\n", name );
        exit( 42 );
    }

    if( !TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &h ) ) {
        fprintf( flog, "Raster16FromTif16: Missing image length tag in [%s].\n", name );
        exit( 42 );
    }

    TIFFGetFieldDefaulted( tif, TIFFTAG_BITSPERSAMPLE, &bps );
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLESPERPIXEL, &spp );
    TIFFGetFieldDefaulted( tif, TIFFTAG_SAMPLEFORMAT, &fmt );
    npixels = w * h;

    fprintf( flog,
    "Raster16FromTif16: Bits per sample %d, samples per pixel %d\n",
    bps, spp );

    fprintf( flog,
    "Raster16FromTif16: Picture is %d by %d, %ld pixels total.\n",
    w, h, npixels );

    raster = (uint16*)_TIFFmalloc( npixels * sizeof(uint16) );

    if( raster ) {

        if( spp == 1 && bps == 16 ) {

            int		NB, nb, j;

            // read encoded strips
            for( j = 0, NB = 0; NB < npixels * 2; NB += nb, ++j ) {

                nb = TIFFReadEncodedStrip( tif, j,
                        raster + NB / 2, (tsize_t)-1 );
            }

            if( fmt == SAMPLEFORMAT_INT ) {

                for( int i = 0; i < npixels; ++i )
                    raster[i] = ((int16*)raster)[i] - (int)SHRT_MIN;
            }
        }
        else {
            fprintf( flog, "Raster16FromTif16: Unknown format.\n" );
            exit( 42 );
        }
    }
    else {
        fprintf( flog, "Raster16FromTif16: Malloc failed.\n" );
        exit( 42 );
    }

    TIFFClose( tif );

    return raster;
}

/* --------------------------------------------------------------- */
/* Raster8ToTifFlt ----------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster8ToTifFlt(
    const char*		name,
    const uint8*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;
    float	*f;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "TIF(f) Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 32 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    f = (float*)malloc( w * sizeof(float) );

    for( int row = 0; row < h; ++row ) {

        for( int j = 0; j < w; ++j )
            f[j] = raster[j + w * row];

        TIFFWriteScanline( image, f, row, 0 );
    }

    free( f );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* Raster8ToTif8 ------------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster8ToTif8(
    const char*		name,
    const uint8*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "Tif(8) Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 8 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    TIFFWriteEncodedStrip( image, 0,
        (void*)raster, w * h * sizeof(uint8) );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* Raster16ToTif16 ----------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster16ToTif16(
    const char*		name,
    const uint16*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "Tif(16) Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 16 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    TIFFWriteEncodedStrip( image, 0,
        (void*)raster, w * h * sizeof(uint16) );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* Raster32ToTifRGBA --------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster32ToTifRGBA(
    const char*		name,
    const uint32*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "TIF(RGB): Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 8 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 4 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );
    TIFFSetField( image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    TIFFWriteEncodedStrip( image, 0,
        (void*)raster, w * h * sizeof(uint32) );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* Raster16ToTif8 ------------------------------------------------ */
/* --------------------------------------------------------------- */

void Raster16ToTif8(
    const char*		name,
    const uint16*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    int				N = w * h;
    vector<uint8>	buf( N );

    for( int i = 0; i < N; ++i )
        buf[i] = raster[i];

    Raster8ToTif8( name, &buf[0], w, h, flog );
}

/* --------------------------------------------------------------- */
/* RasterFltToTifFlt --------------------------------------------- */
/* --------------------------------------------------------------- */

void RasterFltToTifFlt(
    const char*		name,
    const float*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;
    float	*f;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "TIF(f) Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 32 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    for( int row = 0; row < h; ++row )
        TIFFWriteScanline( image, (void*)&raster[w*row], row, 0 );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* RasterDblToTifFlt --------------------------------------------- */
/* --------------------------------------------------------------- */

void RasterDblToTifFlt(
    const char*		name,
    const double*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    TIFF	*image;
    float	*f;

    if( !(image = TIFFOpen( name, "w" )) ) {
        fprintf( flog,
        "TIF(f) Could not open [%s] for writing.\n", name );
        exit( 42 );
    }

// Set values for basic tags before adding data
    TIFFSetField( image, TIFFTAG_IMAGEWIDTH, w );
    TIFFSetField( image, TIFFTAG_IMAGELENGTH, h );
    TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 32 );
    TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, h );

#if USE_TIF_DEFLATE
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
#else
    TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_NONE );
#endif

    TIFFSetField( image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
    TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    TIFFSetField( image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );

// Write the information to the file
    f = (float*)malloc( w * sizeof(float) );

    for( int row = 0; row < h; ++row ) {

        for( int j = 0; j < w; ++j )
            f[j] = raster[j + w * row];

        TIFFWriteScanline( image, f, row, 0 );
    }

    free( f );

// Close the file
    TIFFClose( image );
}

/* --------------------------------------------------------------- */
/* CorrThmToTif8 ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CorrThmToTif8(
    const char*				name,
    const vector<double>	&I,
    int						iw,
    int						tw,
    int						th,
    FILE*					flog )
{
// Change mean to 127, std dev to 35 (or other as specified)

    vector<uint8>	buf( tw * th );
    const char		*p	= getenv( "Convert16BitStdDev" );
    double			std	= (p == NULL ? 35.0 : atof( p ));

    for( int i = 0; i < th; ++i ) {

        for( int j = 0; j < tw; ++j ) {

            int	pix = 127 + (int)(I[j + iw*i]*std);

            if( pix < 0 )
                pix = 0;
            else if( pix > 255 )
                pix = 255;

            buf[j + tw*i] = pix;
        }
    }

    Raster8ToTif8( name, &buf[0], tw, th, flog );
}

/* --------------------------------------------------------------- */
/* VectorDblToTif8 ----------------------------------------------- */
/* --------------------------------------------------------------- */

void VectorDblToTif8(
    const char*				name,
    const vector<double>	&vals,
    int						w,
    int						h,
    FILE*					flog )
{
    int				nPts = vals.size();
    vector<uint8>	buf( nPts );

    for( int i = 0; i < nPts; ++i ) {

        int	pix = 127 + (int)(vals[i]*35.0);

        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;

        buf[i] = pix;
    }

    Raster8ToTif8( name, &buf[0], w, h, flog );
}

/* --------------------------------------------------------------- */
/* Raster8FromPng ------------------------------------------------ */
/* --------------------------------------------------------------- */

uint8* Raster8FromPng(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog,
    bool		transpose )
{
    FILE*		f			= NULL;
    png_structp	png_ptr		= NULL;
    png_infop	info_ptr	= NULL;
    uint8*		raster		= NULL;
    png_byte	header[8];
    png_uint_32	wi, hi;
    int			bit_depth, color_type, x, y, n, ok = false;

// open
    if( !(f = fopen( name, "rb" )) ) {
        fprintf( flog, "PNG(8): Cannot open [%s] for read.\n", name );
        goto exit;
    }

// verify
    fread( header, sizeof(header), 1, f );

    if( png_sig_cmp( header, 0, 8 ) ) {
        fprintf( flog, "PNG(8): Bad header [%s].\n", name );
        goto exit;
    }

// inits
    png_ptr = png_create_read_struct(
                PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    if( !png_ptr ) {
        fprintf( flog,
        "PNG(8): Failed read_struct alloc [%s].\n", name );
        goto exit;
    }

    info_ptr = png_create_info_struct( png_ptr );

    if( !info_ptr ) {
        fprintf( flog,
        "PNG(8): Failed info_struct alloc [%s].\n", name );
        goto exit;
    }

    png_init_io( png_ptr, f );
    png_set_sig_bytes( png_ptr, 8 );

// read
    png_read_png( png_ptr, info_ptr,
        PNG_TRANSFORM_PACKING | PNG_TRANSFORM_SWAP_ENDIAN,
        NULL );

// format
    png_get_IHDR( png_ptr, info_ptr,
        &wi, &hi, &bit_depth, &color_type, NULL, NULL, NULL );

    w = wi;
    h = hi;

// get rows
    raster = (uint8*)malloc( w * h * sizeof(uint8) );

    if( !raster ) {
        fprintf( flog, "PNG(8): Failed raster alloc [%s].\n", name );
        goto exit;
    }

    if( color_type == PNG_COLOR_TYPE_GRAY ) {

        if( bit_depth == 16 ) {

            png_uint_16**	row_ptrs = (png_uint_16**)
                            png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
        else {	// 1, 2, 4, 8

            png_byte**	row_ptrs = (png_byte**)
                        png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
    }
    else {

        // Note:
        // png_uint_32 is typedef'd as unsigned long
        // which is 64 bits on the cluster

        uint32**	row_ptrs = (uint32**)
                    png_get_rows( png_ptr, info_ptr );

        n = 0;

        for( y = 0; y < h; ++y ) {

            for( x = 0; x < w; ++x )
                raster[n++] = row_ptrs[y][x];
        }
    }

// report success
    ok = true;

// done
exit:
    if( png_ptr )
        png_destroy_read_struct( &png_ptr, &info_ptr, NULL );

    if( f )
        fclose( f );

    if( !ok ) {

        if( raster ) {
            free( raster );
            raster = NULL;
        }
    }
    else if( transpose )
        Transpose( raster, w, h );

    return raster;
}

/* --------------------------------------------------------------- */
/* Raster16FromPng ----------------------------------------------- */
/* --------------------------------------------------------------- */

uint16* Raster16FromPng(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog )
{
    FILE*		f			= NULL;
    png_structp	png_ptr		= NULL;
    png_infop	info_ptr	= NULL;
    uint16*		raster		= NULL;
    png_byte	header[8];
    png_uint_32	wi, hi;
    int			bit_depth, color_type, x, y, n, ok = false;

// open
    if( !(f = fopen( name, "rb" )) ) {
        fprintf( flog, "PNG(16): Cannot open [%s] for read.\n", name );
        goto exit;
    }

// verify
    fread( header, sizeof(header), 1, f );

    if( png_sig_cmp( header, 0, 8 ) ) {
        fprintf( flog, "PNG(16): Bad header [%s].\n", name );
        goto exit;
    }

// inits
    png_ptr = png_create_read_struct(
                PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    if( !png_ptr ) {
        fprintf( flog,
        "PNG(16): Failed read_struct alloc [%s].\n", name );
        goto exit;
    }

    info_ptr = png_create_info_struct( png_ptr );

    if( !info_ptr ) {
        fprintf( flog,
        "PNG(16): Failed info_struct alloc [%s].\n", name );
        goto exit;
    }

    png_init_io( png_ptr, f );
    png_set_sig_bytes( png_ptr, 8 );

// read
    png_read_png( png_ptr, info_ptr,
        PNG_TRANSFORM_PACKING | PNG_TRANSFORM_SWAP_ENDIAN,
        NULL );

// format
    png_get_IHDR( png_ptr, info_ptr,
        &wi, &hi, &bit_depth, &color_type, NULL, NULL, NULL );

    w = wi;
    h = hi;

// get rows
    raster = (uint16*)malloc( w * h * sizeof(uint16) );

    if( !raster ) {
        fprintf( flog, "PNG(16): Failed raster alloc [%s].\n", name );
        goto exit;
    }

    if( color_type == PNG_COLOR_TYPE_GRAY ) {

        if( bit_depth == 16 ) {

            png_uint_16**	row_ptrs = (png_uint_16**)
                            png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
        else {	// 1, 2, 4, 8

            png_byte**	row_ptrs = (png_byte**)
                        png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
    }
    else {

        // Note:
        // png_uint_32 is typedef'd as unsigned long
        // which is 64 bits on the cluster

        uint32**	row_ptrs = (uint32**)
                    png_get_rows( png_ptr, info_ptr );

        n = 0;

        for( y = 0; y < h; ++y ) {

            for( x = 0; x < w; ++x )
                raster[n++] = row_ptrs[y][x] & 0xFF;
        }
    }

// report success
    ok = true;

// done
exit:
    if( png_ptr )
        png_destroy_read_struct( &png_ptr, &info_ptr, NULL );

    if( f )
        fclose( f );

    if( !ok ) {

        if( raster ) {
            free( raster );
            raster = NULL;
        }
    }

    return raster;
}

/* --------------------------------------------------------------- */
/* Raster32FromPng ----------------------------------------------- */
/* --------------------------------------------------------------- */

uint32* Raster32FromPng(
    const char*	name,
    uint32		&w,
    uint32		&h,
    FILE*		flog )
{
    FILE*		f			= NULL;
    png_structp	png_ptr		= NULL;
    png_infop	info_ptr	= NULL;
    uint32*		raster		= NULL;
    png_byte	header[8];
    png_uint_32	wi, hi;
    int			bit_depth, color_type, x, y, n, ok = false;

// open
    if( !(f = fopen( name, "rb" )) ) {
        fprintf( flog, "PNG(32): Cannot open [%s] for read.\n", name );
        goto exit;
    }

// verify
    fread( header, sizeof(header), 1, f );

    if( png_sig_cmp( header, 0, 8 ) ) {
        fprintf( flog, "PNG(32): Bad header [%s].\n", name );
        goto exit;
    }

// inits
    png_ptr = png_create_read_struct(
                PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    if( !png_ptr ) {
        fprintf( flog,
        "PNG(32): Failed read_struct alloc [%s].\n", name );
        goto exit;
    }

    info_ptr = png_create_info_struct( png_ptr );

    if( !info_ptr ) {
        fprintf( flog,
        "PNG(32): Failed info_struct alloc [%s].\n", name );
        goto exit;
    }

    png_init_io( png_ptr, f );
    png_set_sig_bytes( png_ptr, 8 );

// read
    png_read_png( png_ptr, info_ptr,
        PNG_TRANSFORM_PACKING | PNG_TRANSFORM_SWAP_ENDIAN,
        NULL );

// format
    png_get_IHDR( png_ptr, info_ptr,
        &wi, &hi, &bit_depth, &color_type, NULL, NULL, NULL );

    w = wi;
    h = hi;

// get rows
    raster = (uint32*)malloc( w * h * sizeof(uint32) );

    if( !raster ) {
        fprintf( flog, "PNG(32): Failed raster alloc [%s].\n", name );
        goto exit;
    }

    if( color_type == PNG_COLOR_TYPE_GRAY ) {

        if( bit_depth == 16 ) {

            png_uint_16**	row_ptrs = (png_uint_16**)
                            png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
        else {	// 1, 2, 4, 8

            png_byte**	row_ptrs = (png_byte**)
                        png_get_rows( png_ptr, info_ptr );

            n = 0;

            for( y = 0; y < h; ++y ) {

                for( x = 0; x < w; ++x )
                    raster[n++] = row_ptrs[y][x];
            }
        }
    }
    else {

        // Note:
        // png_uint_32 is typedef'd as unsigned long
        // which is 64 bits on the cluster

        uint32**	row_ptrs = (uint32**)
                    png_get_rows( png_ptr, info_ptr );

        n = 0;

        for( y = 0; y < h; ++y ) {

            for( x = 0; x < w; ++x )
                raster[n++] = row_ptrs[y][x];
        }
    }

// report success
    ok = true;

// done
exit:
    if( png_ptr )
        png_destroy_read_struct( &png_ptr, &info_ptr, NULL );

    if( f )
        fclose( f );

    if( !ok ) {

        if( raster ) {
            free( raster );
            raster = NULL;
        }
    }

    return raster;
}

/* --------------------------------------------------------------- */
/* Raster8ToPng8 ------------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster8ToPng8(
    const char*		name,
    const uint8*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    FILE*			f = FileOpenOrDie( name, "w", flog );
    png_structp		png_ptr;
    png_infop		info_ptr;
    png_byte**		prow;

// init I/O
    png_ptr		= png_create_write_struct(
                    PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    info_ptr	= png_create_info_struct( png_ptr );

    png_init_io( png_ptr, f );

// header
    png_set_IHDR( png_ptr, info_ptr,
        w, h, 8, PNG_COLOR_TYPE_GRAY,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
        PNG_FILTER_TYPE_BASE );

    png_write_info( png_ptr, info_ptr );

// data
    prow = (png_byte**)malloc( h * sizeof(png_byte*) );

    for( int i = 0; i < h; ++i )
        prow[i] = (png_byte*)(raster + w * i);

    png_write_image( png_ptr, prow );
    free( prow );

// cleanup
    png_write_end( png_ptr, NULL );
    png_destroy_write_struct( &png_ptr, &info_ptr );
    fclose( f );
}

/* --------------------------------------------------------------- */
/* Raster16ToPng16 ----------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster16ToPng16(
    const char*		name,
    const uint16*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    FILE*			f = FileOpenOrDie( name, "w", flog );
    png_structp		png_ptr;
    png_infop		info_ptr;
    png_byte**		prow;

// init I/O
    png_ptr		= png_create_write_struct(
                    PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    info_ptr	= png_create_info_struct( png_ptr );

    png_init_io( png_ptr, f );

// header
    png_set_IHDR( png_ptr, info_ptr,
        w, h, 16, PNG_COLOR_TYPE_GRAY,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
        PNG_FILTER_TYPE_BASE );

    png_write_info( png_ptr, info_ptr );

// data
    prow = (png_byte**)malloc( h * sizeof(png_byte*) );

    for( int i = 0; i < h; ++i )
        prow[i] = (png_byte*)(raster + w * i);

    png_set_swap( png_ptr );
    png_write_image( png_ptr, prow );
    free( prow );

// cleanup
    png_write_end( png_ptr, NULL );
    png_destroy_write_struct( &png_ptr, &info_ptr );
    fclose( f );
}

/* --------------------------------------------------------------- */
/* Raster32ToPngRGBA --------------------------------------------- */
/* --------------------------------------------------------------- */

void Raster32ToPngRGBA(
    const char*		name,
    const uint32*	raster,
    int				w,
    int				h,
    FILE*			flog )
{
    FILE*			f = FileOpenOrDie( name, "w", flog );
    png_structp		png_ptr;
    png_infop		info_ptr;
    png_byte**		prow;

// init I/O
    png_ptr		= png_create_write_struct(
                    PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

    info_ptr	= png_create_info_struct( png_ptr );

    png_init_io( png_ptr, f );

// header
    png_set_IHDR( png_ptr, info_ptr,
        w, h, 8, PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
        PNG_FILTER_TYPE_BASE );

    png_write_info( png_ptr, info_ptr );

// data
    prow = (png_byte**)malloc( h * sizeof(png_byte*) );

    for( int i = 0; i < h; ++i )
        prow[i] = (png_byte*)(raster + w * i);

    png_write_image( png_ptr, prow );
    free( prow );

// cleanup
    png_write_end( png_ptr, NULL );
    png_destroy_write_struct( &png_ptr, &info_ptr );
    fclose( f );
}


