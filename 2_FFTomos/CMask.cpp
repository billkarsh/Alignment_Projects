

#include	"CMask.h"

#include	"File.h"
#include	"ImageIO.h"

#include	"BK_BMAP.h"
#include	"BK_BMAP_CONVERT.h"
#include	"BK_BMP.h"
#include	"BK_MEM.h"
#include	"BK_RGN.h"
#include	"BK_STAT.h"

#include	<math.h>

#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum CStdDefConsts {
/* maximum pixel value plus one */
    pixelValueRange	= 65536L
};

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	SWAPMAPS( m1, m2 )										\
    do {UInt32 *temp = m1; m1 = m2; m2 = temp;} while( 0 )

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* ThreshFromPeak ------------------------------------------------ */
/* --------------------------------------------------------------- */

void CMask::ThreshFromPeak()
{
    if( threshMethod == stdDefMthdFixed )
        m_thresh = threshConst;
    else {

        m_thresh = m_peak;

        if( threshMethod == stdDefMthdAddFixed )
            m_thresh += threshConst;
        else if( threshMethod == stdDefMthdAddRel )
            m_thresh += threshConst * (FP32)sqrt( m_thresh );
        else
            m_thresh *= threshConst;
    }
}

/* --------------------------------------------------------------- */
/* GetMarkerPeak ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CMask::GetMarkerPeak()
{
    m_peak = (FP32)BMPGetModeBox( m_mrkPix,
                m_hPix, m_vPix, &m_roi );
}

/* --------------------------------------------------------------- */
/* GetThresh ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* GetThresh ---------------------------------------------------------
 *
 * Return marker image threshold.
 *
 * ========================================
 * Original code: Bill Karsh, 1/2009.
 *
 * Modifications
 * -------------
 *
 *
 * Copyright (c) 2009 HHMI.
 * All rights reserved.
 *
 */

UInt32 CMask::GetThresh()
{
    if( threshMethod != stdDefMthdFixed )
        GetMarkerPeak();

    ThreshFromPeak();

    return (UInt32)m_thresh;
}

/* --------------------------------------------------------------- */
/* LocalRethreshold ---------------------------------------------- */
/* --------------------------------------------------------------- */

/* LocalRethreshold --------------------------------------------------
 *
 * Each blob is isolated, its intensities are examined, then
 * the blob is rethresholded based on a selected parameter
 * of the intensity distribution.
 *
 * Note: Holes are also filled here if that option selected.
 *
 * The new blobs are accumulated into dstMap.
 *
 * -----
 * Notes
 * -----
 *
 * Tracking and scanning maps are consumed in the process.
 *
 * ========================================
 * Original code: Bill Karsh, 1/2009.
 *
 * Modifications
 * -------------
 *
 *
 * Copyright (c) 2009 HHMI.
 * All rights reserved.
 *
 */

void CMask::LocalRethreshold(
    UInt32	*dstMap,
    UInt32	*trackingMap,
    UInt32	*scanningMap,
    UInt32	*tmpMap0,
    UInt32	*tmpMap1,
    UInt32	T )
{
    UInt32			*hist;
    SInt32			*pRow;
    SUM_F_A_Rec		sum;
    U16Box			rBlob;
    UInt32			rowBytes, newT;
    int				hPix = m_hPix,
                    vPix = m_vPix,
                    h, v, top, left, bot, right;

/* ------------------ */
/* Clear new blob map */
/* ------------------ */

    MEMZeroMap( dstMap, m_hPix, m_vPix );

/* --------------------- */
/* Alloc histogram space */
/* --------------------- */

    hist = (UInt32*)malloc( pixelValueRange * sizeof(UInt32) );

    if( !hist )
        goto exit;

/* ---------------- */
/* Scan grid points */
/* ---------------- */

    top		= m_roi.top;
    left	= m_roi.left;
    bot		= m_roi.bottom;
    right	= m_roi.right;

    if( bot > vPix )
        bot = vPix;

    if( right > hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    rowBytes	= BitToByte( hPix );
    pRow		= (SInt32*)((char*)scanningMap + rowBytes * top);

    for( v = top; v < bot;
        ++v, pRow = (SInt32*)((char*)pRow + rowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( pRow, h ) ) {

                /* ------------- */
                /* Identify blob */
                /* ------------- */

                BMAPBucketTool( &rBlob, scanningMap,
                    (BMAPBucketPtr)m_scratch,
                    hPix, vPix, h, v );

                /* ------------ */
                /* Isolate blob */
                /* ------------ */

                BMAPXorNewPatch( tmpMap0, scanningMap, trackingMap,
                    hPix, vPix, &rBlob );

                /* ------------------ */
                /* Calc new threshold */
                /* ------------------ */

#if 0
// obsolete average method

                SUMBlob_F_A( &sum, m_mrkPix, tmpMap0,
                    hPix, vPix, &rBlob );

                newT = (UInt32)(sum.F * rethreshFrac / sum.A);

#else

                if( lrtMethod == stdLRTMthdPercentile ) {

                    sum.A = RGNHistogramBlob( NULL, NULL, hist,
                                m_mrkPix, tmpMap0,
                                hPix, vPix, &rBlob );

                    newT = STATPercentileBinUInt32(
                            hist, pixelValueRange,
                            sum.A, rethreshFrac );
                }
                else { // stdLRTMthdPercentPk

                    newT = (UInt32)(rethreshFrac
                            * RGNBlob_Ipk( m_mrkPix,
                                tmpMap0, hPix, vPix, &rBlob ));
                }
#endif

                /* ---------------- */
                /* Adjust this blob */
                /* ---------------- */

                RGNRethresholdBlob( m_mrkPix,
                    tmpMap0, hPix, vPix, newT, &rBlob );

                /* ----------------------------- */
                /* Add blob(s) into new blob map */
                /* ----------------------------- */

                BMAPOrPatch( dstMap, tmpMap0,
                    hPix, vPix, &rBlob );

                /* ------------------------ */
                /* Fill holes (if selected) */
                /* ------------------------ */

                if( fillHoles )
                    RGNGetHolesBlob( dstMap, tmpMap0, tmpMap1,
                        m_scratch, hPix, vPix, &rBlob );

                /* ------------------- */
                /* Update tracking map */
                /* ------------------- */

                BMAPCopyPatch( trackingMap, scanningMap,
                    hPix, vPix, &rBlob );
            }
        }
    }

exit:
    if( hist )
        free( hist );
}

/* --------------------------------------------------------------- */
/* InitObjectDetectionMaps --------------------------------------- */
/* --------------------------------------------------------------- */

/* InitObjectDetectionMaps -------------------------------------------
 *
 * Marker image is used to initialize three identical maps for
 * standard object scanning.
 *
 * Return primary global threshold value.
 *
 * ----------
 * Discussion
 * ----------
 *
 * In our standard scan procedure, two maps: {scanning, tracking}
 * are employed and both are progressively cleared as objects in
 * them are found and removed.
 *
 * The third map, besides serving as a workspace for this function,
 * often serves these purposes:
 *
 * - Recreation of scanning and tracking maps in follow-on
 * movie frames.
 *
 * - Deletion of marker object regions from signal channel regions
 * as protection against channel bleed-through, and to guard against
 * erroneous sampling of neighbor objects.
 *
 * -----
 * Notes
 * -----
 *
 * The expensive part of hole filling is not the actual hole fill
 * operation on each blob (RGNGetHolesBlob). Rather, it is the
 * finding of the blobs using the bucket tool. Noting that local
 * rethresholding already uses a blob finding loop, we can place
 * the hole filler right into that loop. On the other hand, if
 * rethresholding is not selected, we use a separate blob finding
 * loop just for hole filling.
 *
 * ========================================
 * Original code: Bill Karsh, 1/2009.
 *
 * Modifications
 * -------------
 *
 *
 * Copyright (c) 2009 HHMI.
 * All rights reserved.
 *
 */

UInt32 CMask::InitObjectDetectionMaps(
    UInt32	*auxObjectMap,
    UInt32	*trackingMap,
    UInt32	*scanningMap,
    UInt32	*tmpMap0,
    UInt32	*tmpMap1 )
{
    UInt32	*mA, *mB, *mC;
    UInt32	T;

/* ------------------------- */
/* Initial alias assignments */
/* ------------------------- */

    mA = auxObjectMap;	// mA always pointed to current map
    mB = trackingMap;
    mC = scanningMap;

/* ------------------------------- */
/* First pass global threshold map */
/* ------------------------------- */

    T = GetThresh();

    BMPGenerateThreshMap16Bit( m_mrkPix, mA,
        m_hPix, m_vPix, T );

/* --------------------------------------------- */
/* Local rethresholding (and optional hole fill) */
/* --------------------------------------------- */

    if( lrtMethod > stdLRTMthdNone ) {

        MEMCopyMap( mB, mA, m_hPix, m_vPix );

        LocalRethreshold( mC, mA, mB, tmpMap0, tmpMap1, T );

        SWAPMAPS( mA, mC );
    }

/* ---------------------------------- */
/* Fill holes (if not rethresholding) */
/* ---------------------------------- */

    else if( fillHoles ) {

        MEMCopyMap( mB, mA, m_hPix, m_vPix );

        BMPFillHolesBox( mC, mA, mB, tmpMap0, tmpMap1,
            m_scratch, m_hPix, m_vPix, &m_roi );

        SWAPMAPS( mA, mC );
    }

/* ----- */
/* Erode */
/* ----- */

    if( erodePix ) {

        int		mapIdx;

        /* ---------------------- */
        /* Erode into unused maps */
        /* ---------------------- */

        mapIdx = BMAPInsetMap_NPixels( mB, mC, mA,
                    m_hPix, m_vPix, erodePix );

        /* -------------- */
        /* Update aliases */
        /* -------------- */

        if( mapIdx )
            SWAPMAPS( mA, mC );
        else
            SWAPMAPS( mA, mB );
    }

/* -------------------- */
/* Replicate result map */
/* -------------------- */

    MEMCopyMap( mB, mA, m_hPix, m_vPix );
    MEMCopyMap( mC, mA, m_hPix, m_vPix );

/* ---- */
/* Done */
/* ---- */

    return T;
}

/* --------------------------------------------------------------- */
/* BBoxPass ------------------------------------------------------ */
/* --------------------------------------------------------------- */

/* BBoxPass ----------------------------------------------------------
 *
 * Return true if object passes the standard bounding-box filter.
 *
 * ========================================
 * Original code: Bill Karsh, 1/2009.
 *
 * Modifications
 * -------------
 *
 *
 * Copyright (c) 2009 HHMI.
 * All rights reserved.
 *
 */

int CMask::BBoxPass( const U16BoxPtr rBlob )
{
    int		pass = true;

    if( checkBBox ) {

        int		min = rBlob->bottom - rBlob->top,
                max = rBlob->right  - rBlob->left;

        if( min > max ) {

            int		t = min;

            min	= max;
            max	= t;
        }

        if( min < minBBox )
            pass = false;
        else if( max > maxBBox )
            pass = false;
    }

    return pass;
}

/* --------------------------------------------------------------- */
/* ApplyFilters -------------------------------------------------- */
/* --------------------------------------------------------------- */

/* Marker_I_Pass -----------------------------------------------------
 *
 * Return true if blob passes standard marker intensity filter.
 *
 * ========================================
 * Original code: Bill Karsh, 1/2009.
 *
 * Modifications
 * -------------
 *
 *
 * Copyright (c) 2009 HHMI.
 * All rights reserved.
 *
 */

int CMask::Marker_I_Pass( const UInt32 *map, const U16BoxPtr rBlob )
{
    SUM_F_A_Rec		sum;
    int				I, pass = true;

    if( check_I ) {

        SUMBlob_F_A( &sum, m_mrkPix, map,
            m_hPix, m_vPix, rBlob );

        I = (sum.A ? sum.F / sum.A : 0);

        if( I < min_I || I > max_I )
            pass = false;
    }

    return pass;
}

/* --------------------------------------------------------------- */
/* ApplyFilters -------------------------------------------------- */
/* --------------------------------------------------------------- */

/* ApplyFilters ------------------------------------------------------
 *
 * Find each object in the scnMap/trkMap pair...
 *
 * If it passes standard flt criteria, OR it into the fltMap.
 *
 */

void CMask::ApplyFilters()
{
    UInt32	*pRow32;
    U16Box	R;
    UInt32	dRow32;
    int		h, v, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    MEMZeroMap( m_fltMap, m_hPix, m_vPix );

/* ---------------- */
/* Scan grid points */
/* ---------------- */

    top		= m_roi.top;
    left	= m_roi.left;
    bot		= m_roi.bottom;
    right	= m_roi.right;

    dRow32	= BitToByte( m_hPix );
    pRow32	= (UInt32*)((char*)m_scnMap + top * dRow32);

    for( v = top; v < bot;
        ++v, pRow32 = (UInt32*)((char*)pRow32 + dRow32) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( pRow32, h ) ) {

                /* ------------- */
                /* Identify blob */
                /* ------------- */

                BMAPBucketTool( &R, m_scnMap,
                (BMAPBucketPtr)m_scratch, m_hPix, m_vPix, h, v );

                /* ----------- */
                /* Size filter */
                /* ----------- */

                if( !BBoxPass( &R ) )
                    goto updateTracking;

                /* ---------------- */
                /* Intensity filter */
                /* ---------------- */

                BMAPXorNewPatch( m_isoMap, m_scnMap, m_trkMap,
                    m_hPix, m_vPix, &R );

                if( Marker_I_Pass( m_isoMap, &R ) ) {

                    // add to fltMap
                    BMAPOrPatch( m_fltMap, m_isoMap,
                        m_hPix, m_vPix, &R );
                }

                BMAPZeroPatch( m_isoMap, m_hPix, m_vPix, &R );

                /* ------------------- */
                /* Update tracking map */
                /* ------------------- */

updateTracking:
                BMAPCopyPatch( m_trkMap, m_scnMap,
                    m_hPix, m_vPix, &R );
            }
        }
    }
}

/* --------------------------------------------------------------- */
/* ReadParamFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CMask::ReadParamFile( FILE *flog, const char *path )
{
    this->flog = flog;

    FILE*		f = FileOpenOrDie( path, "r" );
    CLineScan	LS;

    while( LS.Get( f ) > 0 ) {

        int		i;
        char	c;

        if( LS.line[1] != 'O' )
            continue;

        if( !strncmp( LS.line, "<Object_Definition>", 19 ) ) {

            LS.Get( f );
            LS.Get( f );
            sscanf( LS.line,
            "<Threshold method=\"%d\" const=\"%f\"",
            &i, &threshConst );
            threshMethod = i;

            LS.Get( f );
            LS.Get( f );
            sscanf( LS.line,
            "<Rethreshold method=\"%d\" frac=\"%f\"",
            &i, &rethreshFrac );
            lrtMethod = i;

            LS.Get( f );
            sscanf( LS.line,
            "<Refinements fill_holes=\"%c\" erode_pix=\"%d\"",
            &c, &erodePix );
            fillHoles = (toupper( c ) == 'T');
        }
        else if( !strncmp( LS.line, "<Object_Filters>", 16 ) ) {

            LS.Get( f );
            sscanf( LS.line,
            "<BoundingBox enabled=\"%c\" min=\"%d\" max=\"%d\"",
            &c, &minBBox, &maxBBox );
            checkBBox = (toupper( c ) == 'T');

            LS.Get( f );
            sscanf( LS.line,
            "<Intensity enabled=\"%c\" min=\"%d\" max=\"%d\"",
            &c, &min_I, &max_I );
            check_I = (toupper( c ) == 'T');

            break;
        }
    }

    fclose( f );
}

/* --------------------------------------------------------------- */
/* MaskFromImage ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CMask::MaskFromImage(
    const char		*outpath,
    const UInt16	*src,
    int				hPix,
    int				vPix )
{
    m_mrkPix		= src;
    m_hPix			= hPix;
    m_vPix			= vPix;
    m_roi.top		= 0;
    m_roi.left		= 0;
    m_roi.bottom	= vPix;
    m_roi.right		= hPix;

    InitObjectDetectionMaps( m_iniMap, m_trkMap, m_scnMap,
        m_tmpMap0, m_tmpMap1 );

    ApplyFilters();

    vector<UInt8>	dst( hPix * vPix );

    BMAPConvertDepth1To8( &dst[0], m_fltMap, hPix, vPix, 1 );

    Raster8ToTif8( outpath, &dst[0], hPix, vPix, flog );
}


