/* BK_BMP ------------------------------------------------------------
 *
 * Bitmap / Image operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_BMP.h"
#include	"BK_RGN.h"

#include	"BK_BMAP.h"
#include	"BK_HST.h"
#include	"BK_MEM.h"
#include	"BK_STAT.h"

#include	<stdio.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum {
/* maximum pixel value plus one */
    pixelValueRange	= 4096L,

/* don't bother looking for mode in top-most range */
    omitTopMost		= 128
};






/* BMPMapValue -------------------------------------------------------
 *
 * Return bit value at given location.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMPMapValue( const UInt32 *map, UInt32 hPix, int h, int v )
{
    return (map[v * BitToWord( hPix ) + BitToWord( h )]
            >> (AllButOneBit - WordRem( h ))) & 1;
}


/* BMPGetMode --------------------------------------------------------
 *
 * Return most probable (mode) image intensity.
 *
 * -----
 * Notes
 * -----
 *
 * This is one of the few places where we make use of the fact that
 * the data are really 12-bit instead of 16-bit. We assume we only
 * need pixelValueRange histogram bins. Even if flat-fielding were
 * to generate values slightly larger than 12 bits. That doesn't
 * bother us here for two reasons:
 *
 * (1) The histogram function has built-in overflow detection
 * to protect us from overrunning the bin array.
 *
 * (2) We don't care about very large values because the purpose
 * here is to find the background value which will always be
 * considerably lower than 4095.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMPGetMode(
    const UInt16	*data16Bit,
    UInt32			hPix,
    UInt32			vPix )
{
    UInt32	hist[pixelValueRange];

    HSTAltImageRowsUInt16( NULL, NULL, hist, pixelValueRange,
        data16Bit, hPix, vPix );

    return STATIndexOfMaxValueUInt32( hist,
            pixelValueRange - omitTopMost );
}


/* BMPGetModeBox -----------------------------------------------------
 *
 * Return most probable (mode) image intensity within box.
 *
 * -----
 * Notes
 * -----
 *
 * This is one of the few places where we make use of the fact that
 * the data are really 12-bit instead of 16-bit. We assume we only
 * need pixelValueRange histogram bins. Even if flat-fielding were
 * to generate values slightly larger than 12 bits. That doesn't
 * bother us here for two reasons:
 *
 * (1) The histogram function has built-in overflow detection
 * to protect us from overrunning the bin array.
 *
 * (2) We don't care about very large values because the purpose
 * here is to find the background value which will always be
 * considerably lower than 4095.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMPGetModeBox(
    const UInt16	*data16Bit,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box )
{
    UInt32	hist[pixelValueRange];

    HSTImageBoxUInt16( NULL, NULL, hist, pixelValueRange,
        data16Bit, hPix, vPix, box );

    return STATIndexOfMaxValueUInt32( hist,
            pixelValueRange - omitTopMost );
}


/* BMPGenerateThreshMap16Bit -----------------------------------------
 *
 * For 16-bit data...
 * Generate a one bit deep mask, setting all pixels
 * above threshold to one and others to zero.
 *
 * -----
 * Notes
 * -----
 *
 * - map will be a packed array of 32-bit words.
 * The storage needed is:
 * 32-bit words	= vPix X hPix / 32,
 * bytes		= vPix X hPix / 8.
 *
 * 16-bit pixel pairs are read into a 32-bit register and extracted
 * separately. On X86 family processors, the 16-bit number at the
 * lower RAM address occupies the lower 16 bits of the register.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMPGenerateThreshMap16Bit(
    const UInt16	*data16Bit,
    UInt32			*map,
    UInt32			hPix,
    UInt32			vPix,
    UInt32			thresh )
{
    UInt32	*pPair;
    UInt32	pair, W, word, nWords;
    int		bitPair;

    nWords	= BitToWord( vPix * hPix );
    pPair	= (UInt32*)data16Bit;
    W		= 0;	/* suppress initialization warning */

    for( word = 0; word < nWords; ++word ) {

        for( bitPair = 0; bitPair < HalfBits; ++bitPair ) {

            pair = *pPair++;

            W <<= 1;
            W += SignExtHiBit32( thresh - (pair & LoMask) ) & 1;

            W <<= 1;
            W += SignExtHiBit32( thresh - (pair >> HalfBits) ) & 1;
        }

        map[word] = W;
    }
}


/* BMPGenerateThreshPatch16Bit ---------------------------------------
 *
 * For 16-bit data...
 * Generate a one bit deep mask, setting all pixels in given patch
 * that are above threshold to one and others in patch to zero.
 *
 * -----
 * Notes
 * -----
 *
 * Pixels outside of patch are not touched. Caller may want to zero
 * the map prior to making this call.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMPGenerateThreshPatch16Bit(
    const UInt16	*data16Bit,
    UInt32			*map,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box,
    UInt32			thresh )
{
    UInt32	mapRowBytes, imgRowBytes, W = 0;
    int		v, h, bit, hL, hR, top, left, bot, right;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

    top		= box->top;
    left	= box->left;
    bot		= box->bottom;
    right	= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    hL	= BitToWord( left );
    hR	= BitToWord( right - 1 );

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top
                    + hL * WordBits * sizeof(UInt16));

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        const UInt16	*p16 = data16Bit;

        for( h = hL; h <= hR; ++h ) {

            for( bit = 0; bit < WordBits; ++bit ) {

                W <<= 1;
                W += SignExtHiBit32( thresh - *p16++ ) & 1;
            }

            map[h] = W;
        }
    }

exit:
    return;
}


/* BMPFillHolesBox ---------------------------------------------------
 *
 * Each blob is isolated and its holes are filled.
 *
 * The new blobs are accumulated in dstMap.
 *
 * -----
 * Notes
 * -----
 *
 * Tracking and scanning maps are consumed in the process.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMPFillHolesBox(
    UInt32			*dstMap,
    UInt32			*trackingMap,
    UInt32			*scanningMap,
    UInt32			*tmpMap0,
    UInt32			*tmpMap1,
    void			*bucketScratch,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box )
{
    UInt32	*pRow;
    U16Box	rBlob;
    UInt32	rowBytes;
    int		h, v, top, left, bot, right;

/* -------------------- */
/* Prepare new blob map */
/* -------------------- */

    MEMCopyMap( dstMap, scanningMap, hPix, vPix );

/* ---------------- */
/* Scan grid points */
/* ---------------- */

/* box limits */

    if( box ) {

        top		= box->top;
        left	= box->left;
        bot		= box->bottom;
        right	= box->right;

        if( bot > (int)vPix )
            bot = vPix;

        if( right > (int)hPix )
            right = hPix;

        if( top >= bot )
            goto exit;

        if( left >= right )
            goto exit;
    }
    else {
        top		= 0;
        left	= 0;
        bot		= vPix;
        right	= hPix;
    }

    rowBytes	= BitToByte( hPix );
    pRow		= (UInt32*)((char*)scanningMap + rowBytes * top);

    for( v = top; v < bot;
        ++v, pRow = (UInt32*)((char*)pRow + rowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( pRow, h ) ) {

                /* ------------- */
                /* Identify blob */
                /* ------------- */

                BMAPBucketTool( &rBlob, scanningMap,
                    (BMAPBucketPtr)bucketScratch,
                    hPix, vPix, h, v );

                /* ------------ */
                /* Isolate blob */
                /* ------------ */

                BMAPXorNewPatch( tmpMap0, scanningMap, trackingMap,
                    hPix, vPix, &rBlob );

                /* ---------- */
                /* Fill holes */
                /* ---------- */

                RGNGetHolesBlob( dstMap, tmpMap0, tmpMap1,
                    bucketScratch, hPix, vPix, &rBlob );

                /* ------------------- */
                /* Update tracking map */
                /* ------------------- */

                BMAPCopyPatch( trackingMap, scanningMap,
                    hPix, vPix, &rBlob );
            }
        }
    }

exit:
    return;
}


