/* BK_RGN ------------------------------------------------------------
 *
 * Operations on masked image regions.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_RGN.h"

#include	"BK_BMAP.h"
#include	"BK_BMAP_SET.h"
#include	"BK_GEOM.h"
#include	"BK_MEM.h"
#include	"BK_STAT.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum {
/* maximum pixel value plus one */
    pixelValueRange	= 4096L
};






/* RGNBox_Ipk --------------------------------------------------------
 *
 * Within the intersection of the given blob and box...
 *
 * Get value and coordinates of the peak-intensity pixel.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNBox_Ipk(
    RGN_Ipk_Ptr			ipk,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes;
    int		v, h, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    ipk->Ipk	= 0;
    ipk->v		= top;
    ipk->h		= left;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* ---------------------- */
/* Loop over boxed pixels */
/* ---------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) ) {

                UInt32	i;

                if( (i = data16Bit[h]) > ipk->Ipk ) {

                    ipk->Ipk	= i;
                    ipk->v		= v;
                    ipk->h		= h;
                }
            }
        }
    }

exit:
    return;
}


/* RGNBox_F_A --------------------------------------------------------
 *
 * Within the intersection of the given blob and box...
 *
 * Get total flux and area.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNBox_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes, F, A;
    int		v, h, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    F = 0;
    A = 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* ---------------------- */
/* Loop over boxed pixels */
/* ---------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) ) {

                F += data16Bit[h];
                ++A;
            }
        }
    }

exit:
    sum->F = F;
    sum->A = A;
}


/* RGNBox_F_F2_A -----------------------------------------------------
 *
 * Within the intersection of the given blob and box...
 *
 * Get total flux, sum(flux squared) and area.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNBox_F_F2_A(
    SUM_F_F2_A_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes, F, F2, A;
    int		v, h, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    F	= 0;
    F2	= 0;
    A	= 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* ---------------------- */
/* Loop over boxed pixels */
/* ---------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) ) {

                UInt32	f = data16Bit[h];

                F	+= f;
                F2	+= f * f;
                ++A;
            }
        }
    }

exit:
    sum->F	= F;
    sum->F2	= F2;
    sum->A	= A;
}


/* RGNBoxPixList -----------------------------------------------------
 *
 * Copy all pixels within box to list array.
 *
 * Return pixel count.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNBoxPixList(
    UInt32				*list,
    const UInt16		*data16Bit,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	imgRowBytes;
    int		N, v, h, top, left, bot, right;

    N			= 0;
    imgRowBytes	= hPix * sizeof(UInt16);

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

    for( v = top; v < bot; ++v,
        data16Bit = (UInt16*)((char*)data16Bit + imgRowBytes) ) {

        for( h = left; h < right; ++h, ++N )
            list[N] = data16Bit[h];
    }

    return N;
}


/* RGNHistogramBox ---------------------------------------------------
 *
 * Within the intersection of the given blob and box...
 *
 * Histogram all pixels.
 *
 * Return pixel count, excluding overflow.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNHistogramBox(
    UInt32				*oFlowCnt,
    UInt32				*oFlowSum,
    UInt32				*binArray,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes, N, x, oflowC, oflowS;
    int		v, h, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    MEMZeroBytes( binArray, pixelValueRange * sizeof(UInt32) );
    N			= 0;
    oflowC		= 0;
    oflowS		= 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* ---------------------- */
/* Loop over boxed pixels */
/* ---------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) ) {

                if( (x = data16Bit[h]) < pixelValueRange ) {
                    ++binArray[x];
                    ++N;
                }
                else {
                    oflowS += x;
                    ++oflowC;
                }
            }
        }
    }

exit:
    if( oFlowCnt )
        *oFlowCnt = oflowC;

    if( oFlowSum )
        *oFlowSum = oflowS;

    return N;
}


/* RGNBoxPerim_F_A ---------------------------------------------------
 *
 * Along perimeter of intersection of given blob and box...
 *
 * Get total flux and area.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNBoxPerim_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes, F, A;
    int		v, h, L, R, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    F = 0;
    A = 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    L		= left;
    R		= right - 1;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* --- */
/* Top */
/* --- */

    for( h = left; h < right; ++h ) {

        if( BMAP_IsRowBit( map, h ) ) {

            F += data16Bit[h];
            ++A;
        }
    }

    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes);
    map			= (UInt32*)((char*)map + mapRowBytes);

/* ----- */
/* Sides */
/* ----- */

    --right;

    for( v = top + 2; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        if( left >= L && BMAP_IsRowBit( map, left ) ) {

            F += data16Bit[left];
            ++A;
        }

        if( right <= R && BMAP_IsRowBit( map, right ) ) {

            F += data16Bit[right];
            ++A;
        }
    }

    ++right;

/* ------ */
/* Bottom */
/* ------ */

    if( bot - top > 1 ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) ) {

                F += data16Bit[h];
                ++A;
            }
        }
    }

exit:
    sum->F = F;
    sum->A = A;
}


/* RGNBoxPerimList ---------------------------------------------------
 *
 * Along perimeter of intersection of given blob and box...
 *
 * Copy intensities into array.
 *
 * Return pixel count.
 *
 * -----
 * Notes
 * -----
 *
 * Caller is responsible for providing list with enough space
 * to hold all perimeter pixels.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNBoxPerimList(
    UInt32				*list,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box )
{
    UInt32	mapRowBytes, imgRowBytes, N;
    int		v, h, L, R, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    N			= 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

/* blob limits */

    top		= rBlob->top;
    left	= rBlob->left;
    bot		= rBlob->bottom;
    right	= rBlob->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    L		= left;
    R		= right - 1;

/* box limits */

    if( box->top > top )
        top = box->top;

    if( box->left > left )
        left = box->left;

    if( box->bottom < bot )
        bot = box->bottom;

    if( box->right < right )
        right = box->right;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    map			= (UInt32*)((char*)map + mapRowBytes * top);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* --- */
/* Top */
/* --- */

    for( h = left; h < right; ++h ) {

        if( BMAP_IsRowBit( map, h ) )
            list[N++] = data16Bit[h];
    }

    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes);
    map			= (UInt32*)((char*)map + mapRowBytes);

/* ----- */
/* Sides */
/* ----- */

    --right;

    for( v = top + 2; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        if( left >= L && BMAP_IsRowBit( map, left ) )
            list[N++] = data16Bit[left];

        if( right <= R && BMAP_IsRowBit( map, right ) )
            list[N++] = data16Bit[right];
    }

    ++right;

/* ------ */
/* Bottom */
/* ------ */

    if( bot - top > 1 ) {

        for( h = left; h < right; ++h ) {

            if( BMAP_IsRowBit( map, h ) )
                list[N++] = data16Bit[h];
        }
    }

exit:
    return N;
}


/* RGNBlob_A ---------------------------------------------------------
 *
 * Return blob area.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNBlob_A(
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    return BMAPPatchArea( map, hPix, vPix, rBlob );
}


/* RGNBlob_F ---------------------------------------------------------
 *
 * Return total blob flux.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNBlob_F(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    SUM_F_A_Rec		sum;

    SUMBlob_F_A( &sum, data16Bit, map, hPix, vPix, rBlob );

    return sum.F;
}


/* RGNBlob_I ---------------------------------------------------------
 *
 * Return blob intensity.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 RGNBlob_I(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    SUM_F_A_Rec		sum;

    SUMBlob_F_A( &sum, data16Bit, map, hPix, vPix, rBlob );

    return SUM_RecIntensity( &sum );
}


/* RGNBlob_Ipk -------------------------------------------------------
 *
 * Return peak blob intensity.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNBlob_Ipk(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, W, val, Ipk;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    Ipk			= 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

    top			= rBlob->top;
    left		= rBlob->left;
    bot			= rBlob->bottom;
    right		= rBlob->right;

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

/* ------------------------------- */
/* Loop over 32-bit words in patch */
/* ------------------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        pPair = (UInt32*)data16Bit;

        for( h = hL; h <= hR; ++h, pPair += HalfBits ) {

            if( W = map[h] ) {

                /* --------------------------- */
                /* Loop over bit pairs in word */
                /* --------------------------- */

                for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                    pair = pPair[pairIdx];

                    if( IsHiBitSet32( W ) ) {

                        val = pair & LoMask;

                        if( val > Ipk )
                            Ipk = val;
                    }

                    W <<= 1;

                    if( IsHiBitSet32( W ) ) {

                        val = pair >> HalfBits;

                        if( val > Ipk )
                            Ipk = val;
                    }

                    W <<= 1;
                }
            }
        }
    }

exit:
    return Ipk;
}


/* RGNHistogramBlob --------------------------------------------------
 *
 * All the pixels composing the current blob are histogrammed.
 *
 * Return pixel count, excluding overflow.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 RGNHistogramBlob(
    UInt32				*oFlowCnt,
    UInt32				*oFlowSum,
    UInt32				*binArray,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, W, N, x, oflowC, oflowS;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    MEMZeroBytes( binArray, pixelValueRange * sizeof(UInt32) );
    N			= 0;
    oflowC		= 0;
    oflowS		= 0;

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

    top			= rBlob->top;
    left		= rBlob->left;
    bot			= rBlob->bottom;
    right		= rBlob->right;

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

/* ------------------------------- */
/* Loop over 32-bit words in patch */
/* ------------------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        pPair = (UInt32*)data16Bit;

        for( h = hL; h <= hR; ++h, pPair += HalfBits ) {

            if( W = map[h] ) {

                /* --------------------------- */
                /* Loop over bit pairs in word */
                /* --------------------------- */

                for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                    pair = pPair[pairIdx];

                    if( IsHiBitSet32( W ) ) {

                        if( (x = pair & LoMask) <
                            pixelValueRange ) {

                            ++binArray[x];
                            ++N;
                        }
                        else {
                            oflowS += x;
                            ++oflowC;
                        }
                    }

                    W <<= 1;

                    if( IsHiBitSet32( W ) ) {

                        if( (x = pair >> HalfBits) <
                            pixelValueRange ) {

                            ++binArray[x];
                            ++N;
                        }
                        else {
                            oflowS += x;
                            ++oflowC;
                        }
                    }

                    W <<= 1;
                }
            }
        }
    }

exit:
    if( oFlowCnt )
        *oFlowCnt = oflowC;

    if( oFlowSum )
        *oFlowSum = oflowS;

    return N;
}


/* RGNBlobGrainsV0 ---------------------------------------------------
 *
 * Return count of granules in blob (version 0).
 *
 * sums
 * ----
 * If non-NULL, gets total flux and area over all grains
 * and (if granule count non-zero) for original blob.
 *
 * visMap
 * ------
 * If non-NULL, the bounding boxes of the blob and of the
 * qualified granules are drawn into visMap.
 *
 * -----
 * Notes
 * -----
 *
 * - granPix must be >= 2.
 *
 * - corPix must be <= granPix;
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 RGNBlobGrainsV0(
    RGN_Grain_Ptr		sums,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    int					granPix,
    int					corePix,
    FP32				minGradient,
    UInt32				*visMap )
{
    U16Box	rGrid;
    int		N;

/* ---- */
/* Init */
/* ---- */

    N = 0;

    if( sums )
        MEMZeroBytes( sums, sizeof(RGN_Grain_Rec) );

    if( visMap )
        BMAPSetBox( visMap, hPix, vPix, rBlob );

/* -------------------- */
/* For each grid box... */
/* -------------------- */

    rGrid.top		= rBlob->top;
    rGrid.bottom	= rGrid.top + granPix;

    for( ; rGrid.bottom <= rBlob->bottom;
        rGrid.top = rGrid.bottom, rGrid.bottom += granPix ) {

        rGrid.left	= rBlob->left;
        rGrid.right	= rGrid.left + granPix;

        for( ; rGrid.right <= rBlob->right;
            rGrid.left = rGrid.right, rGrid.right += granPix ) {

            RGN_Ipk_Rec		I1, I2;
            SUM_F_A_Rec		fPeri, fCore;
            U16Box			rGran, rCore;

            /* --------------------- */
            /* Locate bightest pixel */
            /* --------------------- */

            RGNBox_Ipk( &I1, data16Bit, map, hPix, vPix,
                rBlob, &rGrid );

            /* ---------------------- */
            /* Center rGran box on it */
            /* ---------------------- */

            rGran.top		= I1.v - granPix / 2;
            rGran.left		= I1.h - granPix / 2;
            rGran.bottom	= rGran.top  + granPix;
            rGran.right		= rGran.left + granPix;

            /* ----------------------------------------- */
            /* Verify center pixel is brightest in rGran */
            /* ----------------------------------------- */

            RGNBox_Ipk( &I2, data16Bit, map, hPix, vPix,
                rBlob, &rGran );

            if( I2.v != I1.v || I2.h != I1.h )
                continue;

            /* ----------------------- */
            /* Measure rGran perimeter */
            /* ----------------------- */

            RGNBoxPerim_F_A( &fPeri, data16Bit, map, hPix, vPix,
                rBlob, &rGran );

            if( !fPeri.F )
                continue;

            /* -------------------- */
            /* Measure granule core */
            /* -------------------- */

            rCore.top		= I1.v - corePix / 2;
            rCore.left		= I1.h - corePix / 2;
            rCore.bottom	= rCore.top  + corePix;
            rCore.right		= rCore.left + corePix;

            RGNBox_F_A( &fCore, data16Bit, map, hPix, vPix,
                rBlob, &rCore );

            if( !fCore.A )
                continue;

            /* ------- */
            /* Qualify */
            /* ------- */

            if( (FP32)(fCore.F * fPeri.A) / (fCore.A * fPeri.F)
                >= minGradient ) {

                ++N;

                if( sums ) {

                    SUM_F_A_Rec	s;

                    RGNBox_F_A( &s, data16Bit, map, hPix, vPix,
                        rBlob, &rGran );

                    sums->grains.F += s.F;
                    sums->grains.A += s.A;
                }

                if( visMap )
                    BMAPSetBox( visMap, hPix, vPix, &rGran );
            }
        }
    }

    if( sums && N )
        SUMBlob_F_A( &sums->blob, data16Bit, map, hPix, vPix, rBlob );

    return (FP32)N;
}


/* RGNRethresholdBlob ------------------------------------------------
 *
 * Given a patch containing an isolated blob...
 *
 * Clear any blob pixel in the map whose intensity is below thresh.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNRethresholdBlob(
    const UInt16		*data16Bit,
    UInt32				*map,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				thresh,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, W, bit, keep;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    mapRowBytes	= BitToByte( hPix );
    imgRowBytes	= hPix * sizeof(UInt16);

    top			= rBlob->top;
    left		= rBlob->left;
    bot			= rBlob->bottom;
    right		= rBlob->right;

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

/* ------------------------------- */
/* Loop over 32-bit words in patch */
/* ------------------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        map			= (UInt32*)((char*)map + mapRowBytes) ) {

        pPair = (UInt32*)data16Bit;

        for( h = hL; h <= hR; ++h, pPair += HalfBits ) {

            if( W = map[h] ) {

                bit		= HiBit;
                keep	= 0;

                /* --------------------------- */
                /* Loop over bit pairs in word */
                /* --------------------------- */

                for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                    pair = pPair[pairIdx];

                    if( IsHiBitSet32( W ) ) {

                        if( (pair & LoMask) >= thresh )
                            keep |= bit;
                    }

                    W	<<= 1;
                    bit	>>= 1;

                    if( IsHiBitSet32( W ) ) {

                        if( (pair >> HalfBits) >= thresh )
                            keep |= bit;
                    }

                    W	<<= 1;
                    bit	>>= 1;
                }

                /* ---------------------------- */
                /* Update map word with keepers */
                /* ---------------------------- */

                map[h] = keep;
            }
        }
    }

exit:
    return;
}


/* RGNErodeBlob ------------------------------------------------------
 *
 * Erode blob N pixels using two caller-provided workspaces.
 *
 * Return pointer to result workspace.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * The storage spaces: srcMap, erMap0, erMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 *RGNErodeBlob(
    UInt32				*erMap0,
    UInt32				*erMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    UInt32				N )
{
    int		i;

/* -------------------------- */
/* Prepare erosion workspaces */
/* -------------------------- */

    BMAPZeroPatch( erMap0, hPix, vPix, rBlob );
    BMAPZeroPatch( erMap1, hPix, vPix, rBlob );

/* ----- */
/* Erode */
/* ----- */

    i = BMAPInsetPatch_NPixels( erMap0, erMap1, srcMap,
            hPix, vPix, N, rBlob );

    return (i ? erMap1 : erMap0);
}


/* RGNDilateBlob -----------------------------------------------------
 *
 * Dilate blob N pixels using two caller-provided workspaces.
 *
 * Return pointer to result workspace.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * (1) The effect of this function is to grow the source object
 * outward in all directions. The bounding box surrounding the
 * dilated object will be larger than that surrounding the source
 * object. Therefore, the dilated object will live in a larger
 * patch than would be needed for the undilated source. This
 * function will need the larger box/patch as its input...
 *
 * BEFORE calling this function, the caller must prepare the
 * source object as follows:
 *
 * (1a) If the object's undilated bounding box is denoted by
 * origBox, then the caller must calculate the expanded rBlob
 * (to be used by this function) as follows:
 *
 *		GEOMOutsetU16Box( &rBlob, &origBox, hPix, vPix, N ).
 *
 * (1b) The patch containing the source object must correspond to
 * the expanded rBlob parameter. That is, the source object should
 * be isolated into a patch (or map) that has been pre-zeroed out
 * to the larger size.
 *
 * (2) The storage spaces: srcMap, dilMap0, dilMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 *RGNDilateBlob(
    UInt32				*dilMap0,
    UInt32				*dilMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    UInt32				N )
{
    int		i;

    i = BMAPTracePatch_NPixels( dilMap0, dilMap1, srcMap,
            hPix, vPix, N, rBlob );

    return (i ? dilMap1 : dilMap0);
}


/* RGNGetHolesBlob ---------------------------------------------------
 *
 * Calculate holes in source blob and OR them into dstMap.
 *
 * A pixel is in a hole if it has no path to the blob's
 * bounding box.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNGetHolesBlob(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				*tmpMap,
    void				*bucketScratch,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    U16Box	R;
    int		h, v;

/* ------------------------------------- */
/* One-pixel-larger box: blob's exterior */
/* ------------------------------------- */

    GEOMOutsetU16Box( &R, rBlob, hPix, vPix, 1 );

/* -------------------- */
/* Perimeter seed point */
/* -------------------- */

    h = R.left;
    v = R.top;

    if( !rBlob->left ) {

        if( (h = rBlob->right) >= (int)hPix )
            goto exit;
    }

    if( !rBlob->top ) {

        if( (v = rBlob->bottom) >= (int)vPix )
            goto exit;
    }

/* ------------------------------ */
/* Create inverted blob in tmpMap */
/* ------------------------------ */

    BMAPFillPatch( tmpMap, hPix, vPix, &R );

    BMAPBicPatch( tmpMap, srcMap, hPix, vPix, rBlob );

/* ------------------------------- */
/* Remove connections to perimeter */
/* ------------------------------- */

    BMAPBoxedBucketTool( &R, tmpMap,
        (BMAPBucketPtr)bucketScratch,
        hPix, vPix, h, v, &R );

/* --------------------------------- */
/* OR remaining holes to destination */
/* --------------------------------- */

    BMAPOrBox( dstMap, tmpMap, hPix, vPix, rBlob );

exit:
    return;
}


/* RGNDelTouchingRegion ----------------------------------------------
 *
 * Delete from the patch specified by {dstMap, rDst} all objects
 * that touch the killing region specified by {rgnMap, rRgn}.
 *
 * -----
 * Notes
 * -----
 *
 * This is a patch-type operation.
 *
 * ------
 * Method
 * ------
 *
 * The basic idea here is straightforward: The bucket tools can be
 * used to delete logical objects. By ORing a killing object into
 * the destination and then using a bucket tool on a seed point in
 * the killing object we can delete the killer, and by extension,
 * any destination objects that touch it.
 *
 * There are several interesting subtleties in the implementation.
 *
 * To allow the caller to use a destination patch that is not in
 * an entirely clean map we have to use a constrained form of the
 * bucket tool: BMAPBoxedBucketTool. What box should we use as the
 * constraint? The following will explain...
 *
 * One novelty is that this function can be used on two patches that
 * are non-coincident. We only really care about that part of the
 * killer that overlaps the destination. That leads us to work with
 * the intersection between the two regions.
 *
 * A patch is always specified by a box that it covers. Most of the
 * time, such boxes are bounding boxes. If that were always the
 * case for the inputs to this function, then we could use those
 * boxes for calculating the intersections, and to constrain the
 * boxed bucket tool.
 *
 * In general, to allow that either box may not be a bounding box,
 * we work with full patch boxes.
 *
 * The killing region may be discontiguous. We need a loop because
 * several seed points may be required to exhaust the killer pixels.
 * BMAPGetSeedPatch is used to do an efficient search for the next
 * seed point. Each time we delete something from the intersection
 * map D, we BIC the deleted killing object from map K.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void RGNDelTouchingRegion(
    UInt32				*dstMap,
    const U16BoxPtr		rDst,
    const UInt32		*rgnMap,
    const U16BoxPtr		rRgn,
    UInt32				*tmpMap0,
    UInt32				*tmpMap1,
    void				*bucketScratch,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	*K = tmpMap0,
            *D = tmpMap1;
    U16Box	rK, rDstPatch, rRgnPatch;
    UInt32	hSeed, vSeed;

/* ---------------------------------------------------- */
/* Work with entire patch content, not just BBox subset */
/* ---------------------------------------------------- */

    GEOMU16PatchBox( &rDstPatch, rDst, hPix, vPix );
    GEOMU16PatchBox( &rRgnPatch, rRgn, hPix, vPix );

/* --------------------------------------------------- */
/* Limit killing region to intersection of two patches */
/* --------------------------------------------------- */

    if( !GEOMU16BoxIntersection( &rK, &rDstPatch, &rRgnPatch ) )
        goto exit;

/* --------------------- */
/* Killing region empty? */
/* --------------------- */

    if( !BMAPGetSeedPatch( &hSeed, &vSeed, rgnMap,
            hPix, vPix, &rK, rK.top ) ) {

        goto exit;
    }

/* ------------------------------------------ */
/* Store killing region for tracking purposes */
/* ------------------------------------------ */

    BMAPCopyPatch( K, rgnMap, hPix, vPix, &rK );

/* ---------------------------------------- */
/* OR region into dstMap to locate touchers */
/* ---------------------------------------- */

    BMAPOrPatch( dstMap, K, hPix, vPix, &rK );

/* ---------------------------------- */
/* Delete touchers until region empty */
/* ---------------------------------- */

    do {

        U16Box	rDum;

        /* ----------------------------------------- */
        /* Save copy of intersection before deletion */
        /* ----------------------------------------- */

        BMAPCopyPatch( D, dstMap, hPix, vPix, &rK );

        /* ------------------------ */
        /* Delete using object find */
        /* ------------------------ */

        BMAPBoxedBucketTool( &rDum, dstMap,
            (BMAPBucketPtr)bucketScratch,
            hPix, vPix, hSeed, vSeed, &rDstPatch );

        /* ---------------------- */
        /* Update tracking region */
        /* ---------------------- */

        BMAPBicDifPatch( K, D, dstMap, hPix, vPix, &rK );

        /* --------------------- */
        /* Killing region empty? */
        /* --------------------- */

    } while( BMAPGetSeedPatch( &hSeed, &vSeed, K,
                hPix, vPix, &rK, vSeed ) );

exit:
    return;
}


