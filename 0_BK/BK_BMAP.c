/* BK_BMAP -----------------------------------------------------------
 *
 * Bitmap operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_BMAP.h"
#include	"BK_MEM.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	Mask32_ZerosLeftOf( left )								\
    (((UInt32)-1) >> WordRem( left ))

#define	Mask32_ZerosRightOf( right )							\
    (((UInt32)-1) << (AllButOneBit - WordRem( right )))






/* BMAPOrMaps --------------------------------------------------------
 *
 * Generate the logical OR of two maps.
 *
 * -----
 * Notes
 * -----
 *
 * orMap is permitted to occupy the same storage space
 * as either mapA or mapB.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPOrMaps(
    UInt32				*orMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	word, nWords;

    nWords = BitToWord( vPix * hPix );

    for( word = 0; word < nWords; ++word )
        *orMap++ = *mapA++ | *mapB++;
}


/* BMAPXorMaps -------------------------------------------------------
 *
 * Generate the logical XOR of two maps.
 *
 * -----
 * Notes
 * -----
 *
 * xorMap is permitted to occupy the same storage space
 * as either mapA or mapB.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPXorMaps(
    UInt32				*xorMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	word, nWords;

    nWords = BitToWord( vPix * hPix );

    for( word = 0; word < nWords; ++word )
        *xorMap++ = *mapA++ ^ *mapB++;
}


/* BMAPAndMaps -------------------------------------------------------
 *
 * Generate the logical AND of two maps.
 *
 * -----
 * Notes
 * -----
 *
 * andMap is permitted to occupy the same storage space
 * as either mapA or mapB.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPAndMaps(
    UInt32				*andMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	word, nWords;

    nWords = BitToWord( vPix * hPix );

    for( word = 0; word < nWords; ++word )
        *andMap++ = *mapA++ & *mapB++;
}


/* BMAPBicMaps -------------------------------------------------------
 *
 * Generate mapA AND NOT(mapB), also called known as a bit clear
 * (BIC) operation.
 *
 * -----
 * Notes
 * -----
 *
 * bicMap is permitted to occupy the same storage space
 * as either mapA or mapB.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPBicMaps(
    UInt32				*bicMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	word, nWords;

    nWords = BitToWord( vPix * hPix );

    for( word = 0; word < nWords; ++word )
        *bicMap++ = *mapA++ & ~*mapB++;
}


/* BMAPPatchArea -----------------------------------------------------
 *
 * Return count of all one-bits within whole
 * words covering the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPPatchArea(
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pRow;
    UInt32	rowBytes, w;
    int		v, h, hL, hR, top, left, bot, right, area;

    area		= 0;
    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    rowBytes	= BitToByte( hPix );
    hL			= BitToWord( left );
    hR			= BitToWord( right - 1 );
    pRow		= (UInt32*)((char*)map + rowBytes * top);

    for( v = top; v < bot;
        ++v,
        pRow = (UInt32*)((char*)pRow + rowBytes) ) {

        for( h = hL; h <= hR; ++h ) {

            w = pRow[h];

            while( w ) {

                area += w & 1;
                w >>= 1;
            }
        }
    }

exit:
    return area;
}


/* BMAPZeroPatch -----------------------------------------------------
 *
 * Zero whole words that cover the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPZeroPatch(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst;
    UInt32	rowBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + rowBytes * top);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] = 0;
    }

exit:
    return;
}


/* BMAPFillPatch -----------------------------------------------------
 *
 * Set whole words that cover the given rectangular region to
 * 0xFFFFFFFF -- That is, set all patch pixels to ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPFillPatch(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst;
    UInt32	rowBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + rowBytes * top);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] = -1;
    }

exit:
    return;
}


/* BMAPCopyPatch -----------------------------------------------------
 *
 * Copy whole words that cover the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPCopyPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] = pSrc[h];
    }

exit:
    return;
}


/* BMAPInvertPatch ---------------------------------------------------
 *
 * Invert whole words that cover the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPInvertPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] = ~pSrc[h];
    }

exit:
    return;
}


/* BMAPOrPatch -------------------------------------------------------
 *
 * Logically OR whole words that cover the given rectangular region.
 * Result placed in dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPOrPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] |= pSrc[h];
    }

exit:
    return;
}


/* BMAPXorPatch ------------------------------------------------------
 *
 * Logically XOR whole words that cover the given rectangular region.
 * Result placed in dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPXorPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] ^= pSrc[h];
    }

exit:
    return;
}


/* BMAPXorNewPatch ---------------------------------------------------
 *
 * Logically XOR whole words that cover the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPXorNewPatch(
    UInt32				*xorMap,
    const UInt32		*mapA,
    const UInt32		*mapB,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc1, *pSrc2;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)xorMap	+ startBytes);
    pSrc1	= (UInt32*)((char*)mapA		+ startBytes);
    pSrc2	= (UInt32*)((char*)mapB		+ startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst  = (UInt32*)((char*)pDst + rowBytes),
        pSrc1 = (UInt32*)((char*)pSrc1 + rowBytes),
        pSrc2 = (UInt32*)((char*)pSrc2 + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] = pSrc1[h] ^ pSrc2[h];
    }

exit:
    return;
}


/* BMAPAndPatch ------------------------------------------------------
 *
 * dstMap becomes dstMap AND srcMap. Operation acts on whole
 * words covering the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPAndPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] &= pSrc[h];
    }

exit:
    return;
}


/* BMAPBicPatch ------------------------------------------------------
 *
 * dstMap becomes dstMap AND NOT(srcMap). Operation acts on whole
 * words covering the given rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPBicPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + startBytes);
    pSrc	= (UInt32*)((char*)srcMap + startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst = (UInt32*)((char*)pDst + rowBytes),
        pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] &= ~pSrc[h];
    }

exit:
    return;
}


/* BMAPBicDifPatch ---------------------------------------------------
 *
 * dstMap becomes dstMap AND NOT(mapA ^ mapB)...That is, the
 * difference between mapA and mapB is deleted from dstMap.
 *
 * Operation acts on whole words covering the given
 * rectangular region.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPBicDifPatch(
    UInt32				*dstMap,
    const UInt32		*mapA,
    const UInt32		*mapB,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc1, *pSrc2;
    UInt32	rowBytes, startBytes;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    startBytes	= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap	+ startBytes);
    pSrc1	= (UInt32*)((char*)mapA		+ startBytes);
    pSrc2	= (UInt32*)((char*)mapB		+ startBytes);
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    for( v = top; v < bot;
        ++v,
        pDst  = (UInt32*)((char*)pDst + rowBytes),
        pSrc1 = (UInt32*)((char*)pSrc1 + rowBytes),
        pSrc2 = (UInt32*)((char*)pSrc2 + rowBytes) ) {

        for( h = hL; h <= hR; ++h )
            pDst[h] &= ~(pSrc1[h] ^ pSrc2[h]);
    }

exit:
    return;
}


/* BMAPZeroBox -------------------------------------------------------
 *
 * Zero a rectangular region within dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPZeroBox(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst;
    UInt32	rowBytes, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + rowBytes * top);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );
    LMask	= ~LMask;
    RMask	= ~RMask;

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes) ) {

            pDst[hL] &= LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] = 0;

            pDst[hR] &= RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        LMask	|= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes) ) {

            *pDst &= LMask;
        }
    }

exit:
    return;
}


/* BMAPFillBox -------------------------------------------------------
 *
 * Fill a rectangular region within dstMap with ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPFillBox(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst;
    UInt32	rowBytes, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + rowBytes * top);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes) ) {

            pDst[hL] |= LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] = -1;

            pDst[hR] |= RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        LMask	&= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes) ) {

            *pDst |= LMask;
        }
    }

exit:
    return;
}


/* BMAPCopyBox -------------------------------------------------------
 *
 * Copy a rectangular region from srcMap to dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPCopyBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, W, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    W			= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + W);
    pSrc	= (UInt32*)((char*)srcMap + W);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            W			= pDst[hL];
            pDst[hL] 	= (W & ~LMask) | (pSrc[hL] & LMask);

            for( h = hL + 1; h < hR; ++h )
                pDst[h] = pSrc[h];

            W			= pDst[hR];
            pDst[hR]	= (W & ~RMask) | (pSrc[hR] & RMask);
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        pSrc	+= hL;
        LMask	&= RMask;
        RMask	= ~LMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            W		= *pDst;
            *pDst	= (W & RMask) | (*pSrc & LMask);
        }
    }

exit:
    return;
}


/* BMAPOrBox ---------------------------------------------------------
 *
 * Logically OR a rectangular region from srcMap into dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPOrBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, W, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    W			= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + W);
    pSrc	= (UInt32*)((char*)srcMap + W);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            pDst[hL] |= pSrc[hL] & LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] |= pSrc[h];

            pDst[hR] |= pSrc[hR] & RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        pSrc	+= hL;
        LMask	&= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            *pDst |= *pSrc & LMask;
        }
    }

exit:
    return;
}


/* BMAPXorBox --------------------------------------------------------
 *
 * Logically XOR a rectangular region from srcMap into dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPXorBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, W, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    W			= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + W);
    pSrc	= (UInt32*)((char*)srcMap + W);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            pDst[hL] ^= pSrc[hL] & LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] ^= pSrc[h];

            pDst[hR] ^= pSrc[hR] & RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        pSrc	+= hL;
        LMask	&= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            *pDst ^= *pSrc & LMask;
        }
    }

exit:
    return;
}


/* BMAPAndBox --------------------------------------------------------
 *
 * Logically AND a rectangular region from srcMap into dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPAndBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, W, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    W			= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + W);
    pSrc	= (UInt32*)((char*)srcMap + W);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );
    LMask	= ~LMask;
    RMask	= ~RMask;

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            pDst[hL] &= pSrc[hL] | LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] &= pSrc[h];

            pDst[hR] &= pSrc[hR] | RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        pSrc	+= hL;
        LMask	|= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            *pDst &= *pSrc | LMask;
        }
    }

exit:
    return;
}


/* BMAPBicBox --------------------------------------------------------
 *
 * Logically BIC a rectangular region from srcMap into dstMap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPBicBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*pDst, *pSrc;
    UInt32	rowBytes, W, LMask, RMask;
    int		v, h, hL, hR, top, left, bot, right;

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    W			= rowBytes * top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    pDst	= (UInt32*)((char*)dstMap + W);
    pSrc	= (UInt32*)((char*)srcMap + W);

    --right;
    hL		= BitToWord( left );
    hR		= BitToWord( right );
    LMask	= Mask32_ZerosLeftOf( left );
    RMask	= Mask32_ZerosRightOf( right );
    LMask	= ~LMask;
    RMask	= ~RMask;

    if( hL < hR ) {

        /* ------------------ */
        /* At least two words */
        /* ------------------ */

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            pDst[hL] &= ~pSrc[hL] | LMask;

            for( h = hL + 1; h < hR; ++h )
                pDst[h] &= ~pSrc[h];

            pDst[hR] &= ~pSrc[hR] | RMask;
        }
    }
    else {

        /* ------------------------ */
        /* All bits within one word */
        /* ------------------------ */

        pDst	+= hL;
        pSrc	+= hL;
        LMask	|= RMask;

        for( v = top; v < bot;
            ++v,
            pDst = (UInt32*)((char*)pDst + rowBytes),
            pSrc = (UInt32*)((char*)pSrc + rowBytes) ) {

            *pDst &= ~*pSrc | LMask;
        }
    }

exit:
    return;
}


/* BMAPInsetMap_1Pixel -----------------------------------------------
 *
 * The inset map is derived from the original by "killing"
 * any 1-pixel (object pixel) that is immediately adjacent
 * to a 0-pixel (background pixel). Adjacency is limited to
 * the horizontal and vertical directions.
 *
 * Horizontal killing is done by bitwise ANDing each word in the
 * original map with a right- and a left-shifted copy of itself.
 *
 * Vertical killing is done by ANDing each word with those above
 * and below.
 *
 * The original map is treated as if surrounded by a border of
 * zero-bits. This will always cause a white border to develop
 * around the image. In fact, we explicitly zero the entire top
 * and bottom rows.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * The storage space occupied by the two maps must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPInsetMap_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	*dst, *src1, *src2, *src3;
    UInt32	prevLSB, W;
    int		h, v, nWords;

/* --------------- */
/* Initializations */
/* --------------- */

    nWords	= BitToWord( hPix );
    dst		= dstMap;
    src1	= (UInt32*)srcMap;
    src2	= src1 + nWords;
    src3	= src2 + nWords;

/* ------------ */
/* Zero top row */
/* ------------ */

    for( h = 0; h < nWords; ++h )
        *dst++ = 0;

/* ------------- */
/* Interior rows */
/* ------------- */

    --nWords;

    for( v = 2; v < (int)vPix; ++v ) {

        /* ------------------------ */
        /* All but last word on row */
        /* ------------------------ */

        prevLSB = 0;

        for( h = 0; h < nWords; ++h ) {

            W		= *src2++;
            *dst++	= W
                        & ((W >> 1) | (prevLSB << AllButOneBit))
                        & ((W << 1) | (*src2 >> AllButOneBit))
                        & *src1++
                        & *src3++;
            prevLSB	= W & 1;
        }

        /* --------- */
        /* Last word */
        /* --------- */

        W		= *src2++;
        *dst++	= W
                    & ((W >> 1) | (prevLSB << AllButOneBit))
                    &  (W << 1)
                    & *src1++
                    & *src3++;
    }

/* --------------- */
/* Zero bottom row */
/* --------------- */

    for( h = 0; h <= nWords; ++h )
        *dst++ = 0;
}


/* BMAPInsetMap_NPixels ----------------------------------------------
 *
 * Call BMAPInsetMap_1Pixel nPixels times.
 *
 * If the returned value is zero, the result map is in dstMap0.
 * If the returned value is one,  the result map is in dstMap1.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * The storage spaces: srcMap, dstMap0, dstMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPInsetMap_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels )
{
    if( !nPixels-- ) {
        MEMCopyMap( dstMap1, srcMap, hPix, vPix );
        goto exit;
    }

    BMAPInsetMap_1Pixel( dstMap0, srcMap, hPix, vPix );

    while( nPixels >= 2 ) {
        BMAPInsetMap_1Pixel( dstMap1, dstMap0, hPix, vPix );
        BMAPInsetMap_1Pixel( dstMap0, dstMap1, hPix, vPix );
        nPixels -= 2;
    }

    if( nPixels )
        BMAPInsetMap_1Pixel( dstMap1, dstMap0, hPix, vPix );

exit:
    return (nPixels & 1);
}


/* BMAPTraceMap_1Pixel -----------------------------------------------
 *
 * The traced map is derived from the original by setting
 * any 0-pixel (background) that is immediately adjacent
 * to a 1-pixel (object pixel). Adjacency includes horizontal,
 * vertical and diagonal directions.
 *
 * Each word is bitwise ORed with left- and right-shifted copies of
 * itself, and with three copies each of the words above and below:
 * left-, un- and right-shifted.
 *
 * The original map is treated as if surrounded by a border of
 * zero-bits.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * The storage space occupied by the two maps must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPTraceMap_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	*dst, *src1, *src2, *src3;
    UInt32	W1, W2, W3, prvLSB1, prvLSB2, prvLSB3;
    int		v, h, nWords;

/* --------------- */
/* Initializations */
/* --------------- */

    nWords	= BitToWord( hPix );
    dst		= dstMap;
    src1	= (UInt32*)srcMap;
    src2	= src1;
    src3	= src1 + nWords;

    --nWords;

/* ------- */
/* Top row */
/* ------- */

    /* ------------------------ */
    /* All but last word on row */
    /* ------------------------ */

    prvLSB2 = 0;
    prvLSB3 = 0;

    for( h = 0; h < nWords; ++h ) {

        W2		= *src2++;
        W3		= *src3++;
        *dst++	= W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1) | (*src2 >> AllButOneBit)
                    | W3
                    | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                    | (W3 << 1) | (*src3 >> AllButOneBit);
        prvLSB2	= W2 & 1;
        prvLSB3	= W3 & 1;
    }

    /* --------- */
    /* Last word */
    /* --------- */

    W2		= *src2++;
    W3		= *src3++;
    *dst++	= W2
                | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                | (W2 << 1)
                | W3
                | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                | (W3 << 1);

/* ------------- */
/* Interior rows */
/* ------------- */

    for( v = 2; v < (int)vPix; ++v ) {

        /* ------------------------ */
        /* All but last word on row */
        /* ------------------------ */

        prvLSB1 = 0;
        prvLSB2 = 0;
        prvLSB3 = 0;

        for( h = 0; h < nWords; ++h ) {

            W1		= *src1++;
            W2		= *src2++;
            W3		= *src3++;
            *dst++	= W1
                        | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                        | (W1 << 1) | (*src1 >> AllButOneBit)
                        | W2
                        | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                        | (W2 << 1) | (*src2 >> AllButOneBit)
                        | W3
                        | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                        | (W3 << 1) | (*src3 >> AllButOneBit);
            prvLSB1	= W1 & 1;
            prvLSB2	= W2 & 1;
            prvLSB3	= W3 & 1;
        }

        /* --------- */
        /* Last word */
        /* --------- */

        W1		= *src1++;
        W2		= *src2++;
        W3		= *src3++;
        *dst++	= W1
                    | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                    | (W1 << 1)
                    | W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1)
                    | W3
                    | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                    | (W3 << 1);
    }

/* ---------- */
/* Bottom row */
/* ---------- */

    /* ------------------------ */
    /* All but last word on row */
    /* ------------------------ */

    prvLSB1 = 0;
    prvLSB2 = 0;

    for( h = 0; h < nWords; ++h ) {

        W1		= *src1++;
        W2		= *src2++;
        *dst++	= W1
                    | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                    | (W1 << 1) | (*src1 >> AllButOneBit)
                    | W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1) | (*src2 >> AllButOneBit);
        prvLSB1	= W1 & 1;
        prvLSB2	= W2 & 1;
    }

    /* --------- */
    /* Last word */
    /* --------- */

    W1		= *src1;
    W2		= *src2;
    *dst	= W1
                | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                | (W1 << 1)
                | W2
                | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                | (W2 << 1);
}


/* BMAPTraceMap_NPixels ----------------------------------------------
 *
 * Call BMAPTraceMap_1Pixel nPixels times.
 *
 * If the returned value is zero, the result map is in dstMap0.
 * If the returned value is one,  the result map is in dstMap1.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * The storage spaces: srcMap, dstMap0, dstMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPTraceMap_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels )
{
    if( !nPixels-- ) {
        MEMCopyMap( dstMap1, srcMap, hPix, vPix );
        goto exit;
    }

    BMAPTraceMap_1Pixel( dstMap0, srcMap, hPix, vPix );

    while( nPixels >= 2 ) {
        BMAPTraceMap_1Pixel( dstMap1, dstMap0, hPix, vPix );
        BMAPTraceMap_1Pixel( dstMap0, dstMap1, hPix, vPix );
        nPixels -= 2;
    }

    if( nPixels )
        BMAPTraceMap_1Pixel( dstMap1, dstMap0, hPix, vPix );

exit:
    return (nPixels & 1);
}


/* BMAPInsetPatch_1Pixel ---------------------------------------------
 *
 * Similar to BMAPInsetMap_1Pixel but acts on a submap defined by
 * the box parameter.
 *
 * The affected region is a patch covering all words overlapping
 * the box.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * (1) BEFORE calling this function the caller must clear the
 * destination patch, for example, with:
 *
 *		BMAPZeroPatch().
 *
 * (2) The storage space occupied by the two maps must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPInsetPatch_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*dst, *src1, *src2, *src3;
    UInt32	rowBytes, W, prevLSB;
    int		v, h, hL, hR, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    ++top;
    --bot;

    W		= rowBytes * top;
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    dst		= (UInt32*)((char*)dstMap + W);
    src2	= (UInt32*)((char*)srcMap + W);
    src1	= (UInt32*)((char*)src2 - rowBytes);
    src3	= (UInt32*)((char*)src2 + rowBytes);

/* ------------- */
/* Interior rows */
/* ------------- */

    for( v = top; v < bot;
        ++v,
        dst  = (UInt32*)((char*)dst + rowBytes),
        src1 = (UInt32*)((char*)src1 + rowBytes),
        src2 = (UInt32*)((char*)src2 + rowBytes),
        src3 = (UInt32*)((char*)src3 + rowBytes) ) {

        /* ------------------------ */
        /* All but last word on row */
        /* ------------------------ */

        prevLSB = 0;

        for( h = hL; h < hR; ++h ) {

            W		= src2[h];
            dst[h]	= W
                        & ((W >> 1) | (prevLSB << AllButOneBit))
                        & ((W << 1) | (src2[h + 1] >> AllButOneBit))
                        & src1[h]
                        & src3[h];
            prevLSB	= W & 1;
        }

        /* --------- */
        /* Last word */
        /* --------- */

        W		= src2[hR];
        dst[hR]	= W
                    & ((W >> 1) | (prevLSB << AllButOneBit))
                    &  (W << 1)
                    & src1[hR]
                    & src3[hR];
    }

exit:
    return;
}


/* BMAPInsetPatch_NPixels --------------------------------------------
 *
 * Call BMAPInsetPatch_1Pixel nPixels times.
 *
 * If the returned value is zero, the result map is in dstMap0.
 * If the returned value is one,  the result map is in dstMap1.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * (1) BEFORE calling this function the caller must clear the
 * destination patches, for example, with:
 *
 *		BMAPZeroPatch().
 *
 * (2) The storage spaces: srcMap, dstMap0, dstMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPInsetPatch_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels,
    const U16BoxPtr		box )
{
    if( !nPixels-- ) {
        BMAPCopyPatch( dstMap1, srcMap, hPix, vPix, box );
        goto exit;
    }

    BMAPInsetPatch_1Pixel( dstMap0, srcMap, hPix, vPix, box );

    while( nPixels >= 2 ) {
        BMAPInsetPatch_1Pixel( dstMap1, dstMap0, hPix, vPix, box );
        BMAPInsetPatch_1Pixel( dstMap0, dstMap1, hPix, vPix, box );
        nPixels -= 2;
    }

    if( nPixels )
        BMAPInsetPatch_1Pixel( dstMap1, dstMap0, hPix, vPix, box );

exit:
    return (nPixels & 1);
}


/* BMAPTracePatch_1Pixel ---------------------------------------------
 *
 * Similar to BMAPTraceMap_1Pixel but acts on a submap defined by
 * the box parameter.
 *
 * The affected region is a patch covering all words overlapping
 * the dilated box (see notes below).
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
 * origBox, then the caller must calculate the expanded box
 * (to be used by this function) as follows:
 *
 *		GEOMOutsetU16Box( &box, &origBox, hPix, vPix, 1 ).
 *
 * (1b) The patch containing the source object must correspond to
 * the expanded box parameter. That is, the source object should
 * be isolated into a patch (or map) that has been pre-zeroed out
 * to the larger size.
 *
 * (2) The storage space occupied by the two maps must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPTracePatch_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box )
{
    UInt32	*dst, *src1, *src2, *src3;
    UInt32	rowBytes, W1, W2, W3, prvLSB1, prvLSB2, prvLSB3;
    int		v, h, hL, hR, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    rowBytes	= BitToByte( hPix );

    top			= box->top;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    --bot;

    W1		= rowBytes * top;
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    dst		= (UInt32*)((char*)dstMap + W1);
    src1	= (UInt32*)((char*)srcMap + W1);
    src2	= (UInt32*)((char*)src1 + rowBytes);
    src3	= (UInt32*)((char*)src2 + rowBytes);

/* ------- */
/* Top row */
/* ------- */

    /* ------------------------ */
    /* All but last word on row */
    /* ------------------------ */

    prvLSB1 = 0;
    prvLSB2 = 0;

    for( h = hL; h < hR; ++h ) {

        W1		= src1[h];
        W2		= src2[h];
        dst[h]	= W1
                    | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                    | (W1 << 1) | (src1[h + 1] >> AllButOneBit)
                    | W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1) | (src2[h + 1] >> AllButOneBit);
        prvLSB1	= W1 & 1;
        prvLSB2	= W2 & 1;
    }

    /* --------- */
    /* Last word */
    /* --------- */

    W1		= src1[hR];
    W2		= src2[hR];
    dst[hR]	= W1
                | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                | (W1 << 1)
                | W2
                | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                | (W2 << 1);

    dst = (UInt32*)((char*)dst + rowBytes);

/* ------------- */
/* Interior rows */
/* ------------- */

    for( v = top + 1; v < bot; ++v,
        dst  = (UInt32*)((char*)dst + rowBytes),
        src1 = (UInt32*)((char*)src1 + rowBytes),
        src2 = (UInt32*)((char*)src2 + rowBytes),
        src3 = (UInt32*)((char*)src3 + rowBytes) ) {

        /* ------------------------ */
        /* All but last word on row */
        /* ------------------------ */

        prvLSB1 = 0;
        prvLSB2 = 0;
        prvLSB3 = 0;

        for( h = hL; h < hR; ++h ) {

            W1		= src1[h];
            W2		= src2[h];
            W3		= src3[h];
            dst[h]	= W1
                        | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                        | (W1 << 1) | (src1[h + 1] >> AllButOneBit)
                        | W2
                        | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                        | (W2 << 1) | (src2[h + 1] >> AllButOneBit)
                        | W3
                        | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                        | (W3 << 1) | (src3[h + 1] >> AllButOneBit);
            prvLSB1	= W1 & 1;
            prvLSB2	= W2 & 1;
            prvLSB3	= W3 & 1;
        }

        /* --------- */
        /* Last word */
        /* --------- */

        W1		= src1[hR];
        W2		= src2[hR];
        W3		= src3[hR];
        dst[hR]	= W1
                    | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                    | (W1 << 1)
                    | W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1)
                    | W3
                    | (W3 >> 1) | (prvLSB3 << AllButOneBit)
                    | (W3 << 1);
    }

/* ---------- */
/* Bottom row */
/* ---------- */

    /* ------------------------ */
    /* All but last word on row */
    /* ------------------------ */

    prvLSB1 = 0;
    prvLSB2 = 0;

    for( h = hL; h < hR; ++h ) {

        W1		= src1[h];
        W2		= src2[h];
        dst[h]	= W1
                    | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                    | (W1 << 1) | (src1[h + 1] >> AllButOneBit)
                    | W2
                    | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                    | (W2 << 1) | (src2[h + 1] >> AllButOneBit);
        prvLSB1	= W1 & 1;
        prvLSB2	= W2 & 1;
    }

    /* --------- */
    /* Last word */
    /* --------- */

    W1		= src1[hR];
    W2		= src2[hR];
    dst[hR]	= W1
                | (W1 >> 1) | (prvLSB1 << AllButOneBit)
                | (W1 << 1)
                | W2
                | (W2 >> 1) | (prvLSB2 << AllButOneBit)
                | (W2 << 1);

exit:
    return;
}


/* BMAPTracePatch_NPixels --------------------------------------------
 *
 * Call BMAPTracePatch_1Pixel nPixels times.
 *
 * If the returned value is zero, the result map is in dstMap0.
 * If the returned value is one,  the result map is in dstMap1.
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
 * origBox, then the caller must calculate the expanded box
 * (to be used by this function) as follows:
 *
 *		GEOMOutsetU16Box( &box, &origBox, hPix, vPix, nPixels ).
 *
 * (1b) The patch containing the source object must correspond to
 * the expanded box parameter. That is, the source object should
 * be isolated into a patch (or map) that has been pre-zeroed out
 * to the larger size.
 *
 * (2) The storage spaces: srcMap, dstMap0, dstMap1 must not overlap.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPTracePatch_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels,
    const U16BoxPtr		box )
{
    if( !nPixels-- ) {
        BMAPCopyPatch( dstMap1, srcMap, hPix, vPix, box );
        goto exit;
    }

    BMAPTracePatch_1Pixel( dstMap0, srcMap, hPix, vPix, box );

    while( nPixels >= 2 ) {
        BMAPTracePatch_1Pixel( dstMap1, dstMap0, hPix, vPix, box );
        BMAPTracePatch_1Pixel( dstMap0, dstMap1, hPix, vPix, box );
        nPixels -= 2;
    }

    if( nPixels )
        BMAPTracePatch_1Pixel( dstMap1, dstMap0, hPix, vPix, box );

exit:
    return (nPixels & 1);
}


/* BMAPGetSeedPatch --------------------------------------------------
 *
 * Scan map for coordinates of any seed point and return true
 * if found.
 *
 * Patch is scanned from top to bottom, left to right. The first
 * black pixel found terminates the search.
 *
 * The top of the search box is specified separately, allowing
 * caller to make faster repeated searches.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int BMAPGetSeedPatch(
    UInt32				*hSeed,
    UInt32				*vSeed,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box,
    UInt32				top )
{
    UInt32	*pRow;
    UInt32	rowBytes, w;
    int		v, h, hL, hR, left, bot, right, found;

    found		= false;
    left		= box->left;
    bot			= box->bottom;
    right		= box->right;
    *hSeed		= left;
    *vSeed		= top;

    if( bot > (int)vPix )
        bot = vPix;

    if( right > (int)hPix )
        right = hPix;

    if( (int)top >= bot )
        goto exit;

    if( left >= right )
        goto exit;

    rowBytes	= BitToByte( hPix );
    hL			= BitToWord( left );
    hR			= BitToWord( right - 1 );
    pRow		= (UInt32*)((char*)map + rowBytes * top);

    for( v = top; v < bot;
        ++v, pRow = (UInt32*)((char*)pRow + rowBytes) ) {

        for( h = hL; h <= hR; ++h ) {

            if( w = pRow[h] ) {

                *vSeed = v;

                for( v = 0; v < WordBits; ++v, w <<= 1 ) {

                    if( IsHiBitSet32( w ) ) {

                        *hSeed	= h * WordBits + v;
                        found	= true;
                        goto exit;
                    }
                }
            }
        }
    }

exit:
    return found;
}


/* BMAPBucketTool ----------------------------------------------------
 *
 * If the pixel (hSeed, vSeed) is set:
 *
 * - Clear it...
 *
 * - Continue to clear all touching pixels (immediate neighbors
 * in the vertical and horizontal directions). The set of touching
 * pixels is referred to as a "region" in these comments.
 *
 * If (hSeed, vSeed) is already clear, nothing is done. The region
 * is said to be empty in this case.
 *
 * - rgnBox
 * On exit, rgnBox is set to "enclose" the seeded region:
 *
 *		top		= top-most v-coordinate in region
 *		left	= left-most h-coordinate
 *		bottom	= bottom-most + 1
 *		right	= right-most + 1
 *
 * The addition of one to the bottom and right of the box follows
 * common graphics convention. This simplifies the calculation of
 * region width as (right - left) and height as (bottom - top).
 *
 * If the region was empty, the box will be empty. That is:
 *
 *		bottom <= top, or, right <= left.
 *
 * Return "area" (count of pixels in seeded region), or,
 * return zero if region empty.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * Object extent is constrained only by the edges of the full map.
 * If using this tool on a subregion, the caller must make sure
 * the subregion is isolated.
 *
 * - stackBase:
 * BMAPBucketTool needs a large scratch area in which to store
 * intermediate results. The actual amount of storage needed
 * depends upon the shape and size of the seeded region, but is
 * roughly: sizeof(BMAPBucketRec) * (number of bits in region).
 * It is safest to make sure there is adequate space by taking
 * the worst-case value to be:
 *
 *		hPix * vPix * sizeof(BMAPBucketRec).
 *
 * In the current implementation, sizeof(BMAPBucketRec) = 4 bytes,
 * but you should use the sizeof macro for future compatibility.
 * BMAPBucketRec is defined in the header file.
 *
 * -----
 * Notes
 * -----
 *
 * Each pixel that is determined to belong to the region is
 * pushed onto the stack. Each pixel that is popped off the
 * stack is examined as a neighborhood center.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMAPBucketTool(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed )
{
    BMAPBucketPtr	stack;
    UInt32			*thisRow;
    UInt32			rowBytes, w, mask, area;
    int				hIndx, top, left, bot, right;

    rowBytes	= BitToByte( hPix );
    hIndx		= BitToWord( hSeed );
    mask		= BMAP_HToMask( hSeed );
    thisRow		= (UInt32*)((char*)map + rowBytes * vSeed);
    stack		= stackBase;
    w			= thisRow[hIndx];
    top			= vSeed;
    left		= hSeed;
    bot			= -1;
    right		= -1;
    area		= 0;
    --hPix;
    --vPix;

    if( w & mask ) {

        thisRow[hIndx]	= w ^ mask;
        bot				= vSeed;
        right			= hSeed;
        area			= 1;

        goto start;

        do {

            /* --- */
            /* Pop */
            /* --- */

            --stack;
            ++area;
            vSeed	= stack->vSeed;
            hSeed	= stack->hSeed;
            thisRow	= (UInt32*)((char*)map + rowBytes * vSeed);

            /* ----- */
            /* Above */
            /* ----- */
start:
            if( vSeed ) {

                thisRow = (UInt32*)((char*)thisRow - rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed - 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed <= top )
                        top = vSeed - 1;
                }

                thisRow = (UInt32*)((char*)thisRow + rowBytes);
            }

            /* ---- */
            /* Left */
            /* ---- */

            if( hSeed ) {

                --hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed < left )
                        left = hSeed;
                }

                ++hSeed;
            }

            /* ----- */
            /* Right */
            /* ----- */

            if( hSeed < hPix ) {

                ++hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed > right )
                        right = hSeed;
                }

                --hSeed;
            }

            /* ----- */
            /* Below */
            /* ----- */

            if( vSeed < vPix ) {

                thisRow = (UInt32*)((char*)thisRow + rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed + 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed >= bot )
                        bot = vSeed + 1;
                }
            }

        } while( stack > stackBase );
    }

    rgnBox->top		= top;
    rgnBox->left	= left;
    rgnBox->bottom	= bot + 1;
    rgnBox->right	= right + 1;

    return area;
}


/* BMAPBucketTool8Way ------------------------------------------------
 *
 * This function works exactly like BMAPBucketTool except that
 * all eight neighbors are considered.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMAPBucketTool8Way(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed )
{
    BMAPBucketPtr	stack;
    UInt32			*thisRow;
    UInt32			rowBytes, w, mask, area;
    int				hIndx, top, left, bot, right;

    rowBytes	= BitToByte( hPix );
    hIndx		= BitToWord( hSeed );
    mask		= BMAP_HToMask( hSeed );
    thisRow		= (UInt32*)((char*)map + rowBytes * vSeed);
    stack		= stackBase;
    w			= thisRow[hIndx];
    top			= vSeed;
    left		= hSeed;
    bot			= -1;
    right		= -1;
    area		= 0;
    --hPix;
    --vPix;

    if( w & mask ) {

        thisRow[hIndx]	= w ^ mask;
        bot				= vSeed;
        right			= hSeed;
        area			= 1;

        goto start;

        do {

            /* --- */
            /* Pop */
            /* --- */

            --stack;
            ++area;
            vSeed	= stack->vSeed;
            hSeed	= stack->hSeed;
            thisRow	= (UInt32*)((char*)map + rowBytes * vSeed);

            /* ----- */
            /* Above */
            /* ----- */
start:
            if( vSeed ) {

                thisRow = (UInt32*)((char*)thisRow - rowBytes);
                --vSeed;

                /* ------------ */
                /* Above - left */
                /* ------------ */

                if( hSeed ) {

                    --hSeed;

                    hIndx	= BitToWord( hSeed );
                    mask	= BMAP_HToMask( hSeed );
                    w		= thisRow[hIndx];

                    if( w & mask ) {	/* clear and push */

                        thisRow[hIndx]	= w ^ mask;
                        stack->vSeed	= vSeed;
                        stack->hSeed	= hSeed;
                        ++stack;

                        if( (int)hSeed < left )
                            left = hSeed;

                        if( (int)vSeed < top )
                            top = vSeed;
                    }

                    ++hSeed;
                }

                /* ----------- */
                /* Above - mid */
                /* ----------- */

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed < top )
                        top = vSeed;
                }

                /* ------------- */
                /* Above - right */
                /* ------------- */

                if( hSeed < hPix ) {

                    ++hSeed;

                    hIndx	= BitToWord( hSeed );
                    mask	= BMAP_HToMask( hSeed );
                    w		= thisRow[hIndx];

                    if( w & mask ) {	/* clear and push */

                        thisRow[hIndx]	= w ^ mask;
                        stack->vSeed	= vSeed;
                        stack->hSeed	= hSeed;
                        ++stack;

                        if( (int)hSeed > right )
                            right = hSeed;

                        if( (int)vSeed < top )
                            top = vSeed;
                    }

                    --hSeed;
                }

                ++vSeed;
                thisRow = (UInt32*)((char*)thisRow + rowBytes);
            }

            /* ---- */
            /* Left */
            /* ---- */

            if( hSeed ) {

                --hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed < left )
                        left = hSeed;
                }

                ++hSeed;
            }

            /* ----- */
            /* Right */
            /* ----- */

            if( hSeed < hPix ) {

                ++hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed > right )
                        right = hSeed;
                }

                --hSeed;
            }

            /* ----- */
            /* Below */
            /* ----- */

            if( vSeed < vPix ) {

                thisRow = (UInt32*)((char*)thisRow + rowBytes);
                ++vSeed;

                /* ------------ */
                /* Below - left */
                /* ------------ */

                if( hSeed ) {

                    --hSeed;

                    hIndx	= BitToWord( hSeed );
                    mask	= BMAP_HToMask( hSeed );
                    w		= thisRow[hIndx];

                    if( w & mask ) {	/* clear and push */

                        thisRow[hIndx]	= w ^ mask;
                        stack->vSeed	= vSeed;
                        stack->hSeed	= hSeed;
                        ++stack;

                        if( (int)hSeed < left )
                            left = hSeed;

                        if( (int)vSeed > bot )
                            bot = vSeed;
                    }

                    ++hSeed;
                }

                /* ----------- */
                /* Below - mid */
                /* ----------- */

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed > bot )
                        bot = vSeed;
                }

                /* ------------- */
                /* Below - right */
                /* ------------- */

                if( hSeed < hPix ) {

                    ++hSeed;

                    hIndx	= BitToWord( hSeed );
                    mask	= BMAP_HToMask( hSeed );
                    w		= thisRow[hIndx];

                    if( w & mask ) {	/* clear and push */

                        thisRow[hIndx]	= w ^ mask;
                        stack->vSeed	= vSeed;
                        stack->hSeed	= hSeed;
                        ++stack;

                        if( (int)hSeed > right )
                            right = hSeed;

                        if( (int)vSeed > bot )
                            bot = vSeed;
                    }

                    --hSeed;
                }
            }

        } while( stack > stackBase );
    }

    rgnBox->top		= top;
    rgnBox->left	= left;
    rgnBox->bottom	= bot + 1;
    rgnBox->right	= right + 1;

    return area;
}


/* BMAPBoxedBucketTool -----------------------------------------------
 *
 * This function works exactly like BMAPBucketTool with the addition
 * of being constrained to examine only pixels contained within
 * the given box.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * - Your constraining box is assumed to lie within the map.
 *
 * - stackBase:
 * BMAPBucketTool needs a large scratch area in which to store
 * intermediate results. The actual amount of storage needed
 * depends upon the shape and size of the seeded region, but is
 * roughly: sizeof(BMAPBucketRec) * (number of bits in region).
 * It is safest to make sure there is adequate space by taking
 * the worst-case value to be:
 *
 *		hPix * vPix * sizeof(BMAPBucketRec).
 *
 * In the current implementation, sizeof(BMAPBucketRec) = 4 bytes,
 * but you should use the sizeof macro for future compatibility.
 * BMAPBucketRec is defined in the header file.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMAPBoxedBucketTool(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed,
    const U16BoxPtr		box )
{
    BMAPBucketPtr	stack;
    UInt32			*thisRow;
    UInt32			rowBytes, w, mask, area;
    int				hIndx, top, left, bot, right;

    rowBytes	= BitToByte( hPix );
    hIndx		= BitToWord( hSeed );
    mask		= BMAP_HToMask( hSeed );
    thisRow		= (UInt32*)((char*)map + rowBytes * vSeed);
    stack		= stackBase;
    w			= thisRow[hIndx];
    top			= vSeed;
    left		= hSeed;
    bot			= -1;
    right		= -1;
    area		= 0;
    hPix		= box->right  - 1;
    vPix		= box->bottom - 1;

    if( w & mask ) {

        thisRow[hIndx]	= w ^ mask;
        bot				= vSeed;
        right			= hSeed;
        area			= 1;

        goto start;

        do {

            /* --- */
            /* Pop */
            /* --- */

            --stack;
            ++area;
            vSeed	= stack->vSeed;
            hSeed	= stack->hSeed;
            thisRow	= (UInt32*)((char*)map + rowBytes * vSeed);

            /* ----- */
            /* Above */
            /* ----- */
start:
            if( vSeed > box->top ) {

                thisRow = (UInt32*)((char*)thisRow - rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed - 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed <= top )
                        top = vSeed - 1;
                }

                thisRow = (UInt32*)((char*)thisRow + rowBytes);
            }

            /* ---- */
            /* Left */
            /* ---- */

            if( hSeed > box->left ) {

                --hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed < left )
                        left = hSeed;
                }

                ++hSeed;
            }

            /* ----- */
            /* Right */
            /* ----- */

            if( hSeed < hPix ) {

                ++hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed > right )
                        right = hSeed;
                }

                --hSeed;
            }

            /* ----- */
            /* Below */
            /* ----- */

            if( vSeed < vPix ) {

                thisRow = (UInt32*)((char*)thisRow + rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed + 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed >= bot )
                        bot = vSeed + 1;
                }
            }

        } while( stack > stackBase );
    }

    rgnBox->top		= top;
    rgnBox->left	= left;
    rgnBox->bottom	= bot + 1;
    rgnBox->right	= right + 1;

    return area;
}


/* BMAPBoxCGArea -----------------------------------------------------
 *
 * This function works exactly like BMAPBucketTool with the addition
 * of calculating the centroid of the region.
 *
 * Each pixel is weighted equally in the centroid calculation.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * Object extent is constrained only by the edges of the full map.
 * If using this tool on a subregion, the caller must make sure
 * the subregion is isolated.
 *
 * - stackBase:
 * BMAPBucketTool needs a large scratch area in which to store
 * intermediate results. The actual amount of storage needed
 * depends upon the shape and size of the seeded region, but is
 * roughly: sizeof(BMAPBucketRec) * (number of bits in region).
 * It is safest to make sure there is adequate space by taking
 * the worst-case value to be:
 *
 *		hPix * vPix * sizeof(BMAPBucketRec).
 *
 * In the current implementation, sizeof(BMAPBucketRec) = 4 bytes,
 * but you should use the sizeof macro for future compatibility.
 * BMAPBucketRec is defined in the header file.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 BMAPBoxCGArea(
    BMAPBoxCGPtr		rgn,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed )
{
    BMAPBucketPtr	stack;
    UInt32			*thisRow;
    UInt32			rowBytes, w, mask, area, xC, yC;
    int				hIndx, top, left, bot, right;

    rowBytes	= BitToByte( hPix );
    hIndx		= BitToWord( hSeed );
    mask		= BMAP_HToMask( hSeed );
    thisRow		= (UInt32*)((char*)map + rowBytes * vSeed);
    stack		= stackBase;
    w			= thisRow[hIndx];
    top			= vSeed;
    left		= hSeed;
    bot			= -1;
    right		= -1;
    xC			= 0;
    yC			= 0;
    area		= 0;
    --hPix;
    --vPix;

    if( w & mask ) {

        thisRow[hIndx]	= w ^ mask;
        yC				= vSeed;
        bot				= vSeed;
        xC				= hSeed;
        right			= hSeed;
        area			= 1;

        goto start;

        do {

            /* --- */
            /* Pop */
            /* --- */

            --stack;
            ++area;
            vSeed	= stack->vSeed;
            hSeed	= stack->hSeed;
            thisRow	= (UInt32*)((char*)map + rowBytes * vSeed);

            yC += vSeed;
            xC += hSeed;

            /* ----- */
            /* Above */
            /* ----- */
start:
            if( vSeed ) {

                thisRow = (UInt32*)((char*)thisRow - rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed - 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed <= top )
                        top = vSeed - 1;
                }

                thisRow = (UInt32*)((char*)thisRow + rowBytes);
            }

            /* ---- */
            /* Left */
            /* ---- */

            if( hSeed ) {

                --hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed < left )
                        left = hSeed;
                }

                ++hSeed;
            }

            /* ----- */
            /* Right */
            /* ----- */

            if( hSeed < hPix ) {

                ++hSeed;

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)hSeed > right )
                        right = hSeed;
                }

                --hSeed;
            }

            /* ----- */
            /* Below */
            /* ----- */

            if( vSeed < vPix ) {

                thisRow = (UInt32*)((char*)thisRow + rowBytes);

                hIndx	= BitToWord( hSeed );
                mask	= BMAP_HToMask( hSeed );
                w		= thisRow[hIndx];

                if( w & mask ) {	/* clear and push */

                    thisRow[hIndx]	= w ^ mask;
                    stack->vSeed	= vSeed + 1;
                    stack->hSeed	= hSeed;
                    ++stack;

                    if( (int)vSeed >= bot )
                        bot = vSeed + 1;
                }
            }

        } while( stack > stackBase );
    }

    rgn->box.top	= top;
    rgn->box.left	= left;
    rgn->box.bottom	= bot + 1;
    rgn->box.right	= right + 1;

    if( area ) {
        rgn->centroid.x = (int)(xC / area);
        rgn->centroid.y = (int)(yC / area);
    }

    return area;
}


