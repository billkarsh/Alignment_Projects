/* BK_SUM ------------------------------------------------------------
 *
 * Sums over masked images.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_SUM.h"






/* SUMMap_F_A --------------------------------------------------------
 *
 * Scan entire map.
 *
 * For each 1 in the map, accumulate the corresponding image
 * intensity into F and increment A.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SUMMap_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	*pPair;
    UInt32	pair, word, nWords, W, F, A;
    int		pairIdx;

/* --------------- */
/* Initializations */
/* --------------- */

    F = 0;
    A = 0;

    nWords	= BitToWord( vPix * hPix );
    pPair	= (UInt32*)data16Bit;

/* ----------------------------- */
/* Loop over 32-bit words in map */
/* ----------------------------- */

    for( word = 0; word < nWords; ++word, pPair += HalfBits ) {

        if( W = *map++ ) {

            /* --------------------------- */
            /* Loop over bit pairs in word */
            /* --------------------------- */

            for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                UInt32	M;

                pair = pPair[pairIdx];

                M  = SignExtHiBit32( W );
                F += pair & LoMask & M;
                A += M & 1;

                W <<= 1;

                M  = SignExtHiBit32( W );
                F += (pair >> HalfBits) & M;
                A += M & 1;

                W <<= 1;
            }
        }
    }

    sum->F = F;
    sum->A = A;
}


/* SUMMap_Obj_Bkg ----------------------------------------------------
 *
 * Scan entire map.
 *
 * For each 1 in the object map, accumulate the corresponding image
 * intensity into obj.F and increment obj.A. bkg.F and bkg.A are
 * similarly summed using the background map.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SUMMap_Obj_Bkg(
    SUM_Obj_Bkg_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*objMap,
    const UInt32		*bkgMap,
    UInt32				hPix,
    UInt32				vPix )
{
    UInt32	*pPair;
    UInt32	pair, word, nWords, obW, bgW;
    int		pairIdx;

/* --------------- */
/* Initializations */
/* --------------- */

    sum->obj.F	= 0;
    sum->obj.A	= 0;
    sum->bkg.F	= 0;
    sum->bkg.A	= 0;

    nWords		= BitToWord( vPix * hPix );
    pPair		= (UInt32*)data16Bit;

/* ------------------------------ */
/* Loop over 32-bit words in maps */
/* ------------------------------ */

    for( word = 0; word < nWords; ++word, pPair += HalfBits ) {

        obW	= *objMap++;
        bgW	= *bkgMap++;

        if( obW | bgW ) {

            /* --------------------------- */
            /* Loop over bit pairs in word */
            /* --------------------------- */

            for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                pair = pPair[pairIdx];

                if( IsHiBitSet32( obW ) ) {
                    sum->obj.F += pair & LoMask;
                    ++sum->obj.A;
                }
                else if( IsHiBitSet32( bgW ) ) {
                    sum->bkg.F += pair & LoMask;
                    ++sum->bkg.A;
                }

                obW <<= 1;
                bgW <<= 1;

                if( IsHiBitSet32( obW ) ) {
                    sum->obj.F += pair >> HalfBits;
                    ++sum->obj.A;
                }
                else if( IsHiBitSet32( bgW ) ) {
                    sum->bkg.F += pair >> HalfBits;
                    ++sum->bkg.A;
                }

                obW <<= 1;
                bgW <<= 1;
            }
        }
    }
}


/* SUMBlob_F_A -------------------------------------------------------
 *
 * Scan the patch covering all words overlapping rBlob.
 *
 * For each 1 in the map, accumulate the corresponding image
 * intensity into F and increment A.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SUMBlob_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, W, F, A;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    F = 0;
    A = 0;

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

                    UInt32	M;

                    pair = pPair[pairIdx];

                    M  = SignExtHiBit32( W );
                    F += pair & LoMask & M;
                    A += M & 1;

                    W <<= 1;

                    M  = SignExtHiBit32( W );
                    F += (pair >> HalfBits) & M;
                    A += M & 1;

                    W <<= 1;
                }
            }
        }
    }

exit:
    sum->F = F;
    sum->A = A;
}


/* SUMBlob_F_A_Ipk ---------------------------------------------------
 *
 * Scan the patch covering all words overlapping rBlob.
 *
 * For each 1 in the map, accumulate the corresponding image
 * intensity into F and increment A. Also, set Ipk to the
 * largest intensity found in the patch.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SUMBlob_F_A_Ipk(
    SUM_F_A_Ipk_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, W, F, A, Ipk, I;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    F	= 0;
    Ipk	= 0;
    A	= 0;

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

                        I = pair & LoMask;
                        ++A;
                        F += I;

                        if( I > Ipk )
                            Ipk = I;
                    }

                    W <<= 1;

                    if( IsHiBitSet32( W ) ) {

                        I = pair >> HalfBits;
                        ++A;
                        F += I;

                        if( I > Ipk )
                            Ipk = I;
                    }

                    W <<= 1;
                }
            }
        }
    }

exit:
    sum->F		= F;
    sum->Ipk	= Ipk;
    sum->A		= A;
}


/* SUMBlob_Obj_Bkg ---------------------------------------------------
 *
 * Scan the patch covering all words overlapping rBlob.
 *
 * For each 1 in the object map, accumulate the corresponding image
 * intensity into obj.F and increment obj.A. bkg.F and bkg.A are
 * similarly summed using the background map.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SUMBlob_Obj_Bkg(
    SUM_Obj_Bkg_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*objMap,
    const UInt32		*bkgMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob )
{
    UInt32	*pPair;
    UInt32	mapRowBytes, imgRowBytes, pair, obW, bgW;
    int		v, h, hL, hR, pairIdx,
            top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

    sum->obj.F	= 0;
    sum->obj.A	= 0;
    sum->bkg.F	= 0;
    sum->bkg.A	= 0;

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

    pair	= mapRowBytes * top;
    hL		= BitToWord( left );
    hR		= BitToWord( right - 1 );

    objMap		= (UInt32*)((char*)objMap + pair);
    bkgMap		= (UInt32*)((char*)bkgMap + pair);
    data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes * top
                    + hL * WordBits * sizeof(UInt16));

/* --------------------------------- */
/* Loop over 32-bit words in patches */
/* --------------------------------- */

    for( v = top; v < bot; ++v,
        data16Bit	= (UInt16*)((char*)data16Bit + imgRowBytes),
        objMap		= (UInt32*)((char*)objMap + mapRowBytes),
        bkgMap		= (UInt32*)((char*)bkgMap + mapRowBytes) ) {

        pPair = (UInt32*)data16Bit;

        for( h = hL; h <= hR; ++h, pPair += HalfBits ) {

            obW	= objMap[h];
            bgW	= bkgMap[h];

            if( obW | bgW ) {

                /* --------------------------- */
                /* Loop over bit pairs in word */
                /* --------------------------- */

                for( pairIdx = 0; pairIdx < HalfBits; ++pairIdx ) {

                    pair = pPair[pairIdx];

                    if( IsHiBitSet32( obW ) ) {
                        sum->obj.F += pair & LoMask;
                        ++sum->obj.A;
                    }
                    else if( IsHiBitSet32( bgW ) ) {
                        sum->bkg.F += pair & LoMask;
                        ++sum->bkg.A;
                    }

                    obW <<= 1;
                    bgW <<= 1;

                    if( IsHiBitSet32( obW ) ) {
                        sum->obj.F += pair >> HalfBits;
                        ++sum->obj.A;
                    }
                    else if( IsHiBitSet32( bgW ) ) {
                        sum->bkg.F += pair >> HalfBits;
                        ++sum->bkg.A;
                    }

                    obW <<= 1;
                    bgW <<= 1;
                }
            }
        }
    }

exit:
    return;
}


