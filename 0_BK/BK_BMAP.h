/* BK_BMAP.h ---------------------------------------------------------
 *
 * Bitmap operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#include	"BK_DEFS_GEN.h"


#if defined(__cplusplus)
extern "C" {
#endif


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	BMAP_HToMask( h )										\
    (HiBit >> WordRem( h ))

#define	BMAP_IsRowBit( pRow, h )								\
    ((pRow[BitToWord( h )] & BMAP_HToMask( h )) != 0)

/* --------------------------------------------------------------- */
/* Data structures ----------------------------------------------- */
/* --------------------------------------------------------------- */

#if defined(WIN32)
#pragma pack( push, PACK_S )
#pragma pack( 1 )
#elif defined(__linux__)
#pragma pack( push, 1 )
#pragma pack( 1 )
#endif

typedef struct {

    U16Box		box;
    Centroid	centroid;

} BMAPBoxCGRec, *BMAPBoxCGPtr;

typedef struct {

    UInt16		vSeed,
                hSeed;

} BMAPBucketRec, *BMAPBucketPtr;

#if defined(WIN32)
#pragma pack( pop, PACK_S )
#elif defined(__linux__)
#pragma pack( pop )
#endif






void BMAPOrMaps(
    UInt32				*orMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix );

void BMAPXorMaps(
    UInt32				*xorMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix );

void BMAPAndMaps(
    UInt32				*andMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix );

void BMAPBicMaps(
    UInt32				*bicMap,
    UInt32				*mapA,
    UInt32				*mapB,
    UInt32				hPix,
    UInt32				vPix );

int BMAPPatchArea(
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPZeroPatch(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPFillPatch(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPCopyPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPInvertPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPOrPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPXorPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPXorNewPatch(
    UInt32				*xorMap,
    const UInt32		*mapA,
    const UInt32		*mapB,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPAndPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPBicPatch(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPBicDifPatch(
    UInt32				*dstMap,
    const UInt32		*mapA,
    const UInt32		*mapB,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPZeroBox(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPFillBox(
    UInt32				*dstMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPCopyBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPOrBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPXorBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPAndBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPBicBox(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

void BMAPInsetMap_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix );

int BMAPInsetMap_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels );

void BMAPTraceMap_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix );

int BMAPTraceMap_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels );

void BMAPInsetPatch_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

int BMAPInsetPatch_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels,
    const U16BoxPtr		box );

void BMAPTracePatch_1Pixel(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

int BMAPTracePatch_NPixels(
    UInt32				*dstMap0,
    UInt32				*dstMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				nPixels,
    const U16BoxPtr		box );

int BMAPGetSeedPatch(
    UInt32				*hSeed,
    UInt32				*vSeed,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box,
    UInt32				top );

UInt32 BMAPBucketTool(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed );

UInt32 BMAPBucketTool8Way(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed );

UInt32 BMAPBoxedBucketTool(
    U16BoxPtr			rgnBox,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed,
    const U16BoxPtr		box );

UInt32 BMAPBoxCGArea(
    BMAPBoxCGPtr		rgn,
    UInt32				*map,
    BMAPBucketPtr		stackBase,
    int					hPix,
    int					vPix,
    int					hSeed,
    int					vSeed );


#if defined(__cplusplus)
}
#endif
