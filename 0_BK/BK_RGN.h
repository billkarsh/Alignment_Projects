/* BK_RGN.h ----------------------------------------------------------
 *
 * Operations on masked image regions.
 *
 * -----------
 * Definitions
 * -----------
 *
 * Intensity (I)	= power/unit area
 *
 * IOD				= "integrated object density"
 *
 * Flux (F)			= IOD
 *					= Integral of (I . dA)
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#include	"BK_SUM.h"


#if defined(__cplusplus)
extern "C" {
#endif


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

    UInt32	Ipk;	/* peak intensity	*/
    UInt16	v,
            h;

} RGN_Ipk_Rec, *RGN_Ipk_Ptr;

typedef struct {

    SUM_F_A_Rec		grains,
                    blob;

} RGN_Grain_Rec, *RGN_Grain_Ptr;

#if defined(WIN32)
#pragma pack( pop, PACK_S )
#elif defined(__linux__)
#pragma pack( pop )
#endif






/* --------------------------------------------------------------- */
/* Box measurements ---------------------------------------------- */
/* --------------------------------------------------------------- */

void RGNBox_Ipk(
    RGN_Ipk_Ptr			ipk,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

void RGNBox_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

void RGNBox_F_F2_A(
    SUM_F_F2_A_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

UInt32 RGNBoxPixList(
    UInt32				*list,
    const UInt16		*data16Bit,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		box );

UInt32 RGNHistogramBox(
    UInt32				*oFlowCnt,
    UInt32				*oFlowSum,
    UInt32				*binArray,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

/* --------------------------------------------------------------- */
/* Box perimeter measurements ------------------------------------ */
/* --------------------------------------------------------------- */

void RGNBoxPerim_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

UInt32 RGNBoxPerimList(
    UInt32				*list,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    const U16BoxPtr		box );

/* --------------------------------------------------------------- */
/* Blob measurements --------------------------------------------- */
/* --------------------------------------------------------------- */

UInt32 RGNBlob_A(
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

UInt32 RGNBlob_F(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

FP32 RGNBlob_I(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

UInt32 RGNBlob_Ipk(
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

UInt32 RGNHistogramBlob(
    UInt32				*oFlowCnt,
    UInt32				*oFlowSum,
    UInt32				*binArray,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

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
    UInt32				*visMap );

/* --------------------------------------------------------------- */
/* Blob operations ----------------------------------------------- */
/* --------------------------------------------------------------- */

void RGNRethresholdBlob(
    const UInt16		*data16Bit,
    UInt32				*map,
    UInt32				hPix,
    UInt32				vPix,
    UInt32				thresh,
    const U16BoxPtr		rBlob );

UInt32 *RGNErodeBlob(
    UInt32				*erMap0,
    UInt32				*erMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    UInt32				N );

UInt32 *RGNDilateBlob(
    UInt32				*dilMap0,
    UInt32				*dilMap1,
    const UInt32		*srcMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob,
    UInt32				N );

void RGNGetHolesBlob(
    UInt32				*dstMap,
    const UInt32		*srcMap,
    UInt32				*tmpMap,
    void				*bucketScratch,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

void RGNDelTouchingRegion(
    UInt32				*dstMap,
    const U16BoxPtr		rDst,
    const UInt32		*rgnMap,
    const U16BoxPtr		rRgn,
    UInt32				*tmpMap0,
    UInt32				*tmpMap1,
    void				*bucketScratch,
    UInt32				hPix,
    UInt32				vPix );


#if defined(__cplusplus)
}
#endif
