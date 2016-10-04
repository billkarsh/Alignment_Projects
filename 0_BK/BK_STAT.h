/* BK_STAT.h ---------------------------------------------------------
 *
 * Common statistics calculations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_DEFS_GEN.h"


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

    FP64	icpt,
            slope,
            linCor;

} STATLinFitRec, *STATLinFitPtr;

#if defined(WIN32)
#pragma pack( pop, PACK_S )
#elif defined(__linux__)
#pragma pack( pop )
#endif

/* --------------------------------------------------------------- */
/* Integers ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int STATAverageUInt16( const UInt16 *array, UInt32 N );

UInt32 STATIndexOfMaxValueUInt16( const UInt16 *array, UInt32 N );

UInt32 STATIndexOfMaxValueUInt32( const UInt32 *array, UInt32 N );

UInt32 STATPercentileBinUInt32(
    const UInt32	*hist,
    UInt32			nBins,
    UInt32			nCounts,
    FP64			frac );

/* --------------------------------------------------------------- */
/* Floating-point ------------------------------------------------ */
/* --------------------------------------------------------------- */

FP32 STATSigma32( FP32 N, FP32 meanX, FP32 sumX2 );

FP32 STATAveFP32(
    const FP32		*pFirstVal,
    UInt32			recBytes,
    UInt32			N );

FP32 STATMinMaxFP32(
    FP32			*minVal,
    const FP32		*pFirstVal,
    UInt32			recBytes,
    UInt32			N );

FP32 STATAveAndSDInRangeFP32(
    FP32			*stdDev,
    const FP32		*array,
    UInt32			i0,
    UInt32			iLim );

FP32 STATWindowedMaxFP32(
    const FP32		*array,
    UInt32			i0,
    UInt32			iLim,
    UInt32			windowSize );

FP32 STATMedianFP32( FP32 *array, FP32 *wrkspc, UInt32 N );

void STATUnitWeightsFP32( FP32 *wt, UInt32 N );

FP32 STATSigmaMADFP32( FP32 *r, FP32 *wrkspc, UInt32 N );

int STATTukeyWtsFP32( FP32 *w, FP32 *r, FP32 *wrkspc, UInt32 N );

int STATBisquareWtsFP32( FP32 *w, FP32 *r, FP32 *wrkspc, UInt32 N );

void STATSlopeFit32XY(
    STATLinFitPtr	L,
    const FP32		*x,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim );

void STATLinFit32XY(
    STATLinFitPtr	L,
    const FP32		*x,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim );

void STATLinFit64XY(
    STATLinFitPtr	L,
    const FP64		*x,
    const FP64		*y,
    UInt32			i0,
    UInt32			iLim );

void STATLinFit32Y(
    STATLinFitPtr	L,
    FP32			x0,
    FP32			dx,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim );

#if defined(__cplusplus)
}
#endif
