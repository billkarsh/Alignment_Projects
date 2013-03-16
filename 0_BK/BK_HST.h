/* BK_HST.h ----------------------------------------------------------
 *
 * Histogramming.
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






UInt32 HSTAltImageRowsUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	int				hPix,
	int				vPix );

UInt32 HSTImageBoxUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	int				hPix,
	int				vPix,
	const U16BoxPtr	box );

UInt32 HSTUnitWidthUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	UInt32			nData );

int HSTGeneralFP32(
	FP32			*uFlow,
	FP32			*oFlow,
	UInt16			*binArray,
	UInt32			nBins,
	const FP32		*pFirstVal,
	UInt32			recBytes,
	UInt32			nData,
	FP32			minH,
	FP32			maxH );


#if defined(__cplusplus)
}
#endif


