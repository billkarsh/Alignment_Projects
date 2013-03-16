/* BK_BMAP_CONVERT.h -------------------------------------------------
 *
 * Bitmap format conversion operations.
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






void BMAPConvertDepth8To1(
	UInt32				*dstMap,
	const UInt8			*srcMap,
	UInt32				hPix,
	UInt32				vPix );

void BMAPConvertDepth1To8(
	UInt8				*dstMap,
	const UInt32		*srcMap,
	UInt32				hPix,
	UInt32				vPix,
	UInt8				oneAs8bit );


#if defined(__cplusplus)
}
#endif
