/* BK_BMAP_SET.h -----------------------------------------------------
 *
 * Bitmap Set (draw) operations.
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






void BMAPSetPoint(
	UInt32				*map,
	UInt32				hPix,
	int					h,
	int					v );

void BMAPSetHSeg(
	UInt32				*map,
	int					hPix,
	int					h1,
	int					h2,
	int					v );

void BMAPSetVSeg(
	UInt32				*map,
	int					hPix,
	int					h,
	int					v1,
	int					v2 );

void BMAPSetSeg(
	UInt32				*map,
	UInt32				hPix,
	int					h1,
	int					v1,
	int					h2,
	int					v2 );

void BMAPSetBox(
	UInt32				*map,
	int					hPix,
	int					vPix,
	const U16BoxPtr		box );


#if defined(__cplusplus)
}
#endif
