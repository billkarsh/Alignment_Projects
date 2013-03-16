/* BK_GEOM.h ---------------------------------------------------------
 *
 * Common geometry calculations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#include	"BK_DEFS_GEN.h"

#if defined(WIN32)
#include	<windows.h>
#endif


#if defined(__cplusplus)
extern "C" {
#endif






void GEOMInsetU16Box1Pixel( U16BoxPtr dstBox, U16BoxPtr srcBox );

void GEOMOutsetU16Box(
	U16BoxPtr	dstBox,
	U16BoxPtr	srcBox,
	int			hLim,
	int			vLim,
	int			outset );

void GEOMU16PatchBox(
	U16BoxPtr	dstBox,
	U16BoxPtr	srcBox,
	int			hLim,
	int			vLim );

int GEOMU16BoxIntersection(
	U16BoxPtr		outBox,
	const U16BoxPtr	box1,
	const U16BoxPtr	box2 );

void GEOMU16BoxUnion(
	U16BoxPtr		outBox,
	const U16BoxPtr	box1,
	const U16BoxPtr	box2 );

#if defined(WIN32)
int GEOMPinSubrect( RECT *r, const RECT *R );
#endif

FP32 GEOMSegLen( int dx, int dy );


#if defined(__cplusplus)
}
#endif
