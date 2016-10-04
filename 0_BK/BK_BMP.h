/* BK_BMP.h ----------------------------------------------------------
 *
 * Bitmap / Image operations.
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
/* Tests --------------------------------------------------------- */
/* --------------------------------------------------------------- */

int BMPMapValue( const UInt32 *map, UInt32 hPix, int h, int v );

#define	BMPPixelValue( pdata16Bit, hPix, h, v )					\
    ((pdata16Bit)[(v) * (hPix) + h])

/* --------------------------------------------------------------- */
/* Bitmap generation --------------------------------------------- */
/* --------------------------------------------------------------- */

UInt32 BMPGetMode(
    const UInt16	*data16Bit,
    UInt32			hPix,
    UInt32			vPix );

UInt32 BMPGetModeBox(
    const UInt16	*data16Bit,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box );

void BMPGenerateThreshMap16Bit(
    const UInt16	*data16Bit,
    UInt32			*map,
    UInt32			hPix,
    UInt32			vPix,
    UInt32			thresh );

void BMPGenerateThreshPatch16Bit(
    const UInt16	*data16Bit,
    UInt32			*map,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box,
    UInt32			thresh );

void BMPFillHolesBox(
    UInt32			*dstMap,
    UInt32			*trackingMap,
    UInt32			*scanningMap,
    UInt32			*tmpMap0,
    UInt32			*tmpMap1,
    void			*bucketScratch,
    UInt32			hPix,
    UInt32			vPix,
    const U16BoxPtr	box );


#if defined(__cplusplus)
}
#endif
