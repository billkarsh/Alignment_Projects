/* BK_DEFS_GEN.h -----------------------------------------------------
 *
 * Definitions common to all code.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#if defined(__cplusplus)
extern "C" {
#endif


/* ---------------------------------------------------------------*/
/* Basic types ---------------------------------------------------*/
/* ---------------------------------------------------------------*/

typedef	unsigned int	UInt32;
typedef unsigned short	UInt16;
typedef unsigned char	UInt8;
typedef	int				SInt32;
typedef	short			SInt16;
typedef char			SInt8;

typedef double			FP64;
typedef float			FP32;

/* ---------------------------------------------------------------*/
/* Constants -----------------------------------------------------*/
/* ---------------------------------------------------------------*/

#ifndef	false
#define	false	0
#endif
#ifndef	true
#define	true	1
#endif


#ifndef	NULL
#ifdef	__cplusplus
#define	NULL	0
#else
#define	NULL	((void *)0)
#endif
#endif

/* ---------------------------------------------------------------*/
/* Bits and Pieces -----------------------------------------------*/
/* ---------------------------------------------------------------*/

enum DEFS_GENConsts {
    ByteBits			= 8,
    HalfBits			= 16,
    WordBits			= 32,
    AllButOneBit		= WordBits - 1
};


#define	HiMask			0xFFFF0000L
#define	LoMask			0x0000FFFFL
#define	HiBit			0x80000000L
#define	Frac_One		0x00010000L
#define	Frac_OneHalf	0x00008000L


#define	BitToByte( bit )										\
    ((bit) >> 3)

#define	ByteRem( bit )											\
    ((bit) & (ByteBits - 1))

#define	BitToWord( bit )										\
    ((bit) >> 5)

#define	WordRem( bit )											\
    ((bit) & AllButOneBit)

#define	ByteToWord( byte )										\
    ((byte) >> 2)

#define	SignExtHiBit32( W )										\
    ((long)(W) >> AllButOneBit)

#define	IsHiBitSet32( W )										\
    (((W) & HiBit) != 0)

/* ---------------------------------------------------------------*/
/* Common macros -------------------------------------------------*/
/* ---------------------------------------------------------------*/

#define	MAX( a, b )	((a) >= (b) ? (a) : (b))
#define	MIN( a, b )	((a) <= (b) ? (a) : (b))

#define	PHYSIZE( txtbuf )	sizeof(txtbuf)
#define	LOGSIZE( txtbuf )	(sizeof(txtbuf) - 1)

/* ---------------------------------------------------------------*/
/* Common data structures ----------------------------------------*/
/* ---------------------------------------------------------------*/

#if defined(WIN32)
#pragma pack( push, PACK_S )
#pragma pack( 1 )
#elif defined(__linux__)
#pragma pack( push, 1 )
#pragma pack( 1 )
#endif

typedef struct {

    UInt16	top,
            left,
            bottom,
            right;

} U16Box, *U16BoxPtr;

typedef struct {

    UInt32	top,
            left,
            bottom,
            right;

} U32Box, *U32BoxPtr;

typedef struct {	/* same as Macintosh Rect */

    SInt16	top,
            left,
            bottom,
            right;

} S16Box, *S16BoxPtr;

typedef struct {	/* same as PC RECT */

    SInt32	top,
            left,
            bottom,
            right;

} S32Box, *S32BoxPtr;

typedef struct {

    UInt16	x,
            y;

} Centroid, *CentroidPtr;

typedef struct {

    FP32	x,
            y;

} F32Pnt, *F32PntPtr;

typedef struct {

    FP32	top,
            left,
            bottom,
            right;

} F32Box, *F32BoxPtr;

typedef struct {

    UInt16	y0,
            x0,
            dy,
            dx;

} U16ROI, *U16ROIPtr;

#if defined(WIN32)
#pragma pack( pop, PACK_S )
#elif defined(__linux__)
#pragma pack( pop )
#endif


#if defined(__cplusplus)
}
#endif
