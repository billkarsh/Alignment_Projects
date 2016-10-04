/* BK_MEM.h ----------------------------------------------------------
 *
 * Memory management operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#ifndef	_CRT_SECURE_NO_DEPRECATE
#define	_CRT_SECURE_NO_DEPRECATE
#endif


#include	"BK_DEFS_GEN.h"

#include	<string.h>


#if defined(__cplusplus)
extern "C" {
#endif






/* --------------------------------------------------------------- */
/* Initialize ---------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	MEMZeroBytes( baseAddr, nBytes )						\
    memset( baseAddr, 0, nBytes )

#define	MEMZeroMap( map, hPix, vPix )							\
    MEMZeroBytes( map, BitToByte( (hPix) * (vPix) ) )

/* --------------------------------------------------------------- */
/* Copy ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

#if !TARGET_MAC
#define	MEMCopyBytes( dstBase, srcBase, nBytes )				\
    memcpy( dstBase, srcBase, nBytes )
#endif

#define	MEMCopyData16Bit( dstBase, srcBase, hPix, vPix )		\
    MEMCopyBytes( dstBase, srcBase,								\
    (vPix) * (hPix) * sizeof(UInt16) )

#define	MEMCopyData32Bit( dstBase, srcBase, hPix, vPix )		\
    MEMCopyBytes( dstBase, srcBase,								\
    (vPix) * (hPix) * sizeof(UInt32) )

#define	MEMCopyMap( dstMap, srcMap, hPix, vPix )				\
    MEMCopyBytes( dstMap, srcMap, BitToByte( (hPix) * (vPix) ) )

/* --------------------------------------------------------------- */
/* Comparisons --------------------------------------------------- */
/* --------------------------------------------------------------- */

SInt32 MEMFirstMatch32Bit(
    const UInt32	*list,
    SInt32			n,
    UInt32			value );

SInt32 MEMFirstNonMatch32Bit(
    const UInt32	*list,
    SInt32			n,
    UInt32			value );

int MEMEqual32Bit(
    const UInt32	*blk1,
    const UInt32	*blk2,
    UInt32			n );

/* --------------------------------------------------------------- */
/* Array operations ---------------------------------------------- */
/* --------------------------------------------------------------- */

void MEMInsert1Element(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			insBefore );

void MEMInsertNElements(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			insBefore,
    UInt32			nIns );

void MEMDelete1Element(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			delIdx );

void MEMDeleteNElements(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			delFirst,
    UInt32			nDel );

void MEMRoll1Element(
    void			*oneElemBuf,
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			oldPos,
    UInt32			newPos );

/* --------------------------------------------------------------- */
/* Flip byte order: x86 <-> VRAM format -------------------------- */
/* --------------------------------------------------------------- */

void MEMFlipEndianOrder16Bit(
    UInt16			*dst16Bit,
    UInt16			*src16Bit,
    UInt32			nWords );

void MEMFlipEndianOrder32Bit(
    UInt32			*dst32Bit,
    UInt32			*src32Bit,
    UInt32			nWords );


#if defined(__cplusplus)
}
#endif
