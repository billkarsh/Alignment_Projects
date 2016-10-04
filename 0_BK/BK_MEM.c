/* BK_MEM ------------------------------------------------------------
 *
 * Memory management operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_MEM.h"






/* MEMFirstMatch32Bit ------------------------------------------------
 *
 * Return index of first occurrence of value in list,
 * or, return -1 if not present.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

SInt32 MEMFirstMatch32Bit(
    const UInt32	*list,
    SInt32			n,
    UInt32			value )
{
    SInt32	i;

    for( i = 0; i < n; ++i ) {

        if( *list++ == value )
            return i;
    }

    return -1;
}


/* MEMFirstNonMatch32Bit ---------------------------------------------
 *
 * Return index of first element that differs from value,
 * or, return -1 if all match.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

SInt32 MEMFirstNonMatch32Bit(
    const UInt32	*list,
    SInt32			n,
    UInt32			value )
{
    SInt32	i;

    for( i = 0; i < n; ++i ) {

        if( *list++ != value )
            return i;
    }

    return -1;
}


/* MEMEqual32Bit -----------------------------------------------------
 *
 * Return true if two UInt32-blocks are equal.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int MEMEqual32Bit(
    const UInt32	*blk1,
    const UInt32	*blk2,
    UInt32			n )
{
    UInt32	i;

    for( i = 0; i < n; ++i ) {

        if( *blk1++ != *blk2++ )
            return false;
    }

    return true;
}


/* MEMInsert1Element -------------------------------------------------
 *
 * Insert space for one new element into an array of
 * similar elements.
 *
 * insBefore
 * ---------
 * Current index of element that will immediately follow
 * the inserted item.
 *
 * Nothing is done if appending (insBefore >= nElems).
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMInsert1Element(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			insBefore )
{
    char	*src;

    if( insBefore < nElems ) {

        src = (char*)baseAddr + insBefore * elemBytes;

        memmove( src + elemBytes, src,
            (nElems - insBefore) * elemBytes );
    }
}


/* MEMInsertNElements ------------------------------------------------
 *
 * Insert space for nIns new elements into an array of
 * similar elements.
 *
 * insBefore
 * ---------
 * Current index of element that will immediately follow
 * the inserted items.
 *
 * Nothing is done if appending (insBefore >= nElems).
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMInsertNElements(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			insBefore,
    UInt32			nIns )
{
    char	*src;

    if( nIns && insBefore < nElems ) {

        src = (char*)baseAddr + insBefore * elemBytes;

        memmove( src + nIns * elemBytes, src,
            (nElems - insBefore) * elemBytes );
    }
}


/* MEMDelete1Element -------------------------------------------------
 *
 * - Delete one of an array of similar elements.
 *
 * - Recompact the array to fill the gap.
 *
 * Nothing is done if the last element is specified.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMDelete1Element(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			delIdx )
{
    char	*dst;

    if( nElems && delIdx < --nElems ) {

        dst = (char*)baseAddr + delIdx * elemBytes;

        MEMCopyBytes( dst, dst + elemBytes,
            (nElems - delIdx) * elemBytes );
    }
}


/* MEMDeleteNElements ------------------------------------------------
 *
 * Given an array of similar elements:
 *
 * - Delete a run of nDel elements, beginning with delFirst.
 *
 * - Recompact the array to fill the gap.
 *
 * Nothing is done if the last elements are specified.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMDeleteNElements(
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			delFirst,
    UInt32			nDel )
{
    char	*dst;

    if( nElems && nDel && delFirst + nDel < nElems ) {

        dst = (char*)baseAddr + delFirst * elemBytes;

        MEMCopyBytes( dst, dst + nDel * elemBytes,
            (nElems - delFirst - nDel) * elemBytes );
    }
}


/* MEMRoll1Element ---------------------------------------------------
 *
 * Given an array of similar elements:
 *
 * Roll elements in range [oldPos, newPos] by 1 position.
 *
 * If oldPos > newPos, the element at the far end of the range
 * (oldPos) is rolled forward and everybody else moves back.
 *
 * If oldPos < newPos, the near element (oldPos) is moved back
 * and everybody else comes forward.
 *
 * If oldPos = newPos nothing is done.
 *
 * oneElemBuf
 * ----------
 * Caller-provided workspace able to hold one element.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMRoll1Element(
    void			*oneElemBuf,
    void			*baseAddr,
    UInt32			elemBytes,
    UInt32			nElems,
    UInt32			oldPos,
    UInt32			newPos )
{
/* ------------- */
/* Sanity checks */
/* ------------- */

    if( nElems <= 1 )
        goto exit;

    if( oldPos >= nElems )
        oldPos = nElems - 1;

    if( newPos >= nElems )
        newPos = nElems - 1;

    if( oldPos == newPos )
        goto exit;

/* -------- */
/* Make gap */
/* -------- */

    MEMCopyBytes( oneElemBuf,
        (char*)baseAddr + oldPos * elemBytes, elemBytes );

/* -------------------- */
/* Slide other elements */
/* -------------------- */

    if( oldPos > newPos ) {

        char	*src = (char*)baseAddr + newPos * elemBytes;

        memmove( src + elemBytes, src,
            (oldPos - newPos) * elemBytes );
    }
    else {

        char	*dst = (char*)baseAddr + oldPos * elemBytes;

        MEMCopyBytes( dst, dst + elemBytes,
            (newPos - oldPos) * elemBytes );
    }

/* -------- */
/* Fill gap */
/* -------- */

    MEMCopyBytes( (char*)baseAddr + newPos * elemBytes,
        oneElemBuf, elemBytes );

exit:
    return;
}


/* MEMFlipEndianOrder16Bit -------------------------------------------
 *
 * Change 16-bit numeric array from big- to little-endian order
 * or vice versa.
 *
 * -----
 * Notes
 * -----
 *
 * dst16Bit and src16Bit can be same array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMFlipEndianOrder16Bit(
    UInt16			*dst16Bit,
    UInt16			*src16Bit,
    UInt32			nWords )
{
    UInt32	i;
    UInt16	D;

    for( i = 0; i < nWords; ++i ) {

        D = *src16Bit++;

        *dst16Bit = (D << ByteBits) + (D >> ByteBits);
    }
}


/* MEMFlipEndianOrder32Bit -------------------------------------------
 *
 * Change 32-bit numeric array from big- to little-endian order
 * or vice versa.
 *
 * -----
 * Notes
 * -----
 *
 * dst32Bit and src32Bit can be same array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void MEMFlipEndianOrder32Bit(
    UInt32			*dst32Bit,
    UInt32			*src32Bit,
    UInt32			nWords )
{
    UInt32	i, D;

    for( i = 0; i < nWords; ++i ) {

        D = *src32Bit++;

        *dst32Bit++ =
            (D << (WordBits - ByteBits)) +
            ((D << ByteBits) & 0x00FF0000) +
            ((D >> ByteBits) & 0x0000FF00) +
            (D >> (WordBits - ByteBits));
    }
}


