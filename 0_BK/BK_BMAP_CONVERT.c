/* BK_BMAP_CONVERT ---------------------------------------------------
 *
 * Bitmap format conversion operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_BMAP_CONVERT.h"






/* BMAPConvertDepth8To1 ----------------------------------------------
 *
 * Convert a map of depth 8 bits/pixel to equivalent
 * standard bitmap of depth 1 bit/pixel.
 *
 * -----
 * Notes
 * -----
 *
 * An 8 bit deep map that represents a "bitmap" is a contiguous
 * set of pixels, each of which has value zero, or, either {1 or
 * 0xFF}.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPConvertDepth8To1(
	UInt32				*dstMap,
	const UInt8			*srcMap,
	UInt32				hPix,
	UInt32				vPix )
{
	UInt32	word, nWords, W;
	int		b;

	nWords = BitToWord( vPix * hPix );

	for( word = 0; word < nWords; ++word ) {

		for( b = 0; b < WordBits; ++b, W <<= 1 )
			W |= *srcMap++ & 1;

		dstMap[word] = W;
	}
}


/* BMAPConvertDepth1To8 ----------------------------------------------
 *
 * Convert a standard bitmap of depth 1 bit/pixel to equivalent
 * map of depth 8 bits/pixel.
 *
 * oneAs8bit
 * ---------
 * An 8 bit deep map that represents a "bitmap" is a contiguous
 * set of pixels, each of which has value zero or "one." The
 * caller specifies the value to be used as one, which is usually
 * either 1 or 0xFF.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPConvertDepth1To8(
	UInt8				*dstMap,
	const UInt32		*srcMap,
	UInt32				hPix,
	UInt32				vPix,
	UInt8				oneAs8bit )
{
	UInt32	word, nWords, W;
	int		b;

	nWords	= BitToWord( vPix * hPix );

	for( word = 0; word < nWords; ++word ) {

		W = *srcMap++;

		for( b = 0; b < WordBits; ++b, W <<= 1, ++dstMap )
			*dstMap	= (UInt8)(oneAs8bit & SignExtHiBit32( W ));
	}
}


