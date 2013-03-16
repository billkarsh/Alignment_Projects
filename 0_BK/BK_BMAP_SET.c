/* BK_BMAP_SET -------------------------------------------------------
 *
 * Bitmap Set (draw) operations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_BMAP_SET.h"

#define	_USE_MATH_DEFINES
#include	<math.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	HToMask( h )											\
	(HiBit >> WordRem( h ))






/* BMAPSetPoint ------------------------------------------------------
 *
 * Set (h,v) to one.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPSetPoint(
	UInt32				*map,
	UInt32				hPix,
	int					h,
	int					v )
{
	map += v * BitToWord( hPix );

	map[BitToWord( h )] |= HToMask( h );
}


/* BMAPSetHSeg -------------------------------------------------------
 *
 * Set horizontal segment (h1,v) - (h2,v) to ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPSetHSeg(
	UInt32				*map,
	int					hPix,
	int					h1,
	int					h2,
	int					v )
{
	int		h;

	if( h1 > h2 ) {
		h	= h1;
		h1	= h2;
		h2	= h;
	}

	map += v * BitToWord( hPix );

	for( h = h1; h <= h2; ++h )
		map[BitToWord( h )] |= HToMask( h );
}


/* BMAPSetVSeg -------------------------------------------------------
 *
 * Set vertical segment (h,v1) - (h,v2) to ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPSetVSeg(
	UInt32				*map,
	int					hPix,
	int					h,
	int					v1,
	int					v2 )
{
	UInt32	rowBytes	= BitToByte( hPix ),
			mask		= HToMask( h );
	int		v;

	if( v1 > v2 ) {
		v	= v1;
		v1	= v2;
		v2	= v;
	}

	map += v1 * BitToWord( hPix ) + BitToWord( h );

	for( v = v1; v <= v2;
		++v,
		map = (UInt32*)((char*)map + rowBytes) ) {

		*map |= mask;
	}
}


/* BMAPSetSeg --------------------------------------------------------
 *
 * Set any segment (h1,v1) - (h2,v2) to ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPSetSeg(
	UInt32				*map,
	UInt32				hPix,
	int					h1,
	int					v1,
	int					h2,
	int					v2 )
{
	double	r, dr, T, c, s;
	int		dh, dv, N, i;

	dh	= h2 - h1;
	dv	= v2 - v1;
	T	= atan( (double)dv / dh ) + (dh >= 0 ? 0.0 : M_PI);
	r	= sqrt( dh*dh + dv*dv );
	c	= cos( T );
	N	= (int)(r * 1.05);	/* overfill 5% */
	s	= sin( T );
	dr	= r / N;

	for( i = 0; i < N; ++i ) {

		BMAPSetPoint( map, hPix,
			h1 + (int)(i*dr*c), v1 + (int)(i*dr*s) );
	}
}


/* BMAPSetBox --------------------------------------------------------
 *
 * Set box perimeter to ones.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void BMAPSetBox(
	UInt32				*map,
	int					hPix,
	int					vPix,
	const U16BoxPtr		box )
{
	int		bot, right;

	if( (bot = box->bottom) > vPix )
		bot = vPix;

	if( (right = box->right) > hPix )
		right = hPix;

	if( box->top >= bot )
		goto exit;

	if( box->left >= right )
		goto exit;

	--bot;
	--right;

	BMAPSetHSeg( map, hPix, box->left, right, box->top );
	BMAPSetVSeg( map, hPix, box->left, box->top, bot );
	BMAPSetVSeg( map, hPix, right, box->top, bot );
	BMAPSetHSeg( map, hPix, box->left, right, bot );

exit:
	return;
}


