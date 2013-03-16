/* BK_SORT_FP32 ------------------------------------------------------
 *
 * Heap sort for FP32 arrays.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_SORT_FP32.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	Slot4( k )										\
	(*(FP32*)((char*)base + (k)))

#define	oneX4											\
	sizeof(FP32)






/* SORT_FP32_Ascending -----------------------------------------------
 *
 * Heapsort zero-based array of 4-byte FP32 values
 * into ascending order.
 *
 * The values are moved about in their array.
 *
 * - base
 * Base address of array.
 *
 * - nItems
 * Number of values in array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SORT_FP32_Ascending( FP32 *base, UInt32 nItems )
{
	FP32	temp;
	UInt32	iX4, jX4, kX4, nX4;

	if( nItems < 2 )
		goto exit;

	kX4 = ((nItems >> 1) + 1) * sizeof(FP32);
	nX4 = nItems * sizeof(FP32);
	--base;			/* -1 for zero-offset array */

	for(;;) {

		if( kX4 > oneX4 )
			temp = Slot4( kX4 -= oneX4 );
		else {

			temp = Slot4( nX4 );
			Slot4( nX4 ) = Slot4( oneX4 );

			if( (nX4 -= oneX4) == oneX4 ) {
				Slot4( oneX4 ) = temp;
				break;
			}
		}

		iX4 = kX4;
		jX4 = kX4 + kX4;

		while( jX4 <= nX4 ) {

			if( jX4 < nX4 &&
				Slot4( jX4 ) < Slot4( jX4 + oneX4 ) ) {

				jX4 += oneX4;
			}

			if( temp < Slot4( jX4 ) ) {

				Slot4( iX4 ) = Slot4( jX4 );
				iX4  = jX4;
				jX4 += jX4;
			}
			else
				break;
		}

		Slot4( iX4 ) = temp;
	}

exit:
	return;
}


/* SORT_FP32_Descending ----------------------------------------------
 *
 * Heapsort zero-based array of 4-byte FP32 values
 * into descending order.
 *
 * The values are moved about in their array.
 *
 * - base
 * Base address of array.
 *
 * - nItems
 * Number of values in array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void SORT_FP32_Descending( FP32 *base, UInt32 nItems )
{
	FP32	temp;
	UInt32	iX4, jX4, kX4, nX4;

	if( nItems < 2 )
		goto exit;

	kX4 = ((nItems >> 1) + 1) * sizeof(FP32);
	nX4 = nItems * sizeof(FP32);
	--base;			/* -1 for zero-offset array */

	for(;;) {

		if( kX4 > oneX4 )
			temp = Slot4( kX4 -= oneX4 );
		else {

			temp = Slot4( nX4 );
			Slot4( nX4 ) = Slot4( oneX4 );

			if( (nX4 -= oneX4) == oneX4 ) {
				Slot4( oneX4 ) = temp;
				break;
			}
		}

		iX4 = kX4;
		jX4 = kX4 + kX4;

		while( jX4 <= nX4 ) {

			if( jX4 < nX4 &&
				Slot4( jX4 + oneX4 ) < Slot4( jX4 ) ) {

				jX4 += oneX4;
			}

			if( Slot4( jX4 ) < temp ) {

				Slot4( iX4 ) = Slot4( jX4 );
				iX4  = jX4;
				jX4 += jX4;
			}
			else
				break;
		}

		Slot4( iX4 ) = temp;
	}

exit:
	return;
}


