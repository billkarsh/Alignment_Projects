/* BK_SORT.h ---------------------------------------------------------
 *
 * Heap sort.
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
/* Callbacks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef int (*SORTSInt32Proc)( SInt32 A, SInt32 B );

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void SORTSInt32(
	SInt32				*base,
	SInt32				nItems,
	SORTSInt32Proc		compare );

void SORTImmedAscending( SInt32 *base, UInt32 nItems );

void SORTImmedDescending( SInt32 *base, UInt32 nItems );


#if defined(__cplusplus)
}
#endif
