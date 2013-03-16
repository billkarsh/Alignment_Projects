/* BK_SORT_FP32.h ----------------------------------------------------
 *
 * Heap sort for FP32 arrays.
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
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void SORT_FP32_Ascending( FP32 *base, UInt32 nItems );

void SORT_FP32_Descending( FP32 *base, UInt32 nItems );


#if defined(__cplusplus)
}
#endif
