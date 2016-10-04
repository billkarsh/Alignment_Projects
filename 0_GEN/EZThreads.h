

#pragma once


#include	<pthread.h>
#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// On entry, gets ithr = thread index [0..nthr).
// On exit, you should return NULL.
//
typedef	void* (*EZThreadproc)( void* ithr );

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool EZThreads(
    EZThreadproc	proc,
    int				nthr,
    int				stksize_factor,
    const char		*msgname,
    FILE			*flog = stdout );


