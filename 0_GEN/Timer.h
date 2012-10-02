

#pragma once


#include	<stdio.h>
#include	<sys/times.h>


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void    Yield_usec( int microsec );

clock_t StartTiming();
double  DeltaSeconds( clock_t start );
clock_t	StopTiming( FILE *flog, const char *msg, clock_t start );


