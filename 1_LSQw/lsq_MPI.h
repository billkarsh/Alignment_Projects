

#pragma once


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern int	wkid, nwks;	// my workid, count

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MPIInit( int& argc, char**& argv );

void MPIExit();

void MPIWaitForOthers();

bool MPISend( void* buf, int bytes, int wdst, int tag );
bool MPIRecv( void* buf, int bytes, int wsrc, int tag );


