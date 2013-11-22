

#pragma once


#include	<semaphore.h>
#include	<stdio.h>


/* --------------------------------------------------------------- */
/* IMPORTANT ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/*	Build Notes
	-----------

	(1) Disk.h includes <semaphore.h>. There appears to be a
	conflict if <semaphore.h> follows "numerical_recipes.h".

	(2) To link against posix library, use 'g++ -lrt'.
*/

/* --------------------------------------------------------------- */
/* Semaphores ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CMutex {

private:
	sem_t	*mutex;
	char	m_name[256];

public:
	CMutex() : mutex(SEM_FAILED) {};

	bool Get( const char *name );
	void Release();
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool DskExists( const char *path );

double DskBytes( const char *path );

void DskCreateDir( const char *path, FILE* flog );

bool DskAbsPath( char *out, int bufsize, const char *in, FILE* flog );


