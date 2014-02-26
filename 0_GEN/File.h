

#pragma once


#include	<stdio.h>
#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* class CLineScan ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Wrapper around stdio::getline:
// - Helps document proper usage.
// - Helps manage char buffer memory.
//
// Example with locally declared CLineScan:
//
//	{ // open scope
//		CLineScan	LS;
//		for(;;) {
//			if( LS.Get( myfileptr ) <= 0 )
//				break;
//			// read LS.line
//		}
//	} // close scope - auto calls destructor
//
// Example with dynamically allocated CLineScan:
//
//	CLineScan	*ls = new CLineScan;
//	for(;;) {
//		if( ls->Get( myfileptr ) <= 0 )
//			break;
//		// read ls->line
//	}
//	delete ls;
//
class CLineScan {

public:
	char	*line;		// pointer to buffer with read chars
	size_t	bufsize;	// cur physical size of line buffer

public:
	CLineScan() : line(NULL), bufsize(0) {};
	virtual ~CLineScan()	{Release();};

	void Release()
		{
			if( line ) {
				free( line );
				line = NULL;
				bufsize = 0;
			}
		};

	// Return num chars read including newline
	// ( <= 0 implies exhausted lines or error).
	//
	// Note getline will realloc line buffer and
	// update bufsize if input size requires it.
	//
	int Get( FILE *f )
		{return getline( &line, &bufsize, f );};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

FILE *FileOpenOrDie(
	const char	*name,
	const char	*rw,
	FILE		*flog = NULL );

void FileScriptPerms( const char *path );

const char* FileNamePtr( const char *path );
const char* FileDotPtr( const char *path );
char* FileCloneNamePart( const char *path );

bool FileIsExt( const char *path, const char *ext );


