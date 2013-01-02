

#pragma once


#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CLens {

public:
	TForm	Tf[4],
			Ti[4];
	FILE	*flog;

public:
	bool ReadFile( const char *path, FILE* flog = stdout );
	bool ReadIDB( const string &idb, FILE* flog = stdout );
	void Tdfm( TForm &T, int a, int b );
	int  CamID( const char *name );

	void PrintArg(
		char		*buf,
		const char	*aname,
		const char	*bname );
};


