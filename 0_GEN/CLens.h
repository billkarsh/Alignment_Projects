

#pragma once


#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CLens {

private:
	TForm	Tf[4],
			Ti[4];
	FILE	*flog;

public:
	bool ReadFile( const char *path, FILE* flog = stdout );
	bool ReadIDB( const string &idb, FILE* flog = stdout );

	int  CamID( const char *name );
	void UpdateDoublesRHS( double *D, const char *name, bool inv );
	void UpdateTFormRHS( TForm &T, const char *name, bool inv );
	void UpdateTFormLHS( TForm &T, const char *name, bool inv );

	const TForm& GetTf( const char *name ) {return Tf[CamID(name)];};
};


