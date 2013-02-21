

#pragma once


#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CAffineLens {

private:
	TAffine	Tf[4],
			Ti[4];
	FILE	*flog;

public:
	bool ReadFile( const char *path, FILE* flog = stdout );
	bool ReadIDB( const string &idb, FILE* flog = stdout );

	int  CamID( const char *name );
	void UpdateDoublesRHS( double *D, const char *name, bool inv );
	void UpdateTFormRHS( TAffine &T, const char *name, bool inv );
	void UpdateTFormLHS( TAffine &T, const char *name, bool inv );

	const TAffine& GetTf( const char *name ) {return Tf[CamID(name)];};
};


