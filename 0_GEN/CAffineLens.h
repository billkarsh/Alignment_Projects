

#pragma once


#include	"TAffine.h"

#include	<string>
using namespace std;


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

    void UpdateDoublesRHS( double *D, int cam, bool inv );
    void UpdateTFormRHS( TAffine &T, int cam, bool inv );
    void UpdateTFormLHS( TAffine &T, int cam, bool inv );

    const TAffine& GetTf( int cam ) {return Tf[cam];};
};


