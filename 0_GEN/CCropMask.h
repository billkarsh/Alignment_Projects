

#pragma once


#include	"GenDefs.h"

#include	<stdio.h>

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CCropMask {

private:
    IBox	B[4];
    FILE	*flog;
    int		isfile;	// unset if -1

public:
    CCropMask() : isfile(-1) {};

    bool IsFile( const string &idb );

    bool ReadIDB( const string &idb, FILE* flog = stdout );

    bool GetBox( IBox &b, int cam )
        {b=B[cam]; return isfile;};
};


