

#pragma once


#include	<stdio.h>

#include	<map>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* DIR ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Maps name-tag to z-layer

class DIR {
private:
    map<string,int>	dirTbl;
public:
    void ReadDIRFile( const char *dirfile, FILE *FOUT );
    int  ZFromName( const char *name ) const;
};


