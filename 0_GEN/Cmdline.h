

#pragma once


#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool IsArg( const char *pat, const char *argv );
bool GetArg( void *v, const char *pat, const char *argv );
bool GetArgStr( const char* &s, const char *pat, char *argv );
bool GetArgList( vector<int> &v, const char *pat, char *argv );
bool GetArgList( vector<double> &v, const char *pat, char *argv );


