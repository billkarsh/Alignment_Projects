

#include	"Cmdline.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* IsArg --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if matched command line parameter.
//
// Example usage:
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( IsArg( "-nf", argv[i] ) )
//			NoFolds = true;
//	}
//
bool IsArg( const char *pat, const char *argv )
{
    return !strcmp( argv, pat );
}

/* --------------------------------------------------------------- */
/* GetArg -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read argument from command line.
//
// Example usage:
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArg( &ApproxScale, "-SCALE=%lf", argv[i] ) )
//			;
//		else if( GetArg( &Order, "-ORDER=%d", argv[i] ) )
//			;
//	}
//
bool GetArg( void *v, const char *pat, const char *argv )
{
    return 1 == sscanf( argv, pat, v );
}

/* --------------------------------------------------------------- */
/* GetArgStr -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Point at string argument on command line.
//
// Example usage:
//
//	char	*dirptr;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgStr( dirptr, "-d=", argv[i] ) )
//			;
//	}
//
bool GetArgStr( const char* &s, const char *pat, char *argv )
{
    int	len = strlen( pat );

    if( !strncmp( argv, pat, len ) ) {

        s = argv + len;
        return true;
    }

    return false;
}

/* --------------------------------------------------------------- */
/* GetArgList ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read integer argument list from command line.
//
// Example usage: ... -List=2,5,7 ...
//
//	vector<int>	I;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgList( I, "-List=", argv[i] ) )
//			;
//	}
//
bool GetArgList( vector<int> &v, const char *pat, char *argv )
{
    int	len = strlen( pat );

    if( !strncmp( argv, pat, len ) ) {

        char	*s = strtok( argv + len, ":;, " );

        v.clear();

        while( s ) {
            v.push_back( atoi( s ) );
            s = strtok( NULL, ":;, " );
        }

        return true;
    }

    return false;
}

// Read double argument list from command line.
//
// Example usage: ... -List=2.7,5,1.8e7 ...
//
//	vector<double>	D;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgList( D, "-List=", argv[i] ) )
//			;
//	}
//
bool GetArgList( vector<double> &v, const char *pat, char *argv )
{
    int	len = strlen( pat );

    if( !strncmp( argv, pat, len ) ) {

        char	*s = strtok( argv + len, ":;, " );

        v.clear();

        while( s ) {
            v.push_back( atof( s ) );
            s = strtok( NULL, ":;, " );
        }

        return true;
    }

    return false;
}


