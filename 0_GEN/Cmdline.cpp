

#include	"Cmdline.h"

#include	<string>

#ifdef USE_CURL
#include	<curl/curl.h>
#endif

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

		char	*s = argv + len;

		if( strstr( s, "http:" ) )
			return false;

		s = strtok( s, ":;, " );

		v.clear();

		while( s ) {
			v.push_back( atof( s ) );
			s = strtok( NULL, ":;, " );
		}

		return true;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* GetArgListFromURL --------------------------------------------- */
/* --------------------------------------------------------------- */

#ifdef USE_CURL

static void Tokenize(
	vector<string>	&tokens,
	const string	&str,
	const string	&delimiters = "\n" )
{
	string::size_type	start, stop = 0;

	for(;;) {

		start = str.find_first_not_of( delimiters, stop );

		if( start == string::npos )
			return;

		stop = str.find_first_of( delimiters, start + 1 );

		tokens.push_back( str.substr( start, stop - start ) );
	}
}


static size_t WriteCallback( void *src, size_t size, size_t nmemb, void *dst )
{
	size_t	bytes = size * nmemb;

	((std::string*)dst)->append( (char*)src, bytes );

	return bytes;
}


#define	GETCOEFF( idx, name )											\
	if( tokens2[0].find( name ) != std::string::npos )					\
		{coeffs[idx] = strtod( tokens2[2].c_str(), NULL );}


// These URLs are json encodings wherein:
// - token[0] = name
// - token[1] = ':'
// - token[2] = value
//
bool GetArgListFromURL( vector<double> &v, const char *pat, char *argv )
{
	int	len = strlen( pat );

	if( !strncmp( argv, pat, len ) ) {

		char	*url = argv + len;

		if( !strstr( url, "http:" ) )
			return false;

		std::string	readBuffer;
		CURL		*easy_handle = curl_easy_init();
		CURLcode	res;

//		curl_easy_setopt( easy_handle, CURLOPT_VERBOSE, 1L );
		curl_easy_setopt( easy_handle, CURLOPT_URL, url );
		curl_easy_setopt( easy_handle, CURLOPT_WRITEFUNCTION, WriteCallback );
		curl_easy_setopt( easy_handle, CURLOPT_WRITEDATA, &readBuffer );

		res = curl_easy_perform( easy_handle );
		curl_easy_cleanup( easy_handle );

		vector<string>	tokens;
		double			coeffs[] = {0., 0., 0., 0., 0., 0.};

		Tokenize( tokens, readBuffer );

		for( int i = 0, n = tokens.size(); i < n; ++i )  {

			vector<string>	tokens2;
			Tokenize( tokens2, tokens[i], " " );

			GETCOEFF( 3, "imageRow" );
			GETCOEFF( 0, "imageCol" );
			GETCOEFF( 1, "minX" );
			GETCOEFF( 4, "minY" );
			GETCOEFF( 2, "stageX" );
			GETCOEFF( 5, "stageY" );
		}

		v.assign( coeffs, coeffs + 6 );
		return true;
	}

	return false;
}

#endif	// USE_CURL


