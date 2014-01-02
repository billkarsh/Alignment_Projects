

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"GenDefs.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*meta;

public:
	CArgs() : meta(NULL) {};

	void SetCmdLine( int argc, char* argv[] );
};

class Meta {
public:
	int	Hz,
		nchan;
public:
	Meta() : Hz(0), nchan(0) {};
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static char		path[2048];
static char		*name;
static Meta		M;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// Parse command line args

	if( argc < 2 ) {
		printf( "Usage: ephystxt <meta-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			meta = argv[i];
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}
}

/* --------------------------------------------------------------- */
/* ScanMeta ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static int ScanMeta()
{
	char	buf[2048];
	sprintf( buf, "%s/%s.meta", path, name );
	FILE		*f = FileOpenOrDie( buf, "r" );
	CLineScan	LS;

	while( LS.Get( f ) > 0 ) {

		if( LS.line[0] == 'n' ) {

			if( 1 == sscanf( LS.line, "nChans = %d", &M.nchan ) )
				printf( "num channels: %d\n", M.nchan );
		}

		if( LS.line[0] == 's' ) {

			if( 1 == sscanf( LS.line, "sRateHz = %d", &M.Hz ) )
				printf( "Hertz: %d\n", M.Hz );
		}
	}

	fclose( f );

	return (M.Hz != 0) && (M.nchan != 0);
}

/* --------------------------------------------------------------- */
/* Convert ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Convert()
{
	const int maxsamples = 20000;

	char	buf[2048];
	sprintf( buf, "%s/%s.bin", path, name );

	double	filebytes = DskBytes( buf );
	int		sampbytes = sizeof(int16) * M.nchan,
			nsamp = int(filebytes / sampbytes),
			readbytes;

	if( nsamp > maxsamples )
		nsamp = maxsamples;

	vector<char>	D( readbytes = nsamp * sampbytes );
	FILE			*f;

	f = FileOpenOrDie( buf, "rb" );
	fread( &D[0], 1, readbytes, f );
	fclose( f );

	sprintf( buf, "%s.txt", name );
	f = FileOpenOrDie( buf, "w" );

	for( int is = 0; is < nsamp; ++is ) {

		fprintf( f, "%.4e", double(is)/M.Hz );

		char*	s0 = &D[is * sampbytes];

		for( int ic = 0; ic < M.nchan; ++ic ) {
			fprintf( f, "\t%d", *(int16*)(s0 + sizeof(int16)*ic) );
		}

		fprintf( f, "\n" );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	gArgs.SetCmdLine( argc, argv );

// Path and basename

	strcpy( path, gArgs.meta );
	char *s = strrchr( path, '/' );
	if( s )
		*s = 0;
	else {
		path[0] = '.';
		path[1] = 0;
	}

	name = FileCloneNamePart( gArgs.meta );

// Decode meta

	if( !ScanMeta() ) {
		printf( "Missing fields in meta file.\n" );
		return 0;
	}

	Convert();

	return 0;
}


