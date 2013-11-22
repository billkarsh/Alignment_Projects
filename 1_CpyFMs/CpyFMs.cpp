//
// Copy manually edited fold masks to std alignment hierarchy.
//

#include	"Cmdline.h"
#include	"File.h"

#include	<string.h>
#include	<time.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_ldir ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_ldir {

public:
	const char	*infile,
				*outdir;

public:
	CArgs_ldir() : infile(NULL), outdir("NoSuch") {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_ldir	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_ldir::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "CpyFMs.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: CpyFMs <txt-file> -d=<workdir>.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( outdir, "-d=", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );

	freopen( "CpyFMs.log", "a", stdout );
	freopen( "CpyFMs.log", "a", stderr );
}

/* --------------------------------------------------------------- */
/* Process ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Process()
{
	FILE		*f = FileOpenOrDie( gArgs.infile, "r", flog );
	CLineScan	LS;

	for(;;) {

		char	path[2048], cmd[2048];
		int		z, id;

		if( LS.Get( f ) <= 0 )
			break;

		sscanf( LS.line, "%d;%d;%[^;]", &z, &id, path );

		sprintf( cmd, "cp '%s' '%s/%d/%d/fm.tif'",
			path, gArgs.outdir, z, id );

		fprintf( flog, "%s\n", cmd );

		system( cmd );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ------- */
/* Process */
/* ------- */

	Process();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


