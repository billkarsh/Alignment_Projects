//
// Reset all ThmPair_N^M.txt files to headers-only state
// in given range [zmin,zmax]. Current dir must be 'temp'
// at execution time.
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
/* CArgs_cln ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_cln {

public:
	int		zmin, zmax;

public:
	CArgs_cln() : zmin(-1), zmax(-1) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_cln	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_cln::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "cleanthmpair.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: cleanthmpair -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	if( zmin == -1 || zmax == -1 || zmin > zmax ) {
		printf( "Usage: cleanthmpair -zmin=i -zmax=j.\n" );
		printf( "Bad [zmin, zmax]=[%d, %d].\n", zmin, zmax );
		exit( 42 );
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* WriteFiles ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteOne( int a, int b )
{
	char	name[256];
	FILE	*f;

	sprintf( name, "%d/ThmPair_%d^%d.txt", a, a, b );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "Atl\tAcr\tBtl\tBcr\tErr\tDeg\tQ\tR"
	"\tT0\tT1\tX\tT3\tT4\tY\n" );

	fclose( f );
}


static void WriteFiles()
{
	for( int i = gArgs.zmin; i <= gArgs.zmax; ++i ) {

		WriteOne( i, i );

		if( i > 0 )
			WriteOne( i, i - 1 );
	}
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

/* ----------- */
/* Write files */
/* ----------- */

	WriteFiles();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


