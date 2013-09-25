//
// Merge TAffineTable and idb into new Bill file.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	char	*inpath,
			*idb;

public:
	CArgs()
	{
		inpath	= NULL;
		idb		= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static FILE*	flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "TAffineTable.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: taffinetable table idb\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !inpath )
				inpath = argv[i];
			else
				idb = argv[i];
		}
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Merge --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Merge()
{
	FILE		*fi = FileOpenOrDie( gArgs.inpath, "r" );
	FILE		*fo = FileOpenOrDie( "TAffineTableEx.txt", "w" );
	CLineScan	LS;
	string		idb = gArgs.idb;

	while( LS.Get( fi ) > 0 ) {

		double	t[6];
		int		z, id, rgn;

		sscanf( LS.line, "%d\t%d\t%d"
		"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&z, &id, &rgn,
		&t[0], &t[1], &t[2], &t[3], &t[4], &t[5] );

		const Til2Img	*t2i;

		if( !IDBT2ICacheNGet1( t2i, idb, z, id, flog ) )
			continue;

		fprintf( fo, "%d\t%d"
		"\t%f\t%f\t%f\t%f\t%f\t%f"
		"\t%d\t%d\t%d\t%s\n",
		z, id,
		t[0], t[1], t[2], t[3], t[4], t[5],
		t2i->col, t2i->row, t2i->cam, t2i->path.c_str() );
	}

	fclose( fo );
	fclose( fi );
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

/* ----- */
/* Merge */
/* ----- */

	Merge();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


