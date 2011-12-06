//
// Read CTFormTable.txt file and create table of pairwise XY-offsets.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"CTForm.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Entry {

public:
	int		z, id;
	TForm	T;
};


class CTbl {

public:
	vector<Entry>	vE;

public:
	void GetFile( const char *name );
};

/* --------------------------------------------------------------- */
/* CArgs_xytbl -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xytbl {

public:
	char	*infile;
	int		z;

public:
	CArgs_xytbl()
	{
		infile	= NULL;
		z		= -1;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xytbl	gArgs;
static FILE*		flog	= NULL;
static CTbl			A;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xytbl::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "xytable.txt", "w" );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: xytable fileA [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &z, "-z=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// headers

	fprintf( flog, "Za\tTa\tZb\tTb\tX\tY\n" );
}

/* --------------------------------------------------------------- */
/* GetFile ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTbl::GetFile( const char *name )
{
	FILE		*f = FileOpenOrDie( name, "r" );
	CLineScan	LS;

	if( LS.Get( f ) <= 0 )
		goto close;

	for(;;) {

		Entry	E;
		int		rgn;

		if( LS.Get( f ) <= 0 )
			break;

		sscanf( LS.line,
		"%d\t%d\t%d"
		"\t%lf\t%lf\t%lf"
		"\t%lf\t%lf\t%lf",
		&E.z, &E.id, &rgn,
		&E.T.t[0], &E.T.t[1], &E.T.t[2],
		&E.T.t[3], &E.T.t[4], &E.T.t[5] );

		if( gArgs.z >= 0 && E.z != gArgs.z )
			continue;

		vE.push_back( E );
	}

close:
	fclose( f );
}

/* --------------------------------------------------------------- */
/* Record -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Record()
{
	int	n = A.vE.size();

	for( int i = 0; i < n - 1; ++i ) {

		const Entry&	ei = A.vE[i];

		for( int j = i + 1; j < n; ++j ) {

			const Entry&	ej = A.vE[j];

			TForm	ab;

			AToBTrans( ab, ei.T, ej.T );

			if( fabs( ab.t[2] ) > 4096 ||
				fabs( ab.t[5] ) > 4096 ) {

				continue;
			}

			fprintf( flog,
			"%d\t%d\t%d\t%d"
			"\t%f\t%f\n",
			ei.z, ei.id, ej.z, ej.id,
			ab.t[2], ab.t[5] );
		}
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char **argv )
{
/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ---- */
/* Read */
/* ---- */

	A.GetFile( gArgs.infile );

/* ----- */
/* Print */
/* ----- */

	Record();

/* ---- */
/* Done */
/* ---- */

exit:
	fclose( flog );

	return 0;
}


