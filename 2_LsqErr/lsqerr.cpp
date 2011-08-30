//
// Read in one or two files (A,B) of type 'PostFitErrs.txt' as
// output from LSQ. Entries in these files are:
//
//	LyrA	TilA	RgnA	LyrB	TilB	RgnB	SqrErr
//
// Create output file 'lsqerrhst.txt' recording histogram of
// errors for:
//
//	Aall	Asam	Adwn	Ball	Bsam	Bdwn
//
// All histograms are binwidth=1, 200 bins + overflow.
//

#include	"Cmdline.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	NMAX	200

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CHst {

public:
	int	all[NMAX+2],
		sam[NMAX+2],
		dwn[NMAX+2];

public:
	void GetFile( const char *name );
};

/* --------------------------------------------------------------- */
/* CArgs_lsqerr -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_lsqerr {

public:
	char	*inA, *inB;

public:
	CArgs_lsqerr()
	{
		inA	= NULL;
		inB	= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_lsqerr	gArgs;
static FILE*		flog	= NULL;
static CHst			A, B;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_lsqerr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "lsqerrhst.txt", "w" );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: lsqerr fileA [fileB] [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( !inA )
				inA = argv[i];
			else
				inB = argv[i];
		}
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// headers

	fprintf( flog, "Aall\tAsam\tAdwn" );

	if( inB )
		fprintf( flog, "\tBall\tBsam\tBdwn\n" );
	else
		fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* GetFile ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CHst::GetFile( const char *name )
{
	FILE		*f = FileOpenOrDie( name, "r" );
	CLineScan	LS;

	memset( all, 0, sizeof(all) );
	memset( sam, 0, sizeof(sam) );
	memset( dwn, 0, sizeof(dwn) );

	if( LS.Get( f ) <= 0 )
		goto close;

	for(;;) {

		double	err;
		int		za, ta, ra, zb, tb, rb;

		if( LS.Get( f ) <= 0 )
			break;

		sscanf( LS.line,
		"%d\t%d\t%d\t%d\t%d\t%d\t%lf",
		&za, &ta, &ra, &zb, &tb, &rb, &err );

		err = sqrt( err );

		int	ibin = (err < NMAX ? int( err ) : NMAX);

		++all[ibin];

		if( za == zb )
			++sam[ibin];
		else
			++dwn[ibin];
	}

close:
	fclose( f );
}

/* --------------------------------------------------------------- */
/* Record -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Record()
{
	for( int i = 0; i <= NMAX; ++i ) {

		if( gArgs.inB ) {

			fprintf( flog,
			"%d\t%d\t%d\t%d\t%d\t%d\n",
			A.all[i], A.sam[i], A.dwn[i],
			B.all[i], B.sam[i], B.dwn[i] );
		}
		else {
			fprintf( flog,
			"%d\t%d\t%d\n",
			A.all[i], A.sam[i], A.dwn[i] );
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

	A.GetFile( gArgs.inA );

	if( gArgs.inB )
		B.GetFile( gArgs.inB );

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


