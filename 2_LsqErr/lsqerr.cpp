//
// Read in one or two files (A,B) of type 'PostFitErrs.txt' as
// output from LSQ. Entries in these files are:
//
//	LyrA	TilA	RgnA	LyrB	TilB	RgnB	SqrErr
//
// Create output file 'lsqerrhst.txt' recording histogram of
// errors for:
//
//	Err	Aall	Asam	Adwn	Ball	Bsam	Bdwn
//
// All histograms are binwidth=1/div, nbins = div*lim + 1overflow.
//

#include	"Cmdline.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CHst {

public:
	int	*all, *sam, *dwn;

public:
	void GetFile( const char *name );
};

/* --------------------------------------------------------------- */
/* CArgs_lsqerr -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_lsqerr {

public:
	char	*inA, *inB;
	int		lim, div;

public:
	CArgs_lsqerr()
	{
		inA	= NULL;
		inB	= NULL;
		lim	= 200;
		div	= 1;
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
		else if( GetArg( &lim, "-lim=%d", argv[i] ) )
			;
		else if( GetArg( &div, "-div=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// headers

	fprintf( flog, "Err\tAall\tAsam\tAdwn" );

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
	int			emax	= gArgs.lim * gArgs.div;
	int			bytes	= (emax + 1)*sizeof(int);

	all = (int*)malloc( bytes );
	sam = (int*)malloc( bytes );
	dwn = (int*)malloc( bytes );

	memset( all, 0, bytes );
	memset( sam, 0, bytes );
	memset( dwn, 0, bytes );

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

		if( err < 0 ) {
			// -2 = not used
			// -1 = not inlier
			continue;
		}

		err = gArgs.div * sqrt( err );

		int	ibin = (err < emax ? int( err ) : emax);

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
	int	n = gArgs.lim * gArgs.div + 1;

	for( int i = 0; i < n; ++i ) {

		if( gArgs.inB ) {

			fprintf( flog,
			"%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n",
			double(i + 1)/gArgs.div,
			A.all[i], A.sam[i], A.dwn[i],
			B.all[i], B.sam[i], B.dwn[i] );
		}
		else {
			fprintf( flog,
			"%.2f\t%d\t%d\t%d\n",
			double(i + 1)/gArgs.div,
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


