

#include	"lsq_Bounds.h"
#include	"lsq_Error.h"
#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"
#include	"lsq_MPI.h"
#include	"lsq_Untwist.h"
#include	"lsq_XArray.h"

#include	"Cmdline.h"
#include	"File.h"
#include	"Memory.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*tempdir,		// master workspace
				*cachedir,		// {catalog, pnts} files
				*prior;			// start from these solutions
	int			zilo,			// my output range
				zihi,			// iff nwks == 1
				zolo,			// extended input range
				zohi;
	bool		untwist;		// iff prior are affines

public:
	CArgs()
	{
		tempdir		= NULL;
		cachedir	= NULL;
		prior		= NULL;
		zilo		= 0;
		zihi		= 0;
		zolo		= -1;
		zohi		= -1;
		untwist		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
	void GetRanges();
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// Parse command line args

	vector<int>	vi;

	for( int i = 1; i < argc; ++i ) {

		if( GetArgStr( tempdir, "-temp=", argv[i] ) ) {

			printf( "Temp  dir: '%s'.\n", tempdir );
			GetIDB( tempdir );
		}
		else if( GetArgStr( cachedir, "-cache=", argv[i] ) )
			printf( "Cache dir: '%s'.\n", cachedir );
		else if( GetArgStr( prior, "-prior=", argv[i] ) )
			printf( "Prior solutions: '%s'.\n", prior );
		else if( GetArg( &nwks, "-nwks=%d", argv[i] ) )
			;
		else if( GetArgList( vi, "-zi=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zilo = vi[0];
				zihi = vi[1];
				printf( "zi [%d %d]\n", zilo, zihi );
			}
			else {
				printf( "Bad format in -zi [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else if( GetArgList( vi, "-zo=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zolo = vi[0];
				zohi = vi[1];
				printf( "zo [%d %d]\n", zolo, zohi );
			}
			else {
				printf( "Bad format in -zo [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else if( IsArg( "-untwist", argv[i] ) )
			untwist = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}
}

/* --------------------------------------------------------------- */
/* GetRanges ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read my zi and zo ranges from file if nwks > 1.
//
// Also read the overlap range data for workers
// to my immediate left and right. This is used
// in XArray::Updt for solution coherency.
//
void CArgs::GetRanges()
{
	if( nwks <= 1 )
		return;

	FILE		*f = FileOpenOrDie( "ranges.txt", "r" );
	CLineScan	LS;
	int			code = 0;	// {4=LHS, 2=me, 1=RHS}

// Need entries for my left, for me, and my right,
// unless I am the far left or right in which case
// we preset that success code.

	if( !wkid )
		code += 4;
	else if( wkid == nwks - 1 )
		code += 1;

	while( code != 7 && LS.Get( f ) > 0 ) {

		int	iw = atoi( LS.line );

		if( iw < wkid - 1 )
			continue;
		else if( iw == wkid - 1 ) {	// LHS

			int	dum1, dum2;

			sscanf( LS.line, "%d zi=%d,%d zo=%d,%d",
				&iw, &dum1, &zLlo, &dum2, &zLhi );

			++zLlo;
			code += 4;
		}
		else if( iw == wkid ) {	// me

			sscanf( LS.line, "%d zi=%d,%d zo=%d,%d",
				&iw, &zilo, &zihi, &zolo, &zohi );

			printf( "zi [%d %d]\n", zilo, zihi );
			printf( "zo [%d %d]\n", zolo, zohi );

			code += 2;
		}
		else if( iw == wkid + 1 ) {	// RHS

			int	dum1, dum2;

			sscanf( LS.line, "%d zi=%d,%d zo=%d,%d",
				&iw, &zRhi, &dum1, &zRlo, &dum2 );

			--zRhi;
			code += 1;
			break;
		}
	}

	fclose( f );

	if( code != 7 ) {
		printf( "Ranges.txt: Missing entries for wkid %d.\n", wkid );
		MPIExit();
		exit( 42 );
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char **argv )
{
	clock_t	t0 = StartTiming();

/* ---------- */
/* Parameters */
/* ---------- */

	MPIInit( argc, argv );
	gArgs.SetCmdLine( argc, argv );
	gArgs.GetRanges();

/* ------------ */
/* Initial data */
/* ------------ */

	LayerCat( vL, gArgs.tempdir, gArgs.cachedir,
		gArgs.zolo, gArgs.zohi, false );

	InitTables( gArgs.zilo, gArgs.zihi );

	{
		CLoadPoints	*LP = new CLoadPoints;
		LP->Load( gArgs.tempdir, gArgs.cachedir );
		delete LP;
	}

/* ----- */
/* Start */
/* ----- */

	printf( "\n---- Development ----\n" );

	XArray	A;
	A.Load( gArgs.prior );

/* --------- */
/* Summaries */
/* --------- */

	Error( A );

	if( !wkid ) {
		double	erms, emax;
		GetFinalError( erms, emax );
		printf( "\nFINAL RMS %.2f MAX %.2f\n", erms, emax );
	}

	DBox B;
	Bounds( B, A );

	if( !wkid ) {
		printf( "\n" );
		StopTiming( stdout, "Lsq", t0 );
	}

/* ------- */
/* Cleanup */
/* ------- */

	MPIExit();
	VMStats( stdout );

	return 0;
}


