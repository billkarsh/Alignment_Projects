

#include	"lsq_Bounds.h"
#include	"lsq_Dropout.h"
#include	"lsq_Error.h"
#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"
#include	"lsq_Magnitude.h"
#include	"lsq_MPI.h"
#include	"lsq_Solve.h"
#include	"lsq_Split.h"
#include	"lsq_Untwist.h"

#include	"Cmdline.h"
#include	"File.h"
#include	"Memory.h"
#include	"Timer.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*tempdir,		// master workspace
				*cachedir,		// {catalog, pnts} files
				*prior,			// start from these solutions
				*mode;			// {catalog,eval,split,A2A,A2H,H2H}
	int			zilo,			// my output range
				zihi,
				zolo,			// extended input range
				zohi,
				iters,			// solve iterations
				splitmin;		// separate islands > splitmin tiles
	bool		untwist;		// iff prior are affines

public:
	CArgs()
	{
		tempdir		= NULL;
		cachedir	= NULL;
		prior		= NULL;
		mode		= NULL;
		zilo		= 0;
		zihi		= 0;
		zolo		= -1;
		zohi		= -1;
		iters		= 2000;
		splitmin	= 10;
		untwist		= false;
	};

	bool SetCmdLine( int argc, char* argv[] );
	bool GetRanges();
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CArgs::SetCmdLine( int argc, char* argv[] )
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
		else if( GetArgStr( mode, "-mode=", argv[i] ) )
			printf( "Mode: '%s'.\n", mode );
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
				return false;
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
				return false;
			}
		}
		else if( GetArg( &iters, "-iters=%d", argv[i] ) )
			printf( "Iterations: %d\n", iters );
		else if( GetArg( &splitmin, "-splitmin=%d", argv[i] ) )
			printf( "Split-min:  %d\n", iters );
		else if( GetArg( &maxthreads, "-maxthreads=%d", argv[i] ) )
			printf( "Maxthreads: %d\n", maxthreads );
		else if( IsArg( "-untwist", argv[i] ) )
			untwist = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			return false;
		}
	}

	return true;
}

/* --------------------------------------------------------------- */
/* GetRanges ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read my zi and zo ranges from file if nwks > 1.
//
// Also read the overlap range data for workers
// to my immediate left and right. This is used
// for cross-machine updates.
//
bool CArgs::GetRanges()
{
	if( nwks <= 1 )
		return true;

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

			zLlo;	// InitTables() must set layer after this
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

			zRhi;	// InitTables() must set layer before this
			code += 1;
			break;
		}
	}

	fclose( f );

	if( code != 7 ) {
		printf( "Ranges.txt: Missing entries for wkid %d.\n", wkid );
		return false;
	}

	return true;
}

/* --------------------------------------------------------------- */
/* Evaluate ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void Evaluate( const XArray &X )
{
	DBox B;
	Bounds( B, X );

	Magnitude( X );
	Error( X );

	Dropout	D;
	D.Scan();

	if( !wkid ) {

		double	erms, emax;
		GetFinalError( erms, emax );

		printf(
		"\nFINAL RMS %.2f MAX %.2f"
		" {Mx,Rd,Pt,Ki,Cd} %ld %ld %ld %ld %ld\n",
		erms, emax,
		D.rmax, D.read, D.pnts, D.kill, D.cutd );
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

	if( !gArgs.SetCmdLine( argc, argv ) ||
		!gArgs.GetRanges() ) {

		MPIExit();
		exit( 42 );
	}

/* ------------ */
/* Initial data */
/* ------------ */

	if( !LayerCat( vL, gArgs.tempdir, gArgs.cachedir,
			gArgs.zolo, gArgs.zohi, false ) ) {

		MPIExit();
		exit( 42 );
	}

	InitTables( gArgs.zilo, gArgs.zihi );

	{
		CLoadPoints	*LP = new CLoadPoints;
		LP->Load( gArgs.tempdir, gArgs.cachedir );
		delete LP;
	}

/* ----- */
/* Solve */
/* ----- */

	printf( "\n---- Solve ----\n" );

	XArray	Xevn, Xodd;

	if( !strcmp( gArgs.mode, "A2A" ) ) {

		Xevn.Load( gArgs.prior );

		if( gArgs.untwist )
			UntwistAffines( Xevn );

		Xodd.Resize( 6 );
		Solve( Xevn, Xodd, gArgs.iters );
	}
	else if( !strcmp( gArgs.mode, "A2H" ) ) {

		Xevn.Resize( 8 );

		{	// limit A lifetime
			XArray	*A = new XArray;
			A->Load( gArgs.prior );

			if( gArgs.untwist )
				UntwistAffines( *A );

			Solve( *A, Xevn, 1 );
			delete A;
		}

		Xodd.Resize( 8 );
		Solve( Xevn, Xodd, gArgs.iters );
	}
	else if( !strcmp( gArgs.mode, "H2H" ) ) {

		Xevn.Load( gArgs.prior );
		Xodd.Resize( 8 );
		Solve( Xevn, Xodd, gArgs.iters );
	}
	else if( !strcmp( gArgs.mode, "eval" ) ) {

		Xevn.Load( gArgs.prior );

		if( gArgs.untwist )
			UntwistAffines( Xevn );

		gArgs.iters = 0;
	}
	else {	// split

		Xevn.Load( gArgs.prior );
		gArgs.iters = 0;
	}

	const XArray& Xfinal = ((gArgs.iters & 1) ? Xodd : Xevn);

/* ----------- */
/* Postprocess */
/* ----------- */

	if( gArgs.mode[0] == 's' ) {

		Split	S( Xfinal, gArgs.splitmin );

		S.Run();
	}
	else {

		Evaluate( Xfinal );

		if( gArgs.mode[0] != 'e' )
			Xfinal.Save();
	}

/* ------- */
/* Cleanup */
/* ------- */

	if( !wkid ) {
		printf( "\n" );
		StopTiming( stdout, "Lsq", t0 );
	}

	MPIExit();
	VMStats( stdout );

	return 0;
}


