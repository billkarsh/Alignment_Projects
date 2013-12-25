

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"
#include	"lsq_MPI.h"
#include	"lsq_Untwist.h"
#include	"lsq_XArray.h"

#include	"Cmdline.h"
#include	"File.h"
#include	"Memory.h"


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*tempdir,		// master workspace
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

			printf( "Temp dir: '%s'.\n", tempdir );
			GetIDB( tempdir );
		}
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

// Read zi an zo ranges from file if nwks > 1.
//
void CArgs::GetRanges()
{
	if( nwks <= 1 )
		return;

	FILE		*f = FileOpenOrDie( "ranges.txt", "r" );
	CLineScan	LS;

	while( LS.Get( f ) > 0 ) {

		if( atoi( LS.line ) == wkid ) {

			sscanf( LS.line, "%d zi=%d,%d zo=%d,%d",
				&wkid, &zilo, &zihi, &zolo, &zohi );

			printf( "zi [%d %d]\n", zilo, zihi );
			printf( "zo [%d %d]\n", zolo, zohi );

			fclose( f );
			return;
		}
	}

	printf( "Ranges.txt: No entry for wkid %d.\n", wkid );
	MPIExit();
	exit( 42 );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char **argv )
{
/* ---------- */
/* Parameters */
/* ---------- */

	MPIInit( argc, argv );
	gArgs.SetCmdLine( argc, argv );
	gArgs.GetRanges();

/* ------------ */
/* Initial data */
/* ------------ */

	LayerCat( vL, gArgs.tempdir, gArgs.zolo, gArgs.zohi, false );

	InitTables( gArgs.zilo, gArgs.zihi );

	{
		CLoadPoints	*LP = new CLoadPoints;
		LP->Load( gArgs.tempdir );
		delete LP;
	}

/* ----- */
/* Start */
/* ----- */

	printf( "\n---- Development ----\n" );

	XArray	A;
	A.Load( gArgs.prior );
	UntwistAffines( A );
	A.Save();

/* ------- */
/* Cleanup */
/* ------- */

	MPIExit();
	VMStats( stdout );

	return 0;
}


