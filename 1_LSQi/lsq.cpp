

#include	"lsq_Catalog.h"

#include	"Cmdline.h"
#include	"Disk.h"
#include	"Maths.h"
#include	"Memory.h"
#include	"Timer.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	char	tempdir[2048],	// master workspace
			*priorafftbl;	// start affine model from these
	int		wkid,			// my worker id (main=0)
			zilo,			// my output range
			zihi,
			zolo,			// extended input range
			zohi,
			zpernode;		// max layers per node
	bool	clrcat;			// remake point catalog

public:
	CArgs()
	{
		tempdir[0]	= 0;
		priorafftbl	= NULL;
		wkid		= 0;
		zilo		= 0;
		zihi		= 0;
		zolo		= -1;
		zohi		= -1;
		zpernode	= 200;
		clrcat		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs				gArgs;
static vector<CSubdirCat>	vCat;
static int					gnw = 1;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// Name log by worker

	for( int i = 1; i < argc; ++i ) {
		if( GetArg( &wkid, "-wkid=%d", argv[i] ) )
			break;
	}

	char slog[32];
	sprintf( slog, "lsq_%d.txt", wkid );
	freopen( slog, "w", stdout );

// Parse command line args

	printf( "---- Read params ----\n" );

	if( argc < 3 ) {
		printf(
		"Usage: lsq -temp=path -zi=i,j [options].\n" );
		exit( 42 );
	}

	vector<int>	vi;

	for( int i = 1; i < argc; ++i ) {

		const char	*instr;

		if( GetArgStr( instr, "-temp=", argv[i] ) ) {

			DskAbsPath( tempdir, sizeof(tempdir), instr, stdout );
			printf( "Temp dir: '%s'.\n", tempdir );
		}
		else if( GetArgStr( instr, "-prior=", argv[i] ) ) {

			priorafftbl = strdup( instr );
			printf( "Prior solutions: '%s'.\n", priorafftbl );
		}
		else if( GetArg( &wkid, "-wkid=%d", argv[i] ) )
			printf( "wkid %d\n", wkid );
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
		else if( GetArg( &zpernode, "-zpernode=%d", argv[i] ) )
			;
		else if( IsArg( "-clrcat", argv[i] ) )
			clrcat = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	if( zolo == -1 ) {
		zolo = zilo;
		zohi = zihi;
	}
	
	if( !wkid && zilo != zihi ) {

		if( !priorafftbl ) {
			printf( "Solving a stack requires prior affines.\n" );
			exit( 42 );
		}

		if( !DskExists( priorafftbl ) ) {
			printf( "Prior affines not found [%s].\n", priorafftbl );
			exit( 42 );
		}
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* ZoFromZi ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Determine dependency range [zolo,zohi] and return cat index
// of zohi (so master can truncate it's list).
//
static int ZoFromZi(
	int	&zolo,
	int	&zohi,
	int	zilo_icat,
	int	zihi_icat )
{
	zolo = *vCat[zilo_icat].zdown.begin();

// Look up to 10 layers away (12 for safety)

	int	zihi = vCat[zihi_icat].z;
	int	imax = min( zihi_icat + 12, vCat.size() - 1 );

	for( int icat = imax; icat > zihi_icat + 1; --icat ) {

		if( *vCat[icat].zdown.begin() <= zihi ) {

			zohi = vCat[icat].z;
			return icat;
		}
	}

// Default

	zohi = zihi;
	return zihi_icat;
}

/* --------------------------------------------------------------- */
/* MasterLaunchWorkers ------------------------------------------- */
/* --------------------------------------------------------------- */

static void MasterLaunchWorkers()
{
// How many workers?

	int	nc = vCat.size();

	gnw = nc / gArgs.zpernode;

	if( nc - gnw * gArgs.zpernode > 0 )
		++gnw;

	if( gnw <= 1 )
		return;

// Master will be lowest block

	int	zolo,
		zohi,
		zilo_icat = 0,
		zihi_icat = gArgs.zpernode - 1,
		newcatsiz = ZoFromZi( zolo, zohi, zilo_icat, zihi_icat ) + 1;

		gArgs.zihi = vCat[zihi_icat].z;
		gArgs.zohi = zohi;

		printf( "\nMaster own range zi [%d %d] zo [%d %d]\n",
		gArgs.zilo, gArgs.zihi, gArgs.zolo, gArgs.zohi );

// Make ranges

	for( int iw = 1; iw < gnw; ++iw ) {

		zilo_icat = zihi_icat + 1;
		zihi_icat = min( zilo_icat + gArgs.zpernode, nc ) - 1;
		ZoFromZi( zolo, zohi, zilo_icat, zihi_icat );

		char	buf[1024];

		sprintf( buf, "qsub -N lsq-%d -cwd -V -b y -pe batch 16"
		" 'lsq -temp=%s -prior=%s -wkid=%d -zi=%d,%d -zo=%d,%d'",
		iw,
		gArgs.tempdir, gArgs.priorafftbl,
		iw,
		vCat[zilo_icat].z, vCat[zihi_icat].z,
		zolo, zohi );

		system( buf );
	}

	vCat.resize( newcatsiz );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// This flow is followed by master process (gArgs.wkid==0)
// and all workers process.
//
int main( int argc, char **argv )
{
/* ---------------------- */
/* All parse command line */
/* ---------------------- */

	gArgs.SetCmdLine( argc, argv );

/* ---------------------- */
/* All need their catalog */
/* ---------------------- */

	printf( "---- Cataloging ----\n" );

	clock_t	t0 = StartTiming();

	CatPoints( vCat, gArgs.tempdir,
		gArgs.zolo, gArgs.zohi, gArgs.clrcat );

	t0 = StopTiming( stdout, "Catalog", t0 );

/* ------------------------ */
/* Master partitions layers */
/* ------------------------ */

	if( !gArgs.wkid )
		MasterLaunchWorkers();

/* ---- */
/* Done */
/* ---- */

	VMStats( stdout );
	return 0;
}


