//
// lsq comprises two parts.
//
// (a) lsq.exe (log file 'lsq.txt') scans command line params
//		and looks at layers catalog to determine how the task
//		shall be divided among workers. Finally, it launches
//		those workers. If only one worker is needed, control
//		is passed to lsqw by system(). Otherwise, a distributed
//		node-set is configured for mpi data sharing.
//
// (b) lsqw.exe (log file 'lsqw_i.txt') is the worker code.
//


#include	"lsq_Layers.h"

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"Maths.h"
#include	"Memory.h"


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	char		tempdir[2048],	// master workspace
				cachedir[2048];	// {catalog, pnts} files
	const char	*prior;			// start from these solutions
	int			zilo,			// my output range
				zihi,
				zolo,			// extended input range
				zohi,
				zpernode;		// max layers per node
	bool		catclr,			// remake point catalog
				catonly,
				untwist;		// iff prior are affines

public:
	CArgs()
	{
		tempdir[0]	= 0;
		cachedir[0]	= 0;
		prior		= NULL;
		zilo		= 0;
		zihi		= 0;
		zolo		= -1;
		zohi		= -1;
		zpernode	= 200;
		catclr		= false;
		catonly		= false;
		untwist		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
	void LaunchWorkers( const vector<Layer> &vL );
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
// Name front end log

	freopen( "lsq.txt", "w", stdout );

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
			printf( "Temp  dir: '%s'.\n", tempdir );
		}
		else if( GetArgStr( instr, "-cache=", argv[i] ) ) {

			DskAbsPath( cachedir, sizeof(cachedir), instr, stdout );
			printf( "Cache dir: '%s'.\n", cachedir );
		}
		else if( GetArgStr( prior, "-prior=", argv[i] ) )
			printf( "Prior solutions: '%s'.\n", prior );
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
		else if( IsArg( "-catclr", argv[i] ) )
			catclr = true;
		else if( IsArg( "-catonly", argv[i] ) )
			catonly = true;
		else if( IsArg( "-untwist", argv[i] ) )
			untwist = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// Default cahe folder name

	if( !cachedir[0] )
		DskAbsPath( cachedir, sizeof(cachedir), "lsqcache", stdout );

// Default zo range = zi range

	if( zolo == -1 ) {
		zolo = zilo;
		zohi = zihi;
	}

// Stacks need explicit starting solution

	if( zilo != zihi ) {

		if( !prior ) {
			printf( "Solving a stack requires -prior option.\n" );
			exit( 42 );
		}

		if( !DskExists( prior ) ) {
			printf( "Priors not found [%s].\n", prior );
			exit( 42 );
		}
	}
}

/* --------------------------------------------------------------- */
/* ZoFromZi ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Determine dependency range [zolo,zohi].
//
// Look at most 10 layers from each end (12 for safety).
//
static void ZoFromZi(
	int					&zolo,
	int					&zohi,
	int					zilo_icat,
	int					zihi_icat,
	const vector<Layer>	&vL )
{
	int	imax;

// zolo: lowest z that any interior layer touches

	zolo = vL[zilo_icat].Lowest();
	imax = min( zilo_icat + 12, vL.size() - 1 );

	for( int icat = imax; icat > zilo_icat; --icat ) {

		int z = vL[icat].Lowest();

		if( z < zolo )
			zolo = z;
	}

// zohi: highest z that touches interior

	int	zihi = vL[zihi_icat].z;

	imax = min( zihi_icat + 12, vL.size() - 1 );

	for( int icat = imax; icat > zihi_icat; --icat ) {

		if( *vL[icat].zdown.begin() <= zihi ) {

			zohi = vL[icat].z;
			return;
		}
	}

// Default

	zohi = zihi;
}

/* --------------------------------------------------------------- */
/* LaunchWorkers ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::LaunchWorkers( const vector<Layer> &vL )
{
	printf( "\n---- Launching workers ----\n" );

// How many workers?

	int	nL		= vL.size(),
		nwks	= nL / zpernode,
		res		= nL - nwks * zpernode;

	if( res > 0 )
		++nwks;

// Now balance the load a little. If nwks > 1 then at this point,
// all workers but the last have zpernode layers, and we want to
// check that the last isn't much smaller than that. If it is, we
// shift layers to the last by taking an equal number from each
// of the others.

	if( nwks > 1 ) {

		int	fromeach = (zpernode - res) / (nwks - 1);

		zpernode -= fromeach;
	}

	printf( "Workers %d, z-per-node %d.\n", nwks, zpernode );

// Launch the appropriate worker set.

	char	buf[2048];

	if( nwks <= 1 ) {

		// 1 worker: pass flow in process to lsqw.

		sprintf( buf,
		"lsqw -nwks=%d -temp=%s -cache=%s -prior=%s"
		" -zi=%d,%d -zo=%d,%d"
		"%s",
		nwks, tempdir, cachedir, prior,
		zilo, zihi, zolo, zohi,
		(untwist ? " -untwist" : "") );
	}
	else {

		// Write 'ranges.txt' telling each worker which
		// layers it's responsible for and which it needs.

		FILE	*f = FileOpenOrDie( "ranges.txt", "w" );
		int		zilo_icat = 0,
				zihi_icat = zpernode - 1;

		for( int iw = 0; iw < nwks; ++iw ) {

			ZoFromZi( zolo, zohi, zilo_icat, zihi_icat, vL );

			fprintf( f, "%d zi=%d,%d zo=%d,%d\n",
			iw,  vL[zilo_icat].z,  vL[zihi_icat].z, zolo, zohi );

			zilo_icat = zihi_icat + 1;
			zihi_icat = min( zilo_icat + zpernode, nL ) - 1;
		}

		fclose( f );

		// Write script 'mpigo.sht' to:
		// (1) clean up host list 'sge.txt' --> 'hosts.txt'.
		// (2) call mpirun using 'hosts.txt' file.

		f = FileOpenOrDie( "mpigo.sht", "w" );

		fprintf( f, "#!/bin/sh\n" );
		fprintf( f, "\n" );
		fprintf( f, "tail -n +2 sge.txt > hosts.txt\n" );
		fprintf( f, "\n" );
		fprintf( f, "mpirun -perhost 1 -n %d -machinefile hosts.txt"
		" lsqw -nwks=%d -temp=%s -cache=%s -prior=%s"
		"%s\n",
		nwks,
		nwks, tempdir, cachedir, prior,
		(untwist ? " -untwist" : "") );
		fprintf( f, "\n" );

		fclose( f );
		FileScriptPerms( "mpigo.sht" );

		// Submit request to run 'mpigo.sht' script

		sprintf( buf, "qsub -N lsqw -cwd -V -b y -o sge.txt"
		" -pe impi3 %d ./mpigo.sht", 16 * nwks );
	}

	printf( "Launch <>\n<%s>\n", buf );

	system( buf );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char **argv )
{
	gArgs.SetCmdLine( argc, argv );

	vector<Layer>	vL;

	LayerCat( vL, gArgs.tempdir, gArgs.cachedir,
		gArgs.zolo, gArgs.zohi, gArgs.catclr );

	if( !gArgs.catonly )
		gArgs.LaunchWorkers( vL );

	VMStats( stdout );

	return 0;
}


