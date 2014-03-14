//
// lsq comprises two parts.
//
// (a) lsq.exe (log file 'lsq.txt') scans command line params
//		and looks at layers catalog to determine how the task
//		shall be divided among workers. Finally, it launches
//		those workers. MPI is used for multi-worker cases.
//
// (b) lsqw.exe (log file 'lsqw_i.txt') is the worker code.
//


#include	"lsq_Layers.h"

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"Maths.h"
#include	"Memory.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	double		Wr,				// Aff -> (1-Wr)*Aff + Wr*Rgd
				Etol;			// point error tolerance
	char		tempdir[2048],	// master workspace
				cachedir[2048];	// {catalog, pnts} files
	const char	*prior,			// start from these solutions
				*mode;			// {catalog,eval,split,A2A,A2H,H2H}
	int			zilo,			// my output range
				zihi,
				zolo,			// extended input range
				zohi,
				regtype,		// regularizer {T,R}
				iters,			// solve iterations
				splitmin,		// separate islands > splitmin tiles
				zpernode,		// max layers per node
				maxthreads;		// thr/node if not mpi
	bool		catclr,			// remake point catalog
				untwist;		// iff prior are affines

public:
	CArgs()
	{
		Wr			= 0.001;
		Etol		= 10.0;
		tempdir[0]	= 0;
		cachedir[0]	= 0;
		prior		= NULL;
		mode		= "A2A";
		zilo		= 0;
		zihi		= 0;
		zolo		= -1;
		zohi		= -1;
		regtype		= 'R';
		iters		= 2000;
		splitmin	= 1000;
		zpernode	= 200;
		maxthreads	= 1;
		catclr		= false;
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
		else if( GetArgStr( mode, "-mode=", argv[i] ) )
			printf( "Mode: '%s'.\n", mode );
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
		else if( GetArg( &Wr, "-Wr=T,%lf", argv[i] ) ) {
			regtype = 'T';
			printf( "Rglizer Wr: T, %g\n", Wr );
		}
		else if( GetArg( &Wr, "-Wr=R,%lf", argv[i] ) ) {
			regtype = 'R';
			printf( "Rglizer Wr: R, %g\n", Wr );
		}
		else if( GetArg( &Etol, "-Etol=%lf", argv[i] ) )
			printf( "Error  tol: %g\n", Etol );
		else if( GetArg( &iters, "-iters=%d", argv[i] ) )
			printf( "Iterations: %d\n", iters );
		else if( GetArg( &splitmin, "-splitmin=%d", argv[i] ) )
			;
		else if( GetArg( &zpernode, "-zpernode=%d", argv[i] ) )
			;
		else if( GetArg( &maxthreads, "-maxthreads=%d", argv[i] ) )
			;
		else if( IsArg( "-catclr", argv[i] ) )
			catclr = true;
		else if( IsArg( "-untwist", argv[i] ) )
			untwist = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// Mode valid?

	if( !strcmp( mode, "catalog" ) )
		goto mode_ok;
	if( !strcmp( mode, "eval" ) ) {
		iters = 0;
		goto mode_ok;
	}
	if( !strcmp( mode, "split" ) ) {
		iters = 0;
		goto mode_ok;
	}
	if( !strcmp( mode, "A2A" ) )
		goto mode_ok;
	if( !strcmp( mode, "A2H" ) )
		goto mode_ok;
	if( !strcmp( mode, "H2H" ) )
		goto mode_ok;

	printf( "Invalid -mode string [%s].\n", mode );
	exit( 42 );

// Default cache folder name

mode_ok:
	if( !cachedir[0] )
		DskAbsPath( cachedir, sizeof(cachedir), "lsqcache", stdout );

// Default zo range = zi range

	if( zolo == -1 ) {
		zolo = zilo;
		zohi = zihi;
	}

// Explicit starting solution needed if {solve stack, eval, split}

	if( !strcmp( mode, "catalog" ) )
		;
	else if( zilo != zihi
		|| !strcmp( mode, "eval" )
		|| !strcmp( mode, "split" ) ) {

		if( !prior ) {
			printf( "Missing -prior option.\n" );
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
static void ZoFromZi(
	int					&zolo,
	int					&zohi,
	int					zilo_icat,
	int					zihi_icat,
	const vector<Layer>	&vL )
{
	static int	maxspan = LayerCat_MaxSpan( vL );

	int	imax;

// zolo: lowest z that any interior layer touches

	zolo = vL[zilo_icat].Lowest();
	imax = min( zilo_icat + maxspan, zihi_icat );

	for( int icat = imax; icat > zilo_icat; --icat ) {

		int z = vL[icat].Lowest();

		if( z < zolo )
			zolo = z;
	}

// zohi: highest z that touches interior

	int	zihi = vL[zihi_icat].z;

	imax = min( zihi_icat + maxspan, vL.size() - 1 );

	for( int icat = imax; icat > zihi_icat; --icat ) {

		if( vL[icat].Lowest() <= zihi ) {

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
		nwks	= nL / zpernode;

	if( nL - nwks * zpernode > 0 )
		++nwks;

// Rebalance the load.

	if( nwks > 1 ) {

		int	zper = nL / nwks;

		if( nL - nwks * zper > 0 )
			++zper;

		zpernode = zper;
	}

	printf( "Workers %d, z-per-node %d.\n", nwks, zpernode );

// Launch the appropriate worker set.

	char	buf[2048];

	if( nwks <= 1 ) {

		// 1 worker

		if( maxthreads <= 1 ) { // pass flow in-process

			sprintf( buf,
			"lsqw -nwks=%d -temp=%s"
			" -cache=%s -prior=%s"
			" -mode=%s -Wr=%c,%g -Etol=%g -iters=%d"
			" -splitmin=%d -maxthreads=1"
			" -zi=%d,%d -zo=%d,%d"
			"%s",
			nwks, tempdir,
			cachedir, (prior ? prior : ""),
			mode, regtype, Wr, Etol, iters,
			splitmin,
			zilo, zihi, zolo, zohi,
			(untwist ? " -untwist" : "") );
		}
		else {	// qsub for desired slots

			sprintf( buf,
			"qsub -N lsqw -cwd -V -b y -pe batch %d"
			" lsqw -nwks=%d -temp=%s"
			" -cache=%s -prior=%s"
			" -mode=%s -Wr=%c,%g -Etol=%g -iters=%d"
			" -splitmin=%d -maxthreads=%d"
			" -zi=%d,%d -zo=%d,%d"
			"%s",
			maxthreads,
			nwks, tempdir,
			cachedir, (prior ? prior : ""),
			mode, regtype, Wr, Etol, iters,
			splitmin, maxthreads,
			zilo, zihi, zolo, zohi,
			(untwist ? " -untwist" : "") );
		}
	}
	else {

		// Write 'ranges.txt' telling each worker which
		// layers it's responsible for and which it needs.

		FILE	*f = FileOpenOrDie( "ranges.txt", "w" );
		int		zilo_icat = 0,
				zihi_icat = min( zpernode, nL ) - 1;

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
		" lsqw -nwks=%d -temp=%s"
		" -cache=%s -prior=%s"
		" -mode=%s -Wr=%c,%g -Etol=%g -iters=%d"
		" -splitmin=%d -maxthreads=16"
		"%s\n",
		nwks,
		nwks, tempdir,
		cachedir, (prior ? prior : ""),
		mode, regtype, Wr, Etol, iters,
		splitmin,
		(untwist ? " -untwist" : "") );
		fprintf( f, "\n" );

		fclose( f );
		FileScriptPerms( "mpigo.sht" );

		// Submit request to run 'mpigo.sht' script

		sprintf( buf, "qsub -N lsqmpi -cwd -V -b y -o sge.txt"
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

	if( !LayerCat( vL, gArgs.tempdir, gArgs.cachedir,
			gArgs.zolo, gArgs.zohi, gArgs.catclr ) ) {

		exit( 42 );
	}

	if( strcmp( gArgs.mode, "catalog" ) )
		gArgs.LaunchWorkers( vL );

	VMStats( stdout );

	return 0;
}


