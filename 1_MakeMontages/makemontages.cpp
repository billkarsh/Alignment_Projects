//
// makemontages reads an IDB 'image database' and creates an
// alignment workspace with this structure:
//
//	folder 'alnname'				// top folder
//		imageparams.txt				// IDBPATH, IMAGESIZE tags
//		folder '0'					// folder per layer, here, '0'
//			folder '0'				// output folder per tile, here '0'
//			S0_0					// same layer jobs
//				make.same			// make file for same layer
//				ThmPair_0^0.txt		// table of thumbnail results
//			D0_0					// down layer jobs
//				make.down			// make file for cross layers
//				ThmPair_0^j.txt		// table of thumbnail results
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"Geometry.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Pair {
public:
	int	a, b;
public:
	Pair( int a, int b ) : a(a), b(b) {};
};


class Block {
public:
	vector<Pair>	P;
};


class BlockSet {

private:
	enum {
		klowcount = 12
	};

public:
	vector<Block>	K;
	int				w, h,
					kx, ky,
					dx, dy,
					nb;

private:
	void OrientLayer( int is0, int isN );
	void SetDims();
	void PartitionJobs( int is0, int isN );
	void Consolidate();
	void ReportBlocks( int z );

public:
	void CarveIntoBlocks( int is0, int isN );
	void MakeJobs( const char *lyrdir, int z );
};

/* --------------------------------------------------------------- */
/* CArgs_scr ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scr {

public:
	// minareafrac used for same and cross layer.
	// Typical vals:
	// 0.025	historic typical
	// 0.020	Davi and new typical
	// 0.0003	Davi older Harvard data
	//
	double		minareafrac;
	string		idbpath;
	const char	*outdir,
				*exenam;
	int			zmin,
				zmax,
				blksize,
				xml_type,
				xml_min,
				xml_max;
	bool		NoFolds,
				NoDirs,
				davinocorn;

public:
	CArgs_scr()
	{
		minareafrac	= 0.020;
		outdir		= "NoSuch";	// prevent overwriting real dir
		exenam		= "ptest";
		zmin		= 0;
		zmax		= 32768;
		blksize		= 8;
		xml_type	= 0;
		xml_min		= 0;
		xml_max		= 0;
		NoFolds		= false;
		NoDirs		= false;
		davinocorn	= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scr	gArgs;
static CTileSet		TS;
static FILE*		flog	= NULL;
static char			xmlprms[256]	= {0};
static int			gW		= 0,	// universal pic dims
					gH		= 0;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "makemontages.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make alignment workspace: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: makemontages <idbpath> -d=temp -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			idbpath = argv[i];
		else if( GetArgStr( outdir, "-d=", argv[i] ) )
			;
		else if( GetArgStr( exenam, "-exe=", argv[i] ) )
			;
		else if( GetArg( &minareafrac, "-minareafrac=%lf", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &blksize, "-b=%d", argv[i] ) )
			;
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			;
		else if( GetArg( &xml_min, "-xmlmin=%d", argv[i] ) )
			;
		else if( GetArg( &xml_max, "-xmlmax=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else if( IsArg( "-nd", argv[i] ) )
			NoDirs = true;
		else if( IsArg( "-davinc", argv[i] ) )
			davinocorn = true;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
	char	name[2048];

// create the top dir
	DskCreateDir( gArgs.outdir, flog );

// copy imageparams here
	sprintf( name, "cp %s/imageparams.txt %s",
		gArgs.idbpath.c_str(), gArgs.outdir );
	system( name );

// create stack subdir
	sprintf( name, "%s/stack", gArgs.outdir );
	DskCreateDir( name, flog );

// create mosaic subdir
	//sprintf( name, "%s/mosaic", gArgs.outdir );
	//DskCreateDir( name, flog );
}

/* --------------------------------------------------------------- */
/* WriteRunlsqFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void _WriteRunlsqFile( const char *path, int z, bool final )
{
	FILE	*f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Assign each tile (or each subregion if tiles divided by folds)\n" );
	fprintf( f, "# a transform (default is affine) that best describes the mapping\n" );
	fprintf( f, "# of correspondence point pairs from local images to the shared\n" );
	fprintf( f, "# global system.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > lsq -temp=temp0 -zi=i,j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Required:\n" );
	fprintf( f, "# -temp=path\t\t;master workspace (a.k.a. 'temp')\n" );
	fprintf( f, "# -zi=i,j\t\t\t;output tforms in range z=[i..j]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -cache=lsqcache\t;path to {catalog,pnts} data\n" );
	fprintf( f, "# -catclr\t\t\t;rebuild catalog file\n" );
	fprintf( f, "# -zo=p,q\t\t\t;consider input out to z=[p..q]\n" );
	fprintf( f, "# -prior=path\t\t;starting tforms (required if stack)\n" );
	fprintf( f, "# -untwist\t\t\t;untwist prior affines\n" );
	fprintf( f, "# -mode=A2A\t\t\t;action: {catalog,eval,split,A2A,A2H,H2H}\n" );
	fprintf( f, "# -Wr=R,0.001\t\t;Aff -> (1-Wr)*Aff + Wr*(T=Trans, R=Rgd}\n" );
	fprintf( f, "# -Etol=30\t\t\t;max point error (depends upon system size)\n" );
	fprintf( f, "# -iters=2000\t\t;solve iterations\n" );
	fprintf( f, "# -splitmin=1000\t;separate islands > splitmin tiles\n" );
	fprintf( f, "# -zpernode=200\t\t;max layers per cluster node\n" );
	fprintf( f, "# -maxthreads=1\t\t;thr/node if not mpi (16 if mpi)\n" );
	fprintf( f, "# -local\t\t\t;run locally (no qsub) if 1 worker\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );

	if( !final )
		fprintf( f, "lsq -temp=../../ -zi=%d,%d $1\n", z, z );
	else {
		fprintf( f, "lsq -temp=../ -zi=%d,%d"
		" -prior=../cross_wkspc/X_A_BIN_scaf -untwist"
		" -mode=A2A -Wr=R,0 -Etol=500 -iters=10000"
		" -zpernode=200 -maxthreads=16\n",
		gArgs.zmin, gArgs.zmax );
	}

	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}


static void WriteRunlsqFile()
{
	char	buf[2048];
	sprintf( buf, "%s/stack/runlsq.sht", gArgs.outdir );
	_WriteRunlsqFile( buf, -1, true );
}

/* --------------------------------------------------------------- */
/* WriteXviewFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteXviewFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/stack/xviewgo.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Convert Z-file tform data to viewable text or xml.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Z-file organized data can be:\n" );
	fprintf( f, "# - IDB (inpath not specified)\n" );
	fprintf( f, "# - X_A_TXT, X_H_TXT\n" );
	fprintf( f, "# - X_A_MET, X_H_MET\n" );
	fprintf( f, "# - X_A_BIN, X_H_BIN\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > xview inpath -idb=idbpath -z=i,j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Required:\n" );
	fprintf( f, "# inPath\t\t\t;path of a tform source folder\n" );
	fprintf( f, "# -idb=idbpath\t\t;path to idb folder\n" );
	fprintf( f, "# -z=i,j\t\t\t;convert tforms in range z=[i..j]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# ImagePlus type codes used in xml files:\n" );
	fprintf( f, "#\tAUTO\t\t= -1\n" );
	fprintf( f, "#\tGRAY8\t\t= 0\n" );
	fprintf( f, "#\tGRAY16\t\t= 1\n" );
	fprintf( f, "#\tGRAY32\t\t= 2\n" );
	fprintf( f, "#\tCOLOR_256\t= 3\n" );
	fprintf( f, "#\tCOLOR_RGB\t= 4\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -forceWH=w,h\t\t;force bounds (override degcw)\n" );
	fprintf( f, "# -degcw=0\t\t\t;rotate CW degrees\n" );
	fprintf( f, "# -type=X\t\t\t;{T,M,X,B} = {X_?_TXT,X_?_MET,xml,billfile}\n" );
	fprintf( f, "# -meta=path\t\t;billfile with alt image paths\n" );
	fprintf( f, "# -xmltrim=0.0\t\t;trim this much from xml images\n" );
	fprintf( f, "# -xmltype=0\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "xview X_A_BIN -idb=%s -z=%d,%d -degcw=0 -type=X%s\n",
	gArgs.idbpath.c_str(), gArgs.zmin, gArgs.zmax, xmlprms );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteEviewFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteEviewFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/stack/eviewgo.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Histogram one or two 'Error' folders produced by lsqw.\n" );
	fprintf( f, "# The result is a text file for viewing in Excel.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > eview Error [Error_B] -z=i,j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Required:\n" );
	fprintf( f, "# Error\t\t\t;path to lsqw Error folder\n" );
	fprintf( f, "# -z=i,j\t\t;use data in range z=[i..j]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# [Error]\t\t;second Error folder for comparison\n" );
	fprintf( f, "# -div=10\t\t;bin width = 1/div\n" );
	fprintf( f, "# -lim=500\t\t;nbins = div*lim + 1 (for oflo)\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "eview Error -z=%d,%d -div=10 -lim=100\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteCountsamedirsFile ---------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCountsamedirsFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/countsamedirs.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Count all 'Sx_y' dirs in layer range\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./countsamedirs.sht i j\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "cnt=0\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcnt=$(($cnt+1))\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\techo z= $lyr  cum= $cnt\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSSubNFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSSubNFile( int njobs )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/ssub.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, submit all make.same and use the make option -j <n>\n" );
	fprintf( f, "# to set number of concurrent jobs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./ssub.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "nthr=%d\n", njobs );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcd $jb\n" );
	fprintf( f, "\t\t\tqsub -N q$jb-$lyr -cwd -o /dev/null -V -b y -pe batch $nthr make -f make.same -j $nthr EXTRA='\"\"'\n" );
	fprintf( f, "\t\t\tcd ..\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteDSubNFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteDSubNFile( int njobs )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/dsub.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, submit all make.down and use the make option -j <n>\n" );
	fprintf( f, "# to set number of concurrent jobs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./dsub.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "nthr=%d\n", njobs );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcd $jb\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\t\tif [ -e make.down ]\n" );
	fprintf( f, "\t\t\tthen\n" );
	fprintf( f, "\t\t\t\tqsub -N q$jb-$lyr -cwd -o /dev/null -V -b y -pe batch $nthr make -f make.down -j $nthr EXTRA='\"\"'\n" );
	fprintf( f, "\t\t\tfi\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\t\tcd ..\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportFiles ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportFiles()
{
	char	buf[2048];
	FILE	*f;

// same

	sprintf( buf, "%s/sreport.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For montages in layer range...\n" );
	fprintf( f, "# Tabulate sizes of all cluster stderr logs for quick view of faults.\n" );
	fprintf( f, "# Tabulate sizes of all 'pts.same' files for consistency checking.\n" );
	fprintf( f, "# Tabulate subblocks for which there were no points.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./sreport.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l */S*/qS*.e* > SameErrs.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l */S*/pts.same > SamePts.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -f SameNopts.txt\n" );
	fprintf( f, "touch SameNopts.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tfor jb in $(ls -d $lyr/* | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tif [ ! -e $jb/pts.same ]\n" );
	fprintf( f, "\t\t\tthen\n" );
	fprintf( f, "\t\t\t\techo \"$jb\" >> SameNopts.txt\n" );
	fprintf( f, "\t\t\tfi\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );

// down

	sprintf( buf, "%s/dreport.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For down data in layer range...\n" );
	fprintf( f, "# Tabulate sizes of all cluster stderr logs for quick view of faults.\n" );
	fprintf( f, "# Tabulate sizes of all 'pts.down' files for consistency checking.\n" );
	fprintf( f, "# Tabulate subblocks for which there were no points.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./dreport.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l */D*/qD*.e* > DownErrs.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l */D*/pts.down > DownPts.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -f DownNopts.txt\n" );
	fprintf( f, "touch DownNopts.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tfor jb in $(ls -d $lyr/* | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tif [ -e $jb/make.down -a ! -e $jb/pts.down ]\n" );
	fprintf( f, "\t\t\tthen\n" );
	fprintf( f, "\t\t\t\techo \"$jb\" >> DownNopts.txt\n" );
	fprintf( f, "\t\t\tfi\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteMSubFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMSubFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/msub.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For each layer in range, cd to montage dir, run lsq there.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./msub.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr/montage\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tqsub -N mon-$lyr -cwd -V -b y -pe batch 1 \"./runlsq.sht\"\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ../..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteMReportFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMReportFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/mreport.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, gather list of 'FINAL' lines from lsq.txt files\n" );
	fprintf( f, "# for individual montages. Report these in MonSumy.txt.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./mreport.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -rf MonSumy.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tlog=$lyr/montage/lsqw_0.txt\n" );
	fprintf( f, "\tif [ -f \"$log\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\techo Z $lyr `grep -e \"FINAL*\" $log` >> MonSumy.txt\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteGatherMonsFile ------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteGatherMonsFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/gathermons.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For each layer in range, copy montage results to X_A_BIN_mons.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./gathermons.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "dst=X_A_BIN_mons\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -rf $dst\n" );
	fprintf( f, "mkdir -p $dst\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcp $lyr/montage/X_A_BIN/* $dst\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* Write_Crossgo ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_Crossgo()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/crossgo.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Given folder of Z-file montage tforms (srcmons),\n" );
	fprintf( f, "# create workspace mytemp/cross_wkspc, and scripts\n" );
	fprintf( f, "# governing cross layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_topscripts -zmin=i -zmax=j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -nf\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -scl=50\t\t\t;size reduction factor\n" );
	fprintf( f, "# -lgord=1\t\t\t;Legendre poly max order\n" );
	fprintf( f, "# -sdev=42\t\t\t;scape sdev size\n" );
	fprintf( f, "# -resmask\t\t\t;mask out resin\n" );
	fprintf( f, "# -abwide=15\t\t;strip width in tiles\n" );
	fprintf( f, "# -stpcorr=0.02\t\t;req. min strip corr\n" );
	fprintf( f, "# -blkmincorr=0.45\t;required min corr for alignment\n" );
	fprintf( f, "# -blknomcorr=0.50\t;nominal corr for alignment\n" );
	fprintf( f, "# -xyconf=0.75\t\t;search radius = (1-conf)(blockwide)\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_topscripts X_A_BIN_mons -zmin=%d -zmax=%d%s -scl=50 -lgord=1 -sdev=42 -abwide=15 -stpcorr=0.02 -blkmincorr=0.45 -blknomcorr=0.50 -xyconf=0.75\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : "") );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* CreateLayerDir ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Each layer gets a directory named by its z-index. All content
// pertains to this layer, or to this layer acting as a source
// onto itself or other layers.
//
// For example, make.down will contain ptest jobs aligning this
// layer onto that below (z-1).
//
static void CreateLayerDir( char *lyrdir, int L )
{
	fprintf( flog, "\n\nCreateLayerDir: layer %d\n", L );

// Create layer dir
	sprintf( lyrdir, "%s/%d", gArgs.outdir, L );
	DskCreateDir( lyrdir, flog );

// Create montage subdir
	char	buf[2048];
	int		len;

	len = sprintf( buf, "%s/montage", lyrdir );
	DskCreateDir( buf, flog );

// Create montage script
	sprintf( buf + len, "/runlsq.sht" );
	_WriteRunlsqFile( buf, L, false );
}

/* --------------------------------------------------------------- */
/* CreateTileSubdirs --------------------------------------------- */
/* --------------------------------------------------------------- */

// Each tile gets a directory named by its picture id. All content
// pertains to this tile, or this tile acting as a source onto
// other tiles.
//
// For example, folder 8/10 contains the foldmask fm.png for tile
// 10 in layer 8. If this folder contains file 7.11.tf.txt it
// lists the transforms mapping tile 8/10 onto tile 7/11.
//
static void CreateTileSubdirs( const char *lyrdir, int is0, int isN )
{
	fprintf( flog, "--CreateTileSubdirs: layer %d\n",
		TS.vtil[is0].z );

	for( int i = is0; i < isN; ++i ) {

		char	subdir[2048];

		sprintf( subdir, "%s/%d", lyrdir, TS.vtil[i].id );
		DskCreateDir( subdir, flog );
	}
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P).
//
static void WriteMakeFile(
	const char			*lyrdir,
	int					SD,
	int					ix,
	int					iy,
	const vector<Pair>	&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/%c%d_%d/make.%s",
	lyrdir, SD, ix, iy, (SD == 'S' ? "same" : "down") );

	f = FileOpenOrDie( name, "w", flog );

// write 'all' targets line

	fprintf( f, "all: " );

	for( int i = 0; i < np; ++i ) {

		const CUTile&	A = TS.vtil[P[i].a];
		const CUTile&	B = TS.vtil[P[i].b];

		fprintf( f, "%d/%d.%d.map.tif ", A.id, B.z, B.id );
	}

	fprintf( f, "\n\n" );

// Write each 'target: dependencies' line
//		and each 'rule' line

	const char	*option_nf = (gArgs.NoFolds ? " -nf" : "");

	for( int i = 0; i < np; ++i ) {

		const CUTile&	A = TS.vtil[P[i].a];
		const CUTile&	B = TS.vtil[P[i].b];

		fprintf( f,
		"%d/%d.%d.map.tif:\n",
		A.id, B.z, B.id );

		fprintf( f,
		"\t%s %d.%d^%d.%d%s ${EXTRA}\n\n",
		gArgs.exenam, A.z, A.id, B.z, B.id, option_nf );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* OrientLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rotate layer to have smallest footprint.
//
void BlockSet::OrientLayer( int is0, int isN )
{
// Collect all tile corners in (C)

	vector<Point>	C, cnr;
	Set4Corners( cnr, gW, gH );

	for( int i = is0; i < isN; ++i ) {

		vector<Point>	c( 4 );
		memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
		TS.vtil[i].T.Transform( c );

		for( int i = 0; i < 4; ++i )
			C.push_back( c[i] );
	}

// Rotate layer upright and translate to (0,0)

	TAffine	R;
	DBox	B;
	int		deg = TightestBBox( B, C );

	R.NUSetRot( deg*PI/180 );

	for( int i = is0; i < isN; ++i ) {

		TAffine&	T = TS.vtil[i].T;

		T = R * T;
		T.AddXY( -B.L, -B.B );
	}

	w = int(B.R - B.L) + 1;
	h = int(B.T - B.B) + 1;
}

/* --------------------------------------------------------------- */
/* SetDims ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Divide layer bbox into (kx,ky) cells of size (dx,dy).
//
void BlockSet::SetDims()
{
	dx = gArgs.blksize;
	dy = dx;
	dx *= gW;
	dy *= gH;
	kx = (int)ceil( (double)w / dx );
	ky = (int)ceil( (double)h / dy );
	dx = w / kx;
	dy = h / ky;
	nb = kx * ky;
}

/* --------------------------------------------------------------- */
/* PartitionJobs ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fill pair jobs into blocks according to a-tile center.
//
void BlockSet::PartitionJobs( int is0, int isN )
{
	K.clear();
	K.resize( nb );

	int	W2 = gW/2,
		H2 = gH/2;

	for( int a = is0; a < isN; ++a ) {

		Point	pa( W2, H2 );
		int		ix, iy, rowa, cola;

		TS.vtil[a].T.Transform( pa );

		ix = int(pa.x / dx);

		if( ix < 0 )
			ix = 0;
		else if( ix >= kx )
			ix = kx - 1;

		iy = int(pa.y / dy);

		if( iy < 0 )
			iy = 0;
		else if( iy >= ky )
			iy = ky - 1;

		if( gArgs.davinocorn ) {

			const char	*c, *n;

			n = FileNamePtr( TS.vtil[a].name.c_str() );
			c = strstr( n, "col" );
			sscanf( c, "col%d_row%d", &cola, &rowa );
		}

		for( int b = a + 1; b < isN; ++b ) {

			if( gArgs.davinocorn ) {

				const char	*c, *n;
				int			rowb, colb;

				n = FileNamePtr( TS.vtil[b].name.c_str() );
				c = strstr( n, "col" );
				sscanf( c, "col%d_row%d", &colb, &rowb );

				if( (rowb - rowa) * (colb - cola) != 0 )
					continue;
			}

			if( TS.ABOlap( a, b ) > gArgs.minareafrac )
				K[ix + kx*iy].P.push_back( Pair( a, b ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* Consolidate --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Distribute small job sets into low occupancy neighbors.
//
void BlockSet::Consolidate()
{
	if( nb <= 1 )
		return;

	bool	changed;

	do {

		changed = false;

		for( int i = 0; i < nb; ++i ) {

			int	ic = K[i].P.size();

			if( !ic || ic >= klowcount )
				continue;

			int	iy = i / kx,
				ix = i - kx * iy,
				lowc = 0,
				lowi, c;

			// find lowest count neib

			if( iy > 0 && (c = K[i-kx].P.size()) ) {
				lowc = c;
				lowi = i-kx;
			}

			if( iy < ky-1 && (c = K[i+kx].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i+kx;
			}

			if( ix > 0 && (c = K[i-1].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i-1;
			}

			if( ix < kx-1 && (c = K[i+1].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i+1;
			}

			// merge

			if( !lowc )
				continue;

			changed = true;

			for( int j = 0; j < ic; ++j )
				K[lowi].P.push_back( K[i].P[j] );

			K[i].P.clear();
		}

	} while( changed );
}

/* --------------------------------------------------------------- */
/* ReportBlocks -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print job count array.
//
void BlockSet::ReportBlocks( int z )
{
	int	njobs = 0;

	fprintf( flog, "\nZ %d, Array %dx%d, Jobs(i,j):\n", z, kx, ky );

	for( int i = 0; i < nb; ++i ) {

		int	iy = i / kx,
			ix = i - kx * iy,
			ij = K[i].P.size();

		fprintf( flog, "%d%c", ij, (ix == kx - 1 ? '\n' : '\t') );
		njobs += ij;
	}

	fprintf( flog, "Total = %d\n", njobs );
}

/* --------------------------------------------------------------- */
/* CarveIntoBlocks ----------------------------------------------- */
/* --------------------------------------------------------------- */

void BlockSet::CarveIntoBlocks( int is0, int isN )
{
	OrientLayer( is0, isN );
	SetDims();
	PartitionJobs( is0, isN );
	Consolidate();
	ReportBlocks( TS.vtil[is0].z );
}

/* --------------------------------------------------------------- */
/* MakeJobs ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void BlockSet::MakeJobs( const char *lyrdir, int z )
{
	for( int i = 0; i < nb; ++i ) {

		if( K[i].P.size() ) {

			int	iy = i / kx,
				ix = i - kx * iy;

			CreateJobsDir( lyrdir, ix, iy, z, z, flog );
			WriteMakeFile( lyrdir, 'S', ix, iy, K[i].P );
		}
	}
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer()
{
	int		is0, isN;

	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 ) {

		char		lyrdir[2048];
		BlockSet	BS;
		int			z = TS.vtil[is0].z;

		CreateLayerDir( lyrdir, z );

		if( !gArgs.NoDirs )
			CreateTileSubdirs( lyrdir, is0, isN );

		BS.CarveIntoBlocks( is0, isN );
		BS.MakeJobs( lyrdir, z );

		TS.GetLayerLimits( is0 = isN, isN );
	}
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

	TS.SetLogFile( flog );

	int	L = 0;

	if( gArgs.xml_type != 0 )
		L += sprintf( xmlprms + L, " -xmltype=%d", gArgs.xml_type );

	if( gArgs.xml_min != 0 )
		L += sprintf( xmlprms + L, " -xmlmin=%d", gArgs.xml_min );

	if( gArgs.xml_max != 0 )
		L += sprintf( xmlprms + L, " -xmlmax=%d", gArgs.xml_max );

/* ---------------- */
/* Read source data */
/* ---------------- */

	TS.FillFromIDB( gArgs.idbpath, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", (int)TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SetTileDimsFromIDB( gArgs.idbpath );
	TS.GetTileDims( gW, gH );

/* ------------------------------------------------- */
/* Within each layer, sort tiles by dist from center */
/* ------------------------------------------------- */

	TS.SortAll_z_r();

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteRunlsqFile();
	WriteXviewFile();
	WriteEviewFile();

	WriteCountsamedirsFile();
	WriteSSubNFile( 4 );
	WriteDSubNFile( 4 );
	WriteReportFiles();
	WriteMSubFile();
	WriteMReportFile();
	WriteGatherMonsFile();
	Write_Crossgo();

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


