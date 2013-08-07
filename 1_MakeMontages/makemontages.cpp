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
//				ThmPair_0_@_0.txt	// table of thumbnail results
//			D0_0					// down layer jobs
//				make.down			// make file for cross layers
//				ThmPair_0_@_j.txt	// table of thumbnail results
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"Geometry.h"


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
	Pair( int _a, int _b )	{a = _a; b = _b;};
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
	double	minareafrac;
	string	idbpath;
	char	*outdir,
			*exenam;
	int		zmin,
			zmax,
			blksize,
			xml_type,
			xml_min,
			xml_max;
	bool	lens,
			NoFolds,
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
		xml_type	= -999;
		xml_min		= -999;
		xml_max		= -999;
		lens		= false;
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
		else if( IsArg( "-lens", argv[i] ) )
			lens = true;
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

static void _WriteRunlsqFile( const char *path, bool final )
{
	FILE	*f = FileOpenOrDie( path, "w", flog );

	char	xprms[256] = "";
	int		L = 0;

	if( gArgs.lens )
		L += sprintf( xprms + L, "-lens " );

	if( gArgs.xml_type != -999 )
		L += sprintf( xprms + L, "-xmltype=%d ", gArgs.xml_type );

	if( gArgs.xml_min != -999 )
		L += sprintf( xprms + L, "-xmlmin=%d ", gArgs.xml_min );

	if( gArgs.xml_max != -999 )
		L += sprintf( xprms + L, "-xmlmax=%d ", gArgs.xml_max );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Assign each tile (or each subregion if tiles divided by folds)\n" );
	fprintf( f, "# a transform (default is affine) that best describes the mapping of\n" );
	fprintf( f, "# correspondence points from local images to the shared global system.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > lsq pts.all [options] > lsq.txt\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -model=A\t\t\t;{T=translat,S=simlty,A=affine,H=homography}\n" );
	fprintf( f, "# -minmtglinks=1\t;min montage neibs/tile\n" );
	fprintf( f, "# -all\t\t\t\t;override min pts/tile-pair\n" );
	fprintf( f, "# -davinc\t\t\t;no davi bock same lyr corners\n" );
	fprintf( f, "# -same=1.0\t\t\t;wt same lyr vs cross layer\n" );
	fprintf( f, "# -scale=0.1\t\t;wt for unit magnitude constraint\n" );
	fprintf( f, "# -square=0.1\t\t;wt for equal magnitude cosines,sines constraint\n" );
	fprintf( f, "# -scaf=0.01\t\t;wt for external scaffold solution\n" );
	fprintf( f, "# -tformtol=0.2\t\t;max dev of Tform rot-elems from median\n" );
	fprintf( f, "# -threshold=700.0\t;max inlier error\n" );
	fprintf( f, "# -pass=1\t\t\t;num ransac passes\n" );
	fprintf( f, "# -degcw=0.0\t\t;rotate clockwise degrees (for appearance)\n" );
	fprintf( f, "# -lrbt=a,b,c,d\t\t;forced BBOX, default natural\n" );
	fprintf( f, "# -unite=2,path\t\t;unite blocks using this layer of this TFormTable\n" );
	fprintf( f, "# -prior=path\t\t;start affine model with this TFormTable\n" );
	fprintf( f, "# -nproc=1\t\t\t;num processors to use\n" );
	fprintf( f, "# -viserr=10\t\t;create color coded error stack; yellow value\n" );
	fprintf( f, "# -xmltype=0\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t;intensity scale\n" );
	fprintf( f, "# -strings\t\t\t;using old-style path-labeled data items\n" );
	fprintf( f, "# -p=/N\t\t\t\t;tilename id pattern (for -strings)\n" );
	fprintf( f, "# -trim=0.0\t\t\t;extra margin around each tile\n" );
	fprintf( f, "# -refz=0\t\t\t;deprecated\n" );
	fprintf( f, "# -lens\t\t\t\t;apply external affine software lens\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );

	if( final ) {
		fprintf( f, "lsq pts.all -scale=.1 -square=.1 -nproc=1 -prior=../cross_wkspc/HiRes_TF.txt %s$1 > lsq.txt\n",
		xprms );
	}
	else {
		fprintf( f, "lsq pts.all -scale=.1 -square=.1 %s$1 > lsq.txt\n",
		xprms );
	}

	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}


static void WriteRunlsqFile()
{
	char	buf[2048];
	sprintf( buf, "%s/stack/runlsq.sht", gArgs.outdir );
	_WriteRunlsqFile( buf, true );
}

/* --------------------------------------------------------------- */
/* WriteSubmosFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmosFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/mosaic/submos.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Heal montage seams and update superpixel data.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > mos <simple-file> <xmin,ymin,xsize,ysize> <zmin,zmax> [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Default sizing is 0,0,-1,-1 meaning natural size.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -d\t\t\t\t;debug\n" );
	fprintf( f, "# -strings\t\t\t;simple-file data labeled by path strings\n" );
	fprintf( f, "# -warp\t\t\t\t;heal seams\n" );
	fprintf( f, "# -nf\t\t\t\t;no folds\n" );
	fprintf( f, "# -a\t\t\t\t;annotate\n" );
	fprintf( f, "# -tiles\t\t\t;make raveler tiles\n" );
	fprintf( f, "# -noflat\t\t\t;no flat image ('before')\n" );
	fprintf( f, "# -nomap\t\t\t;no map image (where data from)\n" );
	fprintf( f, "# -matlab\t\t\t;matlab/closeness order\n" );
	fprintf( f, "# -drn\t\t\t\t;don't renumber superpixels\n" );
	fprintf( f, "# -dms=0.01\t\t\t;don't move strength\n" );
	fprintf( f, "# -fold_dir=path\t;prepended fm location, default=CWD\n" );
	fprintf( f, "# -region_dir=path\t;results go here, default=CWD\n" );
	fprintf( f, "# -gray_dir=path\t;gray images go here, default=CWD\n" );
	fprintf( f, "# -grey_dir=path\t;gray images go here, default=CWD\n" );
	fprintf( f, "# -sp_dir=path\t\t;superpixel maps go here, default=CWD\n" );
	fprintf( f, "# -inv_dir=path\t\t;inverse maps go here, default=CWD\n" );
	fprintf( f, "# -rav_dir=path\t\t;raveler tiles go here, default=CWD\n" );
	fprintf( f, "# -bmap_dir=path\t;boundary maps go here, default=CWD\n" );
	fprintf( f, "# -s=1\t\t\t\t;scale down by this integer\n" );
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
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tqsub -N mos-$lyr -cwd -V -b y -pe batch 8 \"mos ../stack/simple 0,0,-1,-1 $lyr,$lyr -warp%s > mos_$lyr.txt\"\n",
	(gArgs.NoFolds ? " -nf" : "") );
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

	sprintf( buf, "%s/ssub%d.sht", gArgs.outdir, njobs );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, submit all make.same and use the make option -j <n>\n" );
	fprintf( f, "# to set number of concurrent jobs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./ssub%d.sht <zmin> [zmax]\n",
	njobs );
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
	fprintf( f, "\t\t\tqsub -N q$jb-$lyr -cwd -V -b y -pe batch 8 make -f make.same -j %d EXTRA='\"\"'\n",
	njobs );
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

	sprintf( buf, "%s/dsub%d.sht", gArgs.outdir, njobs );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, submit all make.down and use the make option -j <n>\n" );
	fprintf( f, "# to set number of concurrent jobs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./dsub%d.sht <zmin> [zmax]\n",
	njobs );
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
	fprintf( f, "\t\t\t\tqsub -N q$jb-$lyr -cwd -V -b y -pe batch 8 make -f make.down -j -j %d EXTRA='\"\"'\n",
	njobs );
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
/* WriteMontage1File --------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMontage1File()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/montage1.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For given layer, gather all point pair files into montage/pts.all,\n" );
	fprintf( f, "# cd there, run lsq solver.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./montage1.sht <z>\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cd $1\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -f pts.all\n" );
	fprintf( f, "\n" );
	fprintf( f, "# get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
	fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < ../imageparams.txt)\n" );
	fprintf( f, "\n" );
	fprintf( f, "cp ../imageparams.txt pts.all\n" );
	fprintf( f, "cat $idb/$1/fm.same >> pts.all\n" );
	fprintf( f, "\n" );
	fprintf( f, "for jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tcat $jb/pts.same >> pts.all\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );
	fprintf( f, "mv pts.all montage\n" );
	fprintf( f, "\n" );
	fprintf( f, "cd montage\n" );
	fprintf( f, "./runlsq.sht \"\"\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmonFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmonFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/submon.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For each layer in range, gather points, cd to layer's montage dir,\n" );
	fprintf( f, "# run lsq there.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./submon.sht <zmin> [zmax]\n" );
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
	fprintf( f, "\t\tqsub -N mon-$lyr -cwd -V -b y -pe batch 8 \"./montage1.sht $lyr\"\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportMonsFile ------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportMonsFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sumymons.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, gather list of 'FINAL' lines from lsq.txt files\n" );
	fprintf( f, "# for individual montages. Report these in MonSumy.txt.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./sumymons.sht <zmin> [zmax]\n" );
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
	fprintf( f, "\tlog=$lyr/montage/lsq.txt\n" );
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
/* WriteCombineFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCombineFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/combine.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Gather all point pair files into stack/pts.all\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./combine.sht <zmin> <zmax>\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "rm -f pts.all\n" );
	fprintf( f, "\n" );
	fprintf( f, "#get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
	fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < imageparams.txt)\n" );
	fprintf( f, "\n" );
	fprintf( f, "cp imageparams.txt pts.all\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $2)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tcat $idb/$lyr/fm.same >> pts.all\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $2)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcat $jb/pts.same >> ../pts.all\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tif (($lyr != $1))\n" );
	fprintf( f, "\t\tthen\n" );
	fprintf( f, "\t\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\t\tdo\n" );
	fprintf( f, "\t\t\t\tcat $jb/pts.down >> ../pts.all\n" );
	fprintf( f, "\t\t\tdone\n" );
	fprintf( f, "\t\tfi\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );
	fprintf( f, "mv pts.all stack\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteFinishFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFinishFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/finish.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Gather all points to stack folder, cd there, run lsq solver.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# ./finish.sht\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "./combine.sht %d %d\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "cd stack\n" );
	fprintf( f, "./runlsq.sht \"\"\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSFinishFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSFinishFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sfinish.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Run finish.sht script on its own cluster node.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./sfinish.sht\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "qsub -N finish -cwd -V -b y -pe batch 8 ./finish.sht\n" );
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
	_WriteRunlsqFile( buf, false );
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
		"\t%s %d/%d@%d/%d%s ${EXTRA}\n\n",
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
// Collect all tile corners

	vector<Point>	C;

	for( int i = is0; i < isN; ++i ) {

		vector<Point>	p( 4 );

		p[0] = Point( 0.0 , 0.0 );
		p[1] = Point( gW-1, 0.0 );
		p[2] = Point( gW-1, gH-1 );
		p[3] = Point( 0.0 , gH-1 );

		TS.vtil[i].T.Transform( p );

		for( int i = 0; i < 4; ++i )
			C.push_back( p[i] );
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

/* ---------------- */
/* Read source data */
/* ---------------- */

	TS.FillFromIDB( gArgs.idbpath, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

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
	//WriteSubmosFile();

//	WriteSSubNFile( 4 );
	WriteSSubNFile( 8 );
//	WriteDSubNFile( 4 );
	WriteDSubNFile( 8 );
	WriteReportFiles();
	WriteMontage1File();
	WriteSubmonFile();
	WriteReportMonsFile();
	WriteCombineFile();
	WriteFinishFile();
	WriteSFinishFile();

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


