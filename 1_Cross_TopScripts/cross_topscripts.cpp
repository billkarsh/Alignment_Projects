//
// Write scripts governing cross layer alignment.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"TrakEM2_UTL.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_cross --------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_cross {

public:
	double	stpcorr,
			blkmincorr,
			blknomcorr,
			xyconf;		// search radius = (1-conf)(blockwide)
	char	mondir[2048];
	int		zmin,
			zmax,
			scl,
			lgord,
			sdev,
			abwide;
	bool	resmask,
			NoFolds;

public:
	CArgs_cross()
	{
		stpcorr		= 0.02;
		blkmincorr	= 0.45;
		blknomcorr	= 0.50;
		xyconf		= 0.75;
		zmin		= 0;
		zmax		= 32768;
		scl			= 50;
		lgord		= 1;	// 1  probably good for Davi EM
		sdev		= 42;	// 42 useful for Davi EM
		abwide		= 15;
		resmask		= false;
		NoFolds		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_cross	gArgs;
static string		idb;
static char			*gtopdir = "cross_wkspc";
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_cross::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "cross_topscripts.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make scapeops scripts: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: cross_topscripts -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		const char	*_outdir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &scl, "-scl=%d", argv[i] ) )
			;
		else if( GetArg( &lgord, "-lgord=%d", argv[i] ) )
			;
		else if( GetArg( &sdev, "-sdev=%d", argv[i] ) )
			;
		else if( GetArg( &abwide, "-abwide=%d", argv[i] ) )
			;
		else if( GetArg( &stpcorr, "-stpcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &blkmincorr, "-blkmincorr=%lf", argv[i] ) )
			;
		else if( GetArg( &blknomcorr, "-blknomcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &xyconf, "-xyconf=%lf", argv[i] ) ) {

			if( xyconf < 0.0 || xyconf > 1.0 )
				xyconf = 0.5;
		}
		else if( IsArg( "-resmask", argv[i] ) )
			resmask = true;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );

	DskAbsPath( mondir, sizeof(mondir), "X_A_BIN_mons", flog );

	if( !DskExists( mondir ) ) {

		fprintf( flog,
		"Can't find [%s].\n"
		"(Did you run gathermons.sht yet?)\n", mondir );

		exit( 42 );
	}

	fflush( flog );
}

/* --------------------------------------------------------------- */
/* GetZList ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetZList( vector<int> &zlist )
{
	for( int z = gArgs.zmin; z <= gArgs.zmax; ++z ) {

		char	path[2048];

		sprintf( path, "%s/X_A_%d.bin", gArgs.mondir, z );

		if( DskExists( path ) )
			zlist.push_back( z );
	}
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
// create top subdir
	DskCreateDir( gtopdir, flog );
}

/* --------------------------------------------------------------- */
/* WriteSubscapes ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubscapes( vector<int> &zlist )
{
// compose common argument string for all but last layer

	char	sopt[2048];

	sprintf( sopt,
	"'%s' -idb=%s"
	" -mb -scl=%d -lgord=%d -sdev=%d%s"
	" -ab -abwide=%d -stpcorr=%g",
	gArgs.mondir, idb.c_str(),
	gArgs.scl, gArgs.lgord, gArgs.sdev,
	(gArgs.resmask ? " -resmask" : ""),
	gArgs.abwide, gArgs.stpcorr );

// open file

	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/subscapes.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

// header

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# First step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > scapeops inpath -idb=idbpath -zb=%%d [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Scapeops does montage drawing and/or strip aligning as follows:\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# If drawing a montage...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t-mb -zb=%%d -scl=%%d\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t\t[-lgord=%%d] [-sdev=%%d] [-resmask]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# If aligning strips...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t-ab -za=%%d -zb=%%d -scl=%%d -abwide=%%d -stpcorr=%%lf\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t\t[-lgord=%%d] [-sdev=%%d] [-resmask] [-abdbg] [-abctr=%%lf]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -mb\t\t\t;make montage for layer zb\n" );
	fprintf( f, "# -ab\t\t\t;align layer za to zb\n" );
	fprintf( f, "# -za\t\t\t;layer za used only with -ab option\n" );
	fprintf( f, "# -zb\t\t\t;required layer zb\n" );
	fprintf( f, "# -scl=50\t\t;integer scale reduction\n" );
	fprintf( f, "# -lgord=1\t\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -sdev=42\t\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -resmask\t\t;mask out resin\n" );
	fprintf( f, "# -abwide=15\t;strips this many tiles wide on short axis\n" );
	fprintf( f, "# -stpcorr=0.02\t;required min corr for alignment\n" );
	fprintf( f, "# -abdbg\t\t;make diagnostic images and exit\n" );
	fprintf( f, "# -abctr=0\t\t;debug at this a-to-b angle\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Create output subdirs\n" );
	fprintf( f, "mkdir -p strips\n" );
	fprintf( f, "mkdir -p montages\n" );
	fprintf( f, "mkdir -p scplogs\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Submit layer pairs\n" );

// write all but last layer

	int	nz = zlist.size();

	for( int iz = 1; iz < nz; ++iz ) {

		fprintf( f,
		"qsub -N sc-%d -j y -o out.txt -b y -cwd -V -pe batch 8"
		" scapeops %s -za=%d -zb=%d\n",
		zlist[iz - 1],
		sopt, zlist[iz], zlist[iz - 1] );
	}

// last layer

	fprintf( f, "\n" );
	fprintf( f, "# Just montage last layer\n" );

	fprintf( f,
	"qsub -N sc-%d -j y -o out.txt -b y -cwd -V -pe batch 8"
	" scapeops '%s' -idb=%s -mb -scl=%d -lgord=%d -sdev=%d%s"
	" -zb=%d\n",
	zlist[nz - 1],
	gArgs.mondir, idb.c_str(),
	gArgs.scl, gArgs.lgord, gArgs.sdev,
	(gArgs.resmask ? " -resmask" : ""),
	zlist[nz - 1] );

	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteLowresgo ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteLowresgo()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/lowresgo.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Second step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Collects output from scapeops into 'LowRes.xml' stack\n" );
	fprintf( f, "# having one reduced scale montage per layer.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# You must subsequently view and edit LowRes.xml using TrakEM2.\n" );
	fprintf( f, "# Correct mistakes using 'Align with Manual Landmarks' feature\n" );
	fprintf( f, "# and then Save result with the same name.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_lowres -zmin=i -zmax=j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -table\t\t;write debugging striptable.txt and exit\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_lowres -zmin=%d -zmax=%d\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteScafgo --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScafgo()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/scafgo.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Third step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Apply the coarse layer-layer tforms from 'LowRes.xml' to the\n" );
	fprintf( f, "# individual tiles in 'mytemp/X_A_BIN_mons'. The new coarsely\n" );
	fprintf( f, "# aligned stack is named 'cross_wkspc/X_A_BIN_scaf'. It serves\n" );
	fprintf( f, "# both as the input for the block-block alignment, and as the\n" );
	fprintf( f, "# starting guess (scaffold) for the final LSQ alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_scaffold -zmin=i -zmax=j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -mons=../X_A_BIN_mons\t\t;alternate target montages\n" );
	fprintf( f, "# -lowres=LowRes.xml\t\t;alternate coarse alignment reference\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_scaffold -zmin=%d -zmax=%d\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteCarvego -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCarvego()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/carvego.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Fourth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Carve each layer into blocks of size bxb tiles and create new script\n" );
	fprintf( f, "# 'subblocks.sht' to distribute block-block alignment jobs to cluster.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_carveblocks -zmin=i -zmax=j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -b=10\t\t\t\t;ea. block roughly bXb tiles in area\n" );
	fprintf( f, "# -nf\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -scl=50\t\t\t;integer scale reduction\n" );
	fprintf( f, "# -lgord=1\t\t\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -sdev=42\t\t\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -resmask\t\t\t;mask out resin\n" );
	fprintf( f, "# -blkmincorr=0.45\t;required min corr for alignment\n" );
	fprintf( f, "# -blknomcorr=0.50\t;nominal corr for alignment\n" );
	fprintf( f, "# -xyconf=0.75\t\t;search radius = (1-conf)(blockwide)\n" );
	fprintf( f, "# -maxDZ=10\t\t\t;layers still correlate at this z-index span\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_carveblocks -zmin=%d -zmax=%d%s -b=10 -scl=%d -lgord=%d -sdev=%d%s -blkmincorr=%g -blknomcorr=%g -xyconf=%g\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""),
	gArgs.scl, gArgs.lgord, gArgs.sdev,
	(gArgs.resmask ? " -resmask" : ""),
	gArgs.blkmincorr, gArgs.blknomcorr, gArgs.xyconf );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<int>	zlist;

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

	IDBFromTemp( idb, ".", flog );

	if( idb.empty() )
		exit( 42 );

/* ---------------- */
/* Read source data */
/* ---------------- */

	GetZList( zlist );

	if( zlist.size() < 2 ) {
		fprintf( flog, "Fewer than 2 layers -- do nothing.\n" );
		goto exit;
	}

/* -------------- */
/* Create content */
/* -------------- */

	CreateTopDir();

	WriteSubscapes( zlist );
	WriteLowresgo();
	WriteScafgo();
	WriteCarvego();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


