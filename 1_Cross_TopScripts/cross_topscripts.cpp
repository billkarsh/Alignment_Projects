//
// Write scripts governing cross layer alignment.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
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
	char	xmlfile[2048],
			outdir[2048];
	int		zmin,
			zmax,
			abwide,
			abscl,
			ablgord,
			absdev,
			xml_type,
			xml_min,
			xml_max;
	bool	NoFolds;

public:
	CArgs_cross()
	{
		stpcorr		= 0.02;
		blkmincorr	= 0.45;
		blknomcorr	= 0.50;
		xyconf		= 0.75;
		xmlfile[0]	= 0;
		zmin		= 0;
		zmax		= 32768;
		abwide		= 15;
		abscl		= 50;
		ablgord		= 1;	// 1  probably good for Davi EM
		absdev		= 42;	// 42 useful for Davi EM
		xml_type	= 0;
		xml_min		= 0;
		xml_max		= 0;
		NoFolds		= false;

		strcpy( outdir, "NoSuch" ); // protect real dirs
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_cross	gArgs;
static char			gtopdir[2048];
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

	if( argc < 5 ) {
		printf(
		"Usage: cross_topscripts <xmlfile> -d=temp -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		const char	*_outdir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			DskAbsPath( xmlfile, sizeof(xmlfile), argv[i], flog );
		else if( GetArgStr( _outdir, "-d=", argv[i] ) )
			DskAbsPath( outdir, sizeof(outdir), _outdir, flog );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &abwide, "-abwide=%d", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			;
		else if( GetArg( &ablgord, "-ablgord=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
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
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			;
		else if( GetArg( &xml_min, "-xmlmin=%d", argv[i] ) )
			;
		else if( GetArg( &xml_max, "-xmlmax=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<int> &zlist )
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.xmlfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		zlist.push_back( z );
	}
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
// create top subdir
	sprintf( gtopdir, "%s/cross_wkspc", gArgs.outdir );
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
	"'%s'"
	" -mb -mbscl=%d -mblgord=%d -mbsdev=%d"
	" -ab -abwide=%d -abscl=%d -ablgord=%d -absdev=%d -stpcorr=%g",
	gArgs.xmlfile,
	gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.abwide, gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.stpcorr );

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
	fprintf( f, "# > scapeops myxml -zb=%%d [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Scapeops does montage drawing and/or strip aligning as follows:\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# If drawing a montage...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t-mb -zb=%%d -mbscl=%%d\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t\t[-mblgord=%%d] [-mbsdev=%%d]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# If aligning strips...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t-ab -za=%%d -zb=%%d -abwide=%%d -abscl=%%d -stpcorr=%%lf\n" );
	fprintf( f, "#\n" );
	fprintf( f, "#\t\t[-ablgord=%%d] [-absdev=%%d] [-abdbg] [-abctr=%%lf]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -mb\t\t\t;make montage for layer zb\n" );
	fprintf( f, "# -ab\t\t\t;align layer za to zb\n" );
	fprintf( f, "# -za\t\t\t;layer za used only with -ab option\n" );
	fprintf( f, "# -zb\t\t\t;required layer zb\n" );
	fprintf( f, "# -mbscl=50\t\t;integer scale reduction\n" );
	fprintf( f, "# -abscl=50\t\t;integer scale reduction\n" );
	fprintf( f, "# -mblgord=1\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -ablgord=1\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -mbsdev=42\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -absdev=42\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -abwide=15\t;strips this many tiles wide on short axis\n" );
	fprintf( f, "# -stpcorr=0.02\t;required min corr for alignment\n" );
	fprintf( f, "# -addbg\t\t;make diagnostic images and exit\n" );
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
	" scapeops '%s' -mb -mbscl=%d -mblgord=%d -mbsdev=%d"
	" -zb=%d\n",
	zlist[nz - 1],
	gArgs.xmlfile, gArgs.abscl,
	gArgs.ablgord, gArgs.absdev, zlist[nz - 1] );

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
	fprintf( f, "# > cross_lowres -zmin=i -zmax=j\n" );
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
/* WriteHiresgo -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteHiresgo()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/hiresgo.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Third step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Distribute the coarse layer-layer transforms from 'LowRes.xml'\n" );
	fprintf( f, "# to the individual tiles in (usually) 'newmons.xml', the stack\n" );
	fprintf( f, "# of precision single-layer montaging results. The new stack is\n" );
	fprintf( f, "# named 'HiRes.xml'.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# The resulting coarsely aligned full resolution 'HiRes.xml'\n" );
	fprintf( f, "# can serve as a scaffold in the final LSQ alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_lowtohires montaged_fullres_xml -zmin=i -zmax=j\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -lowres=LowRes.xml\t;alternate coarse alignment reference\n" );
	fprintf( f, "# -xmltype=0\t\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_lowtohires '%s' -zmin=%d -zmax=%d -xmltype=%d -xmlmin=%d -xmlmax=%d\n",
	gArgs.xmlfile, gArgs.zmin, gArgs.zmax,
	gArgs.xml_type, gArgs.xml_min, gArgs.xml_max );
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
	fprintf( f, "# Fourth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Create folder 'X_A_TXT' of coarse affine transforms,\n" );
	fprintf( f, "# one per tile, to be used as a scaffold in final LSQ runs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > XMLGetTF HiRes.xml -zmin=i -zmax=j\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "XMLGetTF HiRes.xml -zmin=%d -zmax=%d\n",
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
	fprintf( f, "# Fifth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Carve each layer into blocks of size bxb tiles and create new script\n" );
	fprintf( f, "# 'subblocks.sht' to distribute block-block alignment jobs to cluster.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_carveblocks myxml -zmin=i -zmax=j [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -b=10\t\t\t\t;ea. block roughly bXb tiles in area\n" );
	fprintf( f, "# -nf\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -abscl=50\t\t\t;integer scale reduction\n" );
	fprintf( f, "# -ablgord=1\t\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -absdev=42\t\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -blkmincorr=0.45\t;required min corr for alignment\n" );
	fprintf( f, "# -blknomcorr=0.50\t;nominal corr for alignment\n" );
	fprintf( f, "# -xyconf=0.75\t\t;search radius = (1-conf)(blockwide)\n" );
	fprintf( f, "# -maxDZ=10\t\t\t;layers still correlate at this z-index span\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_carveblocks HiRes.xml -zmin=%d -zmax=%d%s -b=10 -abscl=%d -ablgord=%d -absdev=%d -blkmincorr=%g -blknomcorr=%g -xyconf=%g\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""),
	gArgs.abscl, gArgs.ablgord, gArgs.absdev,
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

/* ---------------- */
/* Read source file */
/* ---------------- */

	ParseTrakEM2( zlist );

	if( zlist.size() < 2 )
		goto exit;

/* -------------- */
/* Create content */
/* -------------- */

	CreateTopDir();

	WriteSubscapes( zlist );
	WriteLowresgo();
	WriteHiresgo();
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


