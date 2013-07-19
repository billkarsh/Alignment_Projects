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
	double	abcorr,
			xyconf;		// neib radius = (1-conf)(blockwide)
	char	xmlfile[2048],
			outdir[2048];
	int		zmin,
			zmax,
			abwide,
			abscl,
			ablgord,
			absdev;
	bool	NoFolds;

public:
	CArgs_cross()
	{
		abcorr		= 0.20;
		xyconf		= 0.50;
		xmlfile[0]	= 0;
		zmin		= 0;
		zmax		= 32768;
		abwide		= 5;
		abscl		= 200;
		ablgord		= 1;	// 3  probably good for Davi EM
		absdev		= 0;	// 42 useful for Davi EM
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

		char	*_outdir;

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
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &xyconf, "-xyconf=%lf", argv[i] ) ) {

			if( xyconf < 0.0 || xyconf > 1.0 )
				xyconf = 0.5;
		}
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
	" -ab -abwide=%d -abscl=%d -ablgord=%d -absdev=%d -abcorr=%g",
	gArgs.xmlfile,
	gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.abwide, gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.abcorr );

// open file

	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/subscapes.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

// subdirs

	fprintf( f, "mkdir -p strips\n" );
	fprintf( f, "mkdir -p montages\n" );
	fprintf( f, "mkdir -p scplogs\n\n" );

// write all but last layer

	int	nz = zlist.size();

	for( int iz = 1; iz < nz; ++iz ) {

		fprintf( f,
		"qsub -N rd-%d -j y -o out.txt -b y -cwd -V -pe batch 8"
		" scapeops %s -za=%d -zb=%d\n",
		zlist[iz - 1],
		sopt, zlist[iz], zlist[iz - 1] );
	}

// last layer

	fprintf( f,
	"qsub -N rd-%d -j y -o out.txt -b y -cwd -V -pe batch 8"
	" scapeops '%s' -mb -mbscl=%d -mblgord=%d -mbsdev=%d"
	" -zb=%d\n",
	zlist[nz - 1],
	gArgs.xmlfile, gArgs.abscl,
	gArgs.ablgord, gArgs.absdev, zlist[nz - 1] );

	fprintf( f, "\n" );

	fclose( f );

// make executable

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

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f,
	"cross_lowres -zmin=%d -zmax=%d\n\n",
	gArgs.zmin, gArgs.zmax );

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

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f,
	"cross_lowtohires '%s' -lowres=LowRes.xml"
	" -zmin=%d -zmax=%d"
	" -xmltype=0 -xmlmin=0 -xmlmax=0\n\n",
	gArgs.xmlfile, gArgs.zmin, gArgs.zmax );

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

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f,
	"cross_carveblocks"
	" HiRes.xml -zmin=%d -zmax=%d%s"
	" -b=10 -abscl=%d -ablgord=%d -absdev=%d"
	" -abcorr=%g -xyconf=%g\n\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""),
	gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.abcorr, gArgs.xyconf );

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
	WriteCarvego();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


