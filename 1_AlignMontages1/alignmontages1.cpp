//
// Write script to submit scapeops jobs (subscapes).
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
/* CArgs_alnmon -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_alnmon {

public:
	double	abcorr;
	char	xmlfile[2048],
			outdir[2048];
	char	*pat;
	int		zmin,
			zmax,
			abwide,
			abscl,
			absdev;

public:
	CArgs_alnmon()
	{
		abcorr		= 0.20;
		xmlfile[0]	= 0;
		pat			= "/N";
		zmin		= 0;
		zmax		= 32768;
		abwide		= 5;
		abscl		= 200;
		absdev		= 0;	// 12 useful for Davi EM

		strcpy( outdir, "NoSuch" ); // protect real dirs
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_alnmon	gArgs;
static char			gtopdir[2048];
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_alnmon::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "alignmontages1.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make scapeops scripts: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: alignmontages1 <xmlfile> -dtemp -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		char	*_outdir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			DskAbsPath( xmlfile, sizeof(xmlfile), argv[i], flog );
		else if( GetArgStr( _outdir, "-d", argv[i] ) )
			DskAbsPath( outdir, sizeof(outdir), _outdir, flog );
		else if( GetArgStr( pat, "-p", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &abwide, "-abwide=%d", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
			;
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
	sprintf( gtopdir, "%s/alignmontages", gArgs.outdir );
	DskCreateDir( gtopdir, flog );
}

/* --------------------------------------------------------------- */
/* ScriptPerms --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScriptPerms( const char *path )
{
	char	buf[2048];

	sprintf( buf, "chmod ug=rwx,o=rx %s", path );
	system( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubscapes ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubscapes( vector<int> &zlist )
{
// compose common argument string for all but last layer

	char	sopt[2048];

	sprintf( sopt,
	"'%s' -d%s -p%s"
	" -mb -mbscl=%d"
	" -ab -abwide=%d -abscl=%d -absdev=%d -abcorr=%g",
	gArgs.xmlfile, gArgs.outdir, gArgs.pat,
	gArgs.abscl,
	gArgs.abwide, gArgs.abscl, gArgs.absdev, gArgs.abcorr );

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
	" scapeops '%s' -d%s -p%s -mb -mbscl=%d -zb=%d\n",
	zlist[nz - 1],
	gArgs.xmlfile, gArgs.outdir, gArgs.pat,
	gArgs.abscl, zlist[nz - 1] );

	fprintf( f, "\n" );

	fclose( f );

// make executable

	ScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteGather --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteGather()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/gather.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f,
	"alignmontages2 -d. -zmin=%d -zmax=%d\n\n",
	gArgs.zmin, gArgs.zmax );

	fclose( f );
	ScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteFinal ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFinal()
{
	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/final.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f,
	"alignmontages3 '%s' -d. -p%s -zmin=%d -zmax=%d"
	" -xmltype=0 -xmlmin=0 -xmlmax=0\n\n",
	gArgs.xmlfile, gArgs.pat, gArgs.zmin, gArgs.zmax );

	fclose( f );
	ScriptPerms( path );
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
	WriteGather();
	WriteFinal();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


