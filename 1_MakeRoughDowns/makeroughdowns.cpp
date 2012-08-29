//
// Write script to submit roughdownpair jobs.
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
/* CArgs_rough --------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_rough {

public:
	double	Rstrip;
	char	xmlfile[2048],
			outdir[2048];
	char	*pat;
	int		zmin,
			zmax,
			nstriptiles,
			sclfac,
			sdnorm;

public:
	CArgs_rough()
	{
		Rstrip		= 0.14;
		xmlfile[0]	= 0;
		pat			= "/N";
		zmin		= 0;
		zmax		= 32768;
		nstriptiles	= 5;
		sclfac		= 200;
		sdnorm		= 0;	// 12 useful for Davi EM

		strcpy( outdir, "NoSuch" ); // protect real dirs
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_rough	gArgs;
static char			gtopdir[2048];
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_rough::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "makeroughdowns.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make down scripts: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: makeroughdowns <xmlfile> -dtemp -zmin=i -zmax=j"
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
		else if( GetArg( &nstriptiles, "-nstriptiles=%d", argv[i] ) )
			;
		else if( GetArg( &sclfac, "-scale=%d", argv[i] ) )
			;
		else if( GetArg( &sdnorm, "-sdnorm=%d", argv[i] ) )
			;
		else if( GetArg( &Rstrip, "-Rstrip=%lf", argv[i] ) )
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
	sprintf( gtopdir, "%s/roughdowns", gArgs.outdir );
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
/* WriteScript --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScript( vector<int> &zlist )
{
// compose common argument string

	char	sopt[2048];

	sprintf( sopt,
	"-d%s -p%s -nstriptiles=%d -scale=%d -sdnorm=%d -Rstrip=%g",
	gArgs.outdir, gArgs.pat, gArgs.nstriptiles,
	gArgs.sclfac, gArgs.sdnorm, gArgs.Rstrip );

// open file

	char	path[2048];
	FILE	*f;

	sprintf( path, "%s/subrough.sht", gtopdir );
	f = FileOpenOrDie( path, "w", flog );

// write

	int		nz = zlist.size();

	fprintf( f, "#!/bin/sh\n\n" );

	for( int iz = 1; iz < nz; ++iz ) {

		fprintf( f,
		"qsub -N rd-%d -j y -o out.txt -b y -cwd -V -pe batch 8"
		" roughdownpair '%s' %s -za=%d -zb=%d\n",
		zlist[iz], gArgs.xmlfile, sopt, zlist[iz], zlist[iz - 1] );
	}

	fprintf( f, "\n" );

	fclose( f );

// make executable

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

	WriteScript( zlist );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


