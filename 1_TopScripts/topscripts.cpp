//
// Create the top-level scripts for alignment.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"

#include	<string.h>
#include	<time.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_scr ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scr {

public:
	char	script[2048],
			idb[2048];
	char	*infile;
	int		zmin,
			zmax;

public:
	CArgs_scr()
	{
		script[0]	= 0;
		idb[0]		= 0;
		infile		= NULL;
		zmin		= 0;
		zmax		= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scr	gArgs;
static FILE*		flog	= NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "topscripts.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Create top-level scripts: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: topscripts <infile> -script=scriptpath"
		" -idb=idbpath -z=i,j.\n" );
		exit( 42 );
	}

	vector<int>	vi;
	const char	*pchar;

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( pchar, "-script=", argv[i] ) )
			DskAbsPath( script, sizeof(script), pchar, flog );
		else if( GetArgStr( pchar, "-idb=", argv[i] ) )
			DskAbsPath( idb, sizeof(idb), pchar, flog );
		else if( GetArgList( vi, "-z=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zmin = vi[0];
				zmax = vi[1];
			}
			else {
				fprintf( flog,
				"Bad format in -z [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Write_dbgo ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_dbgo()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "dbgo.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Make image database 'idb' from text 'billfile' or xml file...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > makeidb inpath -script=scriptpath -idb=idbpath -z=i,j\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Required:\n" );
	fprintf( f, "# inPath\t\t\t\t;mylayout.txt or myxml.xml file\n" );
	fprintf( f, "# -script=scriptpath\t;alignment pipeline params file\n" );
	fprintf( f, "# -idb=idbpath\t\t\t;idb folder to create\n" );
	fprintf( f, "# -z=i,j\t\t\t\t;process layers in range z=[i..j]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# NOTE: Omit -idb parameter to generate RawData.xml.\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "makeidb %s -script=%s -idb=%s -z=%d,%d\n",
	gArgs.infile, gArgs.script, gArgs.idb,
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* Write_mongo --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_mongo()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "mongo.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Make alignment workspace 'temp' and create scripts within\n" );
	fprintf( f, "# to drive montaging pipeline steps.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > makemontages temp -script=scriptpath -idb=idbpath -z=i,j\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Required:\n" );
	fprintf( f, "# temp\t\t\t\t\t;workspace directory to create\n" );
	fprintf( f, "# -script=scriptpath\t;alignment pipeline params file\n" );
	fprintf( f, "# -idb=idbpath\t\t\t;path to idb directory\n" );
	fprintf( f, "# -z=i,j\t\t\t\t;align layers in range z=[i..j]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -exe=ptestalt\t\t\t;exe other than 'ptest'\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "wrk=temp0\n" );
	fprintf( f, "\n" );
	fprintf( f, "makemontages $wrk -script=%s -idb=%s -z=%d,%d\n",
	gArgs.script, gArgs.idb,
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );
	fprintf( f, "cp prms/* $wrk\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
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

/* ------- */
/* Scripts */
/* ------- */

	Write_dbgo();
	Write_mongo();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


