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
	char	idb[2048];
	char	*infile;
	int		zmin,
			zmax,
			xml_type,
			xml_min,
			xml_max;
	bool	NoFolds;

public:
	CArgs_scr()
	{
		infile		= NULL;
		zmin		= 0;
		zmax		= 32768;
		xml_type	= 0;
		xml_min		= 0;
		xml_max		= 0;
		NoFolds		= false;

		strcpy( idb, "NoSuch" ); // protect real dirs
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scr	gArgs;
static FILE*		flog			= NULL;
static char			xmlprms[256]	= {0};






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
		"Usage: topscripts <infile> -idb=idbpath -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		const char	*_dir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( _dir, "-idb=", argv[i] ) )
			DskAbsPath( idb, sizeof(idb), _dir, flog );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			;
		else if( GetArg( &xml_min, "-xmlmin=%d", argv[i] ) )
			;
		else if( GetArg( &xml_max, "-xmlmax=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
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
	fprintf( f, "# Make image database 'idb' from Rick or xml file...\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > makeidb <mylayout.txt or myxml.xml> -idb=idb\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -zmin=i -zmax=j\t\t\t;restricts layer range\n" );
	fprintf( f, "# -nf\t\t\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -crop=CropRectsFile.txt\t;{cam,x0,y0,dx,dy} for ea. cam\n" );
	fprintf( f, "# -xmltype=0\t\t\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t\t\t;intensity scale\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# NOTE: Omit -idb parameter to generate RawData.xml.\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "makeidb %s -idb=%s -zmin=%d -zmax=%d%s%s\n",
	gArgs.infile, gArgs.idb, gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""), xmlprms );
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
	fprintf( f, "# Make workspace 'temp' for 'myidb'...\n" );
	fprintf( f, "# Makes everything except prms/matchparams.txt file.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > makemontages myidb -d=temp -zmin=i -zmax=j\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -nf\t\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -nd\t\t\t\t\t;no tile directories\n" );
	fprintf( f, "# -b=8\t\t\t\t\t;subblock 1D size <in tiles>\n" );
	fprintf( f, "# -minareafrac=0.02\t\t;same and cross req. area olap frac\n" );
	fprintf( f, "# -exe=ptestalt\t\t\t;exe other than 'ptest'\n" );
	fprintf( f, "# -davinc\t\t\t\t;no davi bock same lyr corners\n" );
	fprintf( f, "# -xmltype=0\t\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "wrk=temp0\n" );
	fprintf( f, "\n" );
	fprintf( f, "makemontages %s -d=$wrk -zmin=%d -zmax=%d%s -nd -b=8%s\n",
	gArgs.idb, gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""), xmlprms );
	fprintf( f, "\n" );
	fprintf( f, "cp prms/* $wrk\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* XMLHasMRC ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if *.mrc images.
//
static bool XMLHasMRC()
{
	FILE		*f = FileOpenOrDie( gArgs.infile, "r", flog );
	CLineScan	LS;
	int			nfile_path = 0;
	bool		ismrc = false;

	while( LS.Get( f ) > 0 ) {

		if( strstr( LS.line, "ex.mrc\"" ) ) {
			ismrc = true;
			break;
		}

		if( nfile_path >= 10 )
			break;

		if( strstr( LS.line, "file_path" ) ) {
			++nfile_path;
		}
	}

	fclose( f );

	return ismrc;
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

	int	L = 0;

	if( gArgs.xml_type != 0 )
		L += sprintf( xmlprms + L, " -xmltype=%d", gArgs.xml_type );

	if( gArgs.xml_min != 0 )
		L += sprintf( xmlprms + L, " -xmlmin=%d", gArgs.xml_min );

	if( gArgs.xml_max != 0 )
		L += sprintf( xmlprms + L, " -xmlmax=%d", gArgs.xml_max );

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


