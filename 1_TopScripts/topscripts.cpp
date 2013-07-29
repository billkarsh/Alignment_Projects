//
// Create the top-level scripts for alignment.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"


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
		"Usage: topscripts <infile> -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
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
	fprintf( f, "# Make image database 'idb' from 'myxml' file...\n" );
	fprintf( f, "# Makes everything.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > makeidb myxml -d=idb\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -zmin=i -zmax=j\t\t\t;restricts layer range\n" );
	fprintf( f, "# -nf\t\t\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -lens=AffineLensFile.txt\t;external affine software lens\n" );
	fprintf( f, "# -k=MyClicksfile.txt\t\t;align with manual landmarks\n" );
	fprintf( f, "# -xmltype=0\t\t\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "makeidb %s -d=idb0 -zmin=%d -zmax=%d%s%s\n",
	gArgs.infile, gArgs.zmin, gArgs.zmax,
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
	fprintf( f, "# > makemontages myidb -d=temp -zmin=0 -zmax=5\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -nf\t\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -nd\t\t\t\t\t;no tile directories\n" );
	fprintf( f, "# -b=8\t\t\t\t\t;subblock 1D size <in tiles>\n" );
	fprintf( f, "# -minareafrac=0.02\t\t;same and cross req. area olap frac\n" );
	fprintf( f, "# -exe=ptestalt\t\t\t;exe other than 'ptest'\n" );
	fprintf( f, "# -lens\t\t\t\t\t;use lens file\n" );
	fprintf( f, "# -davinc\t\t\t\t;no davi bock same lyr corners\n" );
	fprintf( f, "# -xmltype=0\t\t\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "wrk=temp0\n" );
	fprintf( f, "\n" );
	fprintf( f, "makemontages idb0 -d=$wrk -zmin=%d -zmax=%d%s -nd -b=8%s\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""), xmlprms );
	fprintf( f, "\n" );
	fprintf( f, "cp prms/* $wrk\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* Write_upmongo --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_upmongo()
{
	char	buf[2048], inname[256];
	FILE	*f;

	if( FileExtIsXML( gArgs.infile ) ) {

		const char	*name	= FileNamePtr( gArgs.infile ),
					*dot	= FileDotPtr( name );

		sprintf( inname, "%.*s", dot - name, name );
	}
	else
		strcpy( inname, "PreClicks" );

	sprintf( buf, "upmongo.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Update transforms in given 'myxml' file using TAffineTables in given\n" );
	fprintf( f, "# temp directory. Preserves montage-montage orientation. Output file\n" );
	fprintf( f, "# named 'myxml_v2.xml'.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > updatemontages myxml temp\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -zmin=i -zmax=j\t;restricts layer range\n" );
	fprintf( f, "# -force\t\t\t;overwrite TForms with lsq solutions\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "inname=%s\n",
	inname );
	fprintf( f, "\n" );
	fprintf( f, "updatemontages $inname.xml temp0 -zmin=%d -zmax=%d\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );
	fprintf( f, "mv $inname\"_v2.xml\" \"newmons.xml\"\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* Write_crossgo ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_crossgo()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "crossgo.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Write scripts governing cross layer alignment.\n" );
	fprintf( f, "# Creates folder mytemp/cross_wkspc\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_topscripts myxml -d=temp0 -zmin=0 -zmax=5\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -nf\t\t\t;no foldmasks\n" );
	fprintf( f, "# -abwide=5\t\t;strip width in tiles\n" );
	fprintf( f, "# -abscl=200\t;size reduction factor\n" );
	fprintf( f, "# -ablgord=1\t;Legendre poly max order\n" );
	fprintf( f, "# -absdev=42\t;scape sdev size\n" );
	fprintf( f, "# -abcorr=0.20\t;req. min corr\n" );
	fprintf( f, "# -xyconf=0.5\t;search radius = (1-conf)(blockwide)\n" );
	fprintf( f, "# -xmltype=0\t;ImagePlus type code\n" );
	fprintf( f, "# -xmlmin=0\t\t;intensity scale\n" );
	fprintf( f, "# -xmlmax=0\t\t;intensity scale\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "cross_topscripts newmons.xml -d=temp0 -zmin=%d -zmax=%d%s -abwide=5 -abscl=200 -ablgord=1 -absdev=42 -abcorr=.20 -xyconf=0.5%s\n",
	gArgs.zmin, gArgs.zmax,
	(gArgs.NoFolds ? " -nf" : ""), xmlprms );
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

	int	L = 0;

	if( gArgs.xml_type != 0 )
		L += sprintf( xmlprms + L, " -xmltype=%d", gArgs.xml_type );

	if( gArgs.xml_min != 0 )
		L += sprintf( xmlprms + L, " -xmlmin=%d", gArgs.xml_min );

	if( gArgs.xml_max != 0 )
		L += sprintf( xmlprms + L, " -xmlmax=%d", gArgs.xml_max );

	Write_dbgo();
	Write_mongo();
	Write_upmongo();
	Write_crossgo();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


