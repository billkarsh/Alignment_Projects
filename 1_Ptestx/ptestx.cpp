//
// ptestx is a wrapper function for ptest that first creates a disk
// folder/file structure that ptest expects, and then calls ptest.
//
// E.g.
// > ptestx za/ia@zb/ib -d=temp -prm=matchparams.txt
//	[-idb=xxx]
//	[-ima=xxx -imb=xxx]
//	[-fma=xxx -fmb=xxx]
//
// Required params:
//	za/ia@zb/ib:	z/tileid numerical labels (A onto B)
//	d:				output directory to create
//	prm:			path to parameter file
//
// Options:
//	idb:			path to IDB
//	ima, imb:		paths to the images
//	fma, fmb:		paths to foldmasks
//
// Image specs:
// If any of {ima, imb, fma, fmb} are specified, they override
// any idb-derived specs. However, idb is required if any image
// is to be fetched from a given IDB.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_ptx ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_ptx {

public:
	vector<char*>	passon;
	char			im[2][2048],
					fm[2][2048],
					tmp[2048],
					idb[2048],
					prm[2048];
	int				z[2], tile[2],
					clrtmp;

public:
	CArgs_ptx()
	{
		im[0][0]	= 0;
		im[1][0]	= 0;
		fm[0][0]	= 0;
		fm[1][0]	= 0;
		tmp[0]		= 0;
		idb[0]		= 0;
		prm[0]		= 0;
		clrtmp		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
	void TopLevel();
	void ALayer();
	void Call_ptest();
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_ptx	gArgs;
static FILE*		flog	= NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_ptx::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "ptestx.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 5 ) {
usage:
		printf( "Usage: ptestx <za/ia@zb/ib> -d=temp"
		" -prm=matchparams.txt [-idb=xxx]"
		" [-ima=xxx -imb=xxx -fma=zzz -fmb=xxx]"
		" [options].\n" );
		exit( 42 );
	}

	char	*key, *_arg;
	int		ok = true;

	passon.push_back( argv[1] );	// key

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			key = argv[i];
		else if( GetArgStr( _arg, "-d=", argv[i] ) )
			ok = DskAbsPath( tmp, sizeof(tmp), _arg, flog );
		else if( GetArgStr( _arg, "-idb=", argv[i] ) )
			ok = DskAbsPath( idb, sizeof(idb), _arg, flog );
		else if( GetArgStr( _arg, "-prm=", argv[i] ) )
			ok = DskAbsPath( prm, sizeof(prm), _arg, flog );
		else if( GetArgStr( _arg, "-ima=", argv[i] ) )
			ok = DskAbsPath( im[0], sizeof(im[0]), _arg, flog );
		else if( GetArgStr( _arg, "-imb=", argv[i] ) )
			ok = DskAbsPath( im[1], sizeof(im[1]), _arg, flog );
		else if( GetArgStr( _arg, "-fma=", argv[i] ) )
			ok = DskAbsPath( fm[0], sizeof(fm[0]), _arg, flog );
		else if( GetArgStr( _arg, "-fmb=", argv[i] ) )
			ok = DskAbsPath( fm[1], sizeof(fm[1]), _arg, flog );
		else if( IsArg( "-clr", argv[i] ) )
			clrtmp = true;
		else
			passon.push_back( argv[i] );

		if( !ok ) {
			fprintf( flog, "\nBad arg path '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );

	if( !key || !tmp[0] || !prm[0] )
		goto usage;

// Decode labels in key

	if( 4 != sscanf( key, "%d/%d@%d/%d",
		&z[0], &tile[0], &z[1], &tile[1] ) ) {

		fprintf( flog, "Bad label string '%s'.\n", key );
		exit( 42 );
	}

	fflush( flog );
}

/* --------------------------------------------------------------- */
/* TopLevel ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CArgs_ptx::TopLevel()
{
	char	buf[2048];

// preclear
	if( clrtmp ) {
		sprintf( buf, "rm -rf %s", tmp );
		system( buf );
	}

// create top directory
	DskCreateDir( tmp, flog );

// copy param files
	if( idb[0] ) {
		sprintf( buf, "cp %s/imageparams.txt %s", idb, tmp );
		system( buf );
	}

	sprintf( buf, "cp %s %s/matchparams.txt", prm, tmp );
	system( buf );
}

/* --------------------------------------------------------------- */
/* ALayer -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_ptx::ALayer()
{
	char	buf[2048];
	int		len;

/* ---------------------------- */
/* Create A layer and jobs dirs */
/* ---------------------------- */

	len = sprintf( buf, "%s/%d", tmp, z[0] );
	DskCreateDir( buf, flog );

	CreateJobsDir( buf, 0, 0, z[0], z[1], flog );

/* --------------- */
/* Create tile dir */
/* --------------- */

	sprintf( buf + len, "/%d", tile[0] );
	DskCreateDir( buf, flog );
}

/* --------------------------------------------------------------- */
/* Call_ptest ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_ptx::Call_ptest()
{
	char	buf[2048];
	int		len, narg;

// set jobs dir
	sprintf( buf, "%s/%d/%c0_0", tmp, z[0],
	(z[1] == z[0] ? 'S' : 'D') );

	chdir( buf );

// build command line with any passon args
	len  = sprintf( buf, "ptest" );
	narg = passon.size();

	for( int i = 0; i < narg; ++i )
		len += sprintf( buf + len, " %s", passon[i] );

// add optional image args
	if( im[0][0] )
		len += sprintf( buf + len, " -ima=%s", im[0] );

	if( im[1][0] )
		len += sprintf( buf + len, " -imb=%s", im[1] );

	if( fm[0][0] )
		len += sprintf( buf + len, " -fma=%s", fm[0] );

	if( fm[1][0] )
		len += sprintf( buf + len, " -fmb=%s", fm[1] );

// call ptest
	fprintf( flog, "%s\n.", buf );
	system( buf );
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

/* --------- */
/* Top level */
/* --------- */

	gArgs.TopLevel();

/* ------- */
/* A layer */
/* ------- */

	gArgs.ALayer();

/* --------------------- */
/* Call through to ptest */
/* --------------------- */

	gArgs.Call_ptest();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


