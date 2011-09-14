//
// ptestx is a wrapper function for ptest that first creates a disk
// folder/file structure that ptest expects, and then calls ptest.
//
// E.g.
// > ptestx za/ia@zb/ib -d=temp
//	-thm=thmparams.txt -msh=dmeshparams.txt
//	[-idb=xxx]
//	[-ima=xxx -imb=xxx]
//	[-fma=xxx -fmb=xxx]
//
// Required params:
//	za/ia@zb/ib:	z/tileid numerical labels (A onto B)
//	d:				output directory to create
//	thm, msh:		paths to parameter files
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
					thm[2048],
					msh[2048];
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
		thm[0]		= 0;
		msh[0]		= 0;
		clrtmp		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
	void TopLevel();
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

	if( argc < 6 ) {
usage:
		printf( "Usage: ptestx <za/ia@zb/ib> -d=temp"
		" -thm=thmparams.txt -msh=dmeshparams.txt"
		" [-idb=xxx] [-ima=xxx -imb=xxx -fma=zzz -fmb=xxx]"
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
		else if( GetArgStr( _arg, "-thm=", argv[i] ) )
			ok = DskAbsPath( thm, sizeof(thm), _arg, flog );
		else if( GetArgStr( _arg, "-msh=", argv[i] ) )
			ok = DskAbsPath( msh, sizeof(msh), _arg, flog );
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

	if( !key || !tmp[0] || !thm[0] || !msh[0] )
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

	sprintf( buf, "cp %s %s", thm, tmp );
	system( buf );

	sprintf( buf, "cp %s %s", msh, tmp );
	system( buf );
}

/* --------------------------------------------------------------- */
/* ThmPair ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void ThmPair( const char *sz, int za, int zb )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/ThmPair_%d_@_%d.txt", sz, za, zb );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "Atl\tBtl\tAcr\tBcr\tErr\tDeg\tQ\tR"
	"\tT0\tT1\tX\tT3\tT4\tY\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* TileDir ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void TileDir( const char *sz, int itile )
{
	char	buf[4096];

// create tile dir
	sprintf( buf, "%s/%d", sz, gArgs.tile[itile] );
	DskCreateDir( buf, flog );
}

/* --------------------------------------------------------------- */
/* EachLayer ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void EachLayer()
{
	for( int iz = 0; iz < 2; ++iz ) {

		char	sz[2048], buf[2048];

		/* ---------------- */
		/* Different layer? */
		/* ---------------- */

		if( iz == 1 && (gArgs.z[1] == gArgs.z[0]) )
			break;

		/* ---------------- */
		/* Create layer dir */
		/* ---------------- */

		sprintf( sz, "%s/%d", gArgs.tmp, gArgs.z[iz] );
		DskCreateDir( sz, flog );

		/* ------- */
		/* ThmPair */
		/* ------- */

		ThmPair( sz, gArgs.z[iz], gArgs.z[iz] );

		if( gArgs.z[1] != gArgs.z[0] )
			ThmPair( sz, gArgs.z[iz], gArgs.z[!iz] );

		/* ----------------------------- */
		/* Create tile dir (1 or 2 dirs) */
		/* ----------------------------- */

		TileDir( sz, iz );

		if( gArgs.z[1] == gArgs.z[0] )
			TileDir( sz, 1 );
	}
}

/* --------------------------------------------------------------- */
/* Call_ptest ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_ptx::Call_ptest()
{
	char	buf[2048];
	int		len, narg;

// set working dir to layer A
	sprintf( buf, "%s/%d", tmp, z[0] );
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

/* ---------- */
/* Each layer */
/* ---------- */

	EachLayer();

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


