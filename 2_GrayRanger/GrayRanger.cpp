//
// Write script to submit GraRan1Lyr jobs for:
//
// z range
// chn,
// pct,
// lrbt
//


#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"

#include	"tinyxml.h"

#include	<sys/stat.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_gray ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_gray {

public:
	IBox	roi;
	double	pct;
	char	*infile;
	int		zmin,
			zmax,
			chn;

public:
	CArgs_gray()
	{
		roi.L	= roi.R = 0;
		pct		= 99.5;
		infile	= NULL;
		zmin	= 0;
		zmax	= 32768;
		chn		= -1;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_gray	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_gray::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "GrayRanger.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: GrayRanger <xml-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		vector<int>	vi;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &chn, "-chn=%d", argv[i] ) )
			;
		else if( GetArg( &pct, "-pct=%lf", argv[i] ) )
			;
		else if( GetArgList( vi, "-lrbt=", argv[i] ) && vi.size() == 4 )
			memcpy( &roi, &vi[0], 4*sizeof(int) );
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
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	layer = hDoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

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
/* WriteScript --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScript( vector<int> &zlist )
{
// compose common argument string

	char	sopt[256];
	int		pos = 0;

	if( gArgs.chn >= 0 )
		pos += sprintf( sopt + pos, "-chn=%d ", gArgs.chn );

	pos += sprintf( sopt + pos, "-pct=%g ", gArgs.pct );

	if( gArgs.roi.L != gArgs.roi.R ) {

		pos += sprintf( sopt + pos, "-lrbt=%d,%d,%d,%d ",
				gArgs.roi.L, gArgs.roi.R,
				gArgs.roi.B, gArgs.roi.T );
	}

// open file

	FILE	*f = FileOpenOrDie( "make.scales", "w", flog );

// write

	int		nz = zlist.size();

	fprintf( f, "#!/bin/csh\n\n" );

	for( int iz = 0; iz < nz; ++iz ) {

		fprintf( f,
		"qsub -N gr-%d -j y -o out.txt -b y -cwd -V -pe batch 4"
		" GraRan1Lyr '%s' -z=%d %s\n",
		zlist[iz], gArgs.infile, zlist[iz], sopt );
	}

	fprintf( f, "\n" );

	fclose( f );

// make executable

	chmod( "make.scales", S_IRWXU | S_IRWXG | S_IRWXO );
}

/* --------------------------------------------------------------- */
/* GRTemp -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GRTemp( vector<int> &zlist )
{
	DskCreateDir( "GRTemp", flog );

	FILE	*f = FileOpenOrDie( "GRTemp/range.txt", "w", flog );

	fprintf( f, "%d\t%d\n", zlist[0], zlist.size() );
	fclose( f );
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

	if( !zlist.size() )
		goto exit;

/* ------------ */
/* Write script */
/* ------------ */

	WriteScript( zlist );

/* ------------- */
/* Create GRTemp */
/* ------------- */

	GRTemp( zlist );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


