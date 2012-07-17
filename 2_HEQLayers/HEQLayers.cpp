//
// Write script to submit HEQ1Lyr jobs for:
//
// tag
// z range
// pct,
// lrbt
//


#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

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
	char	*infile,
			*tag;
	int		zmin,
			zmax;

public:
	CArgs_gray()
	{
		roi.L	= roi.R = 0;
		pct		= 99.5;
		infile	= NULL;
		tag		= NULL;
		zmin	= 0;
		zmax	= 32768;
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

	flog = FileOpenOrDie( "HEQLayers.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: HEQLayers <xml-file> <tag> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		vector<int>	vi;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile )
				infile = argv[i];
			else
				tag = argv[i];
		}
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
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
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
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
/* WriteScript --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScript( vector<int> &zlist )
{
// compose common argument string

	char	sopt[256];
	int		pos = 0;

	pos += sprintf( sopt + pos, "-pct=%g ", gArgs.pct );

	if( gArgs.roi.L != gArgs.roi.R ) {

		pos += sprintf( sopt + pos, "-lrbt=%d,%d,%d,%d ",
				gArgs.roi.L, gArgs.roi.R,
				gArgs.roi.B, gArgs.roi.T );
	}

// open file

	FILE	*f = FileOpenOrDie( "make.heq.sh", "w", flog );

// write

	int		nz = zlist.size();

	fprintf( f, "#!/bin/sh\n\n" );

	for( int iz = 0; iz < nz; ++iz ) {

		fprintf( f,
		"qsub -N heq-%d -j y -o out.txt -b y -cwd -V -pe batch 4"
		" HEQ1Lyr '%s' %s -z=%d %s\n",
		zlist[iz], gArgs.infile, gArgs.tag, zlist[iz], sopt );
	}

	fprintf( f, "\n" );

	fclose( f );

// make executable

	chmod( "make.heq.sh", S_IRWXU | S_IRWXG | S_IRWXO );
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

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


