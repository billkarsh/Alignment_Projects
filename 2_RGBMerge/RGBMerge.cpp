//
// Write script to submit RGBM1Lyr jobs for:
//
// tag
// z range
// -R=chn,pct
// -G=chn,pct
// -B=chn,pct
// -spanRGB=LLT
// -lrbt
//
// -spanRGB option needs three character string. Each character is
// either {L=whole layer, T=per tile} setting span of tiles used
// to scale each channel.
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
/* CArgs_rgbm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_rgbm {

public:
	IBox		roi;
	double		pct[3];
	const char	*infile,
				*tag,
				*span;
	int			zmin, zmax,
				RGB[3];

public:
	CArgs_rgbm()
	{
		roi.L	= roi.R = 0;
		pct[0]	= 99.5;
		pct[1]	= 99.5;
		pct[2]	= 99.5;
		infile	= NULL;
		tag		= NULL;
		span	= "LLL";
		zmin	= 0;
		zmax	= 32768;
		RGB[0]	= -1;
		RGB[1]	= -1;
		RGB[2]	= -1;
	};

	bool ScanChan( int chn, const char *pat, char *argv );
	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_rgbm	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* ScanChan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Scan channel/pct arguments like:
// -R=0[,95.5]
//
bool CArgs_rgbm::ScanChan( int chn, const char *pat, char *argv )
{
	int	c;

	if( 1 == sscanf( argv, pat, &c ) ) {

		double	p;

		RGB[chn] = c;

		if( argv[4] == ',' && 1 == sscanf( argv + 5, "%lf", &p ) )
			pct[chn] = p;

		return true;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_rgbm::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "RGBMerge.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: RGBMerge <xml-file> <tag>"
		" <-[R,G,B]=i,pct> [options].\n" );
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
		else if( ScanChan( 0, "-R=%d", argv[i] ) )
			;
		else if( ScanChan( 1, "-G=%d", argv[i] ) )
			;
		else if( ScanChan( 2, "-B=%d", argv[i] ) )
			;
		else if( GetArgStr( span, "-spanRGB=", argv[i] ) )
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

	const char	*sRGB = "RGB";
	char		sopt[256];
	int			pos = 0;

	for( int i = 0; i < 3; ++i ) {

		if( gArgs.RGB[i] >= 0 ) {
			pos += sprintf( sopt + pos, "-%c=%d,%g ",
					sRGB[i], gArgs.RGB[i], gArgs.pct[i] );
		}
	}

	pos += sprintf( sopt + pos, "-spanRGB=%s ", gArgs.span );

	if( gArgs.roi.L != gArgs.roi.R ) {

		pos += sprintf( sopt + pos, "-lrbt=%d,%d,%d,%d ",
				gArgs.roi.L, gArgs.roi.R,
				gArgs.roi.B, gArgs.roi.T );
	}

// open file

	FILE	*f = FileOpenOrDie( "make.merge.sht", "w", flog );

// write

	int		nz = zlist.size();

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );

	for( int iz = 0; iz < nz; ++iz ) {

		fprintf( f,
		"QSUB_1NODE.sht \"rgbm-%d\" \"-j y -o out.txt\" 4"
		" \"RGBM1Lyr '%s' %s -z=%d %s\"\n",
		zlist[iz], gArgs.infile, gArgs.tag, zlist[iz], sopt );
	}

	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( "make.merge.sht" );
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


