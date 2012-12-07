//
// "FF tomography images"
//
// Apply flatfield to all images in a rick file.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	char	*prmfile;

public:
	CArgs()
	{
		prmfile	= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static FILE*	flog = NULL;
static char		tifdir[2048],
				ffdir[2048];
static FILE*	frick = NULL;
static uint16*	ffras[4] = {NULL,NULL,NULL,NULL};
static double	ffave[4];
static uint32	gW = 0,	gH = 0;		// universal pic dims






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "FFTomos.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: fftomos <paramfile>.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			prmfile = argv[i];
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* ReadParams ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Param file:
// TIFpath=./TIF
// rick=fullpath
// chan=0,fullpath
// ... as many as used chans, up to 4
//
static void ReadParams()
{
	FILE*		f = FileOpenOrDie( gArgs.prmfile, "r" );
	CLineScan	LS;
	char		buf[2048];

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "TIFpath=%s", tifdir ) )
		exit( 42 );
	fprintf( flog, "TIFpath=%s\n", tifdir );

	strcpy( ffdir, tifdir );
	char	*s = strrchr( ffdir, '/' );
	strcpy( s + 1, "FF" );
	DskCreateDir( ffdir, flog );

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "rick=%s", buf ) )
		exit( 42 );
	fprintf( flog, "rick=%s\n", buf );
	frick = FileOpenOrDie( buf, "r" );

	while( LS.Get( f ) > 0 ) {

		int	chan, np;

		if( 2 != sscanf( LS.line, "chan=%d,%s", &chan, buf ) )
			break;
		fprintf( flog, "chan=%d,%s\n", chan, buf );
		ffras[chan] = Raster16FromTif16( buf, gW, gH, flog );

		np = gW * gH;

		ffave[chan] = 0;
		for( int i = 0; i < np; ++i ) {
			if( ffras[chan][i] <= 0 )
				ffras[chan][i] = 1;
			ffave[chan] += ffras[chan][i];
		}
		ffave[chan] /= np;
	}

	fprintf( flog, "\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ProcessRick --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ProcessRick()
{
	CLineScan	LS;
	int			chan, np = gW * gH;

	while( LS.Get( frick ) > 0 ) {

		char	path[2048], name[256];

		sscanf( LS.line, "%s", name );

		sprintf( path, "%s/%s", tifdir, name );
		uint16*	ras = Raster16FromTif16( path, gW, gH, flog );

		if( !ras )
			break;

		char	*s = strrchr( name, '_' );
		chan = atoi( s + 1 );

		if( chan < 0 || chan > 3 || !ffras[chan] ) {
			fprintf( flog, "Bad chan # [%s]\n", name );
			goto done_ff;
		}

		for( int i = 0; i < np; ++i )
			ras[i] = uint16(ras[i]*ffave[chan]/ffras[chan][i]);

done_ff:
		sprintf( path, "%s/%.*s.FF.tif",
			ffdir, strlen( name ) - 4, name );
		Raster16ToTif16( path, ras, gW, gH, flog );

		RasterFree( ras );
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	gArgs.SetCmdLine( argc, argv );

	ReadParams();

	ProcessRick();

	if( frick )
		fclose( frick );

	for( int i = 0; i < 4; ++i )
		RasterFree( ffras[i] );

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


