//
// "FF and scale tomography images"
//
// Apply flatfield to all images in a rick file.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"CTForm.h"


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
static TForm	gT[4];
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
// chan=0,scale=1,dx=0,dy=0,fullpath
// ... as many as used chans, up to 4
//
static void ReadParams()
{
	FILE*		f = FileOpenOrDie( gArgs.prmfile, "r" );
	CLineScan	LS;
	char		buf[2048];

// scann the TIF folder path

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "TIFpath=%s", tifdir ) )
		exit( 42 );
	fprintf( flog, "TIFpath=%s\n", tifdir );

// create the FF folder next to it

	strcpy( ffdir, tifdir );
	char	*s = strrchr( ffdir, '/' );
	strcpy( s + 1, "FF" );
	DskCreateDir( ffdir, flog );

// scan rick file name and open it

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "rick=%s", buf ) )
		exit( 42 );
	fprintf( flog, "rick=%s\n", buf );
	frick = FileOpenOrDie( buf, "r" );

// now for each channel directive

	while( LS.Get( f ) > 0 ) {

		double	scl, x, y;
		int		chan, np;

		// scan parameters

		if( 5 != sscanf( LS.line,
			"chan=%d,scale=%lf,dx=%lf,dy=%lf,%s",
			&chan, &scl, &x, &y, buf ) ) {

			break;
		}

		fprintf( flog, "chan=%d,scale=%g,dx=%g,dy=%g,%s\n",
			chan, scl, x, y, buf );

		// compute FF average

		ffras[chan] = Raster16FromTif16( buf, gW, gH, flog );

		np = gW * gH;

		ffave[chan] = 0;
		for( int i = 0; i < np; ++i ) {
			if( ffras[chan][i] <= 0 )
				ffras[chan][i] = 1;
			ffave[chan] += ffras[chan][i];
		}
		ffave[chan] /= np;

		// compute gT

		gT[chan].NUSetScl( scl );
		gT[chan].SetXY( x, y );
	}

	fprintf( flog, "\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* MagChannel ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MagChannel( uint16* ras, const TForm &T )
{
	int				np = gW * gH;
	vector<double>	v( np );
	vector<double>	I( np, 0.0 );

	for( int i = 0; i < np; ++i )
		v[i] = ras[i];

	for( int i = 0; i < np; ++i ) {

		int		y = i / gW,
				x = i - gW * y;
		Point	p( x, y );

		T.Transform( p );
		DistributePixel( p.x, p.y, v[i], I, gW, gH );
	}

	for( int i = 0; i < np; ++i )
		ras[i] = (uint16)I[i];
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

		// read an image name

		sscanf( LS.line, "%s", name );

		// get its channel

		char	*s = strrchr( name, '_' );
		chan = atoi( s + 1 );

		// only interested in given channels

		if( chan < 0 || chan > 3 || !ffras[chan] )
			continue;

		// now read the image

		sprintf( path, "%s/%s", tifdir, name );
		uint16*	ras = Raster16FromTif16( path, gW, gH, flog );

		if( !ras )
			break;

		// flat-field the image

		for( int i = 0; i < np; ++i )
			ras[i] = uint16(ras[i]*ffave[chan]/ffras[chan][i]);

		// apply TForm

		if( gT[chan].t[0] != 1.0 || gT[chan].t[2] || gT[chan].t[5] )
			MagChannel( ras, gT[chan] );

		// write it out

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


