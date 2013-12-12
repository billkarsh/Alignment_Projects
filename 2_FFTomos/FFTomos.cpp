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
#include	"TAffine.h"
#include	"CMask.h"

#include	<string.h>


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
	CArgs() : prmfile(NULL) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static CMask	Mask;
static FILE*	flog = NULL;
static char		tifdir[2048],
				ffdir[2048],
				mskdir[2048],
				rick[2048];
static double	gscale		= 1.0;
static uint16*	ffras[4]	= {NULL,NULL,NULL,NULL};
static double	ffave[4]	= {0,0,0,0};
static double	ffstd[4]	= {0,0,0,0};
static int		fford[4]	= {0,0,0,0};	// leg poly order
static int		ffoff[4]	= {0,0,0,0};	// use vals <= mode+offset
static int		ischn[4]	= {0,0,0,0};
static int		useT[4]		= {0,0,0,0};
static TAffine	gT[4];
static uint32	gW			= 0,
				gH			= 0;	// universal pic dims
static int		gped		= 0,
				gmchn		= -1;






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
// scale=1
// ped=0
// mask=T,0,path/myparams.xml
// chan=0,useAff=T,Aff=[1 0 0 0 1 0],fullpath
// ... as many as used chans, up to 4
//
// If fullpath for flat-field has form LEG:order:offset:mean:std
// then Legendre polys are used instead of external file.
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

// scan rick file name

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "rick=%s", rick ) )
		exit( 42 );
	fprintf( flog, "rick=%s\n", rick );

// scale

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "scale=%lf", &gscale ) )
		exit( 42 );
	fprintf( flog, "scale=%g\n", gscale );

// pedestal

	if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "ped=%d", &gped ) )
		exit( 42 );
	fprintf( flog, "ped=%d\n", gped );

// mask

	char	cUse;

	if( LS.Get( f ) <= 0 ||
		3 != sscanf( LS.line, "mask=%c,%d,%s", &cUse, &gmchn, buf ) ) {

		exit( 42 );
	}
	fprintf( flog, "mask=%c,%d,%s\n", cUse, gmchn, buf );

	if( toupper( cUse ) == 'T' ) {

		Mask.ReadParamFile( flog, buf );

		strcpy( mskdir, tifdir );
		char	*s = strrchr( mskdir, '/' );
		strcpy( s + 1, "Mask" );
		DskCreateDir( mskdir, flog );
	}
	else
		gmchn = -1;

// now for each channel directive

	while( LS.Get( f ) > 0 ) {

		double	A[6];
		int		chan, np;
		char	cUse;

		// scan parameters

		if( 9 != sscanf( LS.line,
			"chan=%d,useAff=%c,Aff=[%lf%lf%lf%lf%lf%lf],%s",
			&chan, &cUse,
			&A[0], &A[1], &A[2], &A[3], &A[4], &A[5],
			buf ) ) {

			break;
		}

		fprintf( flog,
		"chan=%d,useAff=%c,Aff=[%f %f %f %f %f %f],%s\n",
		chan, cUse,
		A[0], A[1], A[2], A[3], A[4], A[5],
		buf );

		ischn[chan] = true;

		// determine FF method

		if( 4 == sscanf( buf, "LEG:%d:%d:%lf:%lf",
			&fford[chan], &ffoff[chan],
			&ffave[chan], &ffstd[chan] ) ) {

			fprintf( flog, "ff chan %d using LEG [%d,%d,%f,%f].\n",
			chan,
			fford[chan], ffoff[chan],
			ffave[chan], ffstd[chan] );

		}
		else {

			// external file: compute FF average

			ffras[chan] = Raster16FromTif16( buf, gW, gH, flog );

			np = gW * gH;

			for( int i = 0; i < np; ++i ) {

				if( ffras[chan][i] >= gped )
					ffras[chan][i] -= gped;

				if( ffras[chan][i] == 0 )
					ffras[chan][i] = 1;

				ffave[chan] += ffras[chan][i];
			}

			ffave[chan] /= np;
		}

		// compute gT

		if( useT[chan] = (toupper( cUse ) == 'T') )
			gT[chan] = TAffine( A );
	}

	fprintf( flog, "\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* FFChannel ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void FFChannel( uint16* ras, int chan )
{
	int	np = gW * gH;

	if( ffras[chan] ) {

		// external file

		for( int i = 0; i < np; ++i ) {

			if( ras[i] >= gped )
				ras[i] -= gped;

			ras[i] = uint16(ras[i]*ffave[chan]/ffras[chan][i]);
		}
	}
	else {

		// ped subtract

		for( int i = 0; i < np; ++i ) {

			if( ras[i] >= gped )
				ras[i] -= gped;
		}

		// Legendre polys

		vector<double>	V;

		LegPolyFlatten( V, ras, gW, gH, fford[chan], ffoff[chan] );

		// rescale to given mean, stddev

		for( int i = 0; i < np; ++i ) {

			int	pix = int(ffave[chan] + V[i] * ffstd[chan]);

			if( pix < 0 )
				pix = 0;
			else if( pix > 65535 )
				pix = 65535;

			ras[i] = pix;
		}
	}
}

/* --------------------------------------------------------------- */
/* MagChannel ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MagChannel( uint16* ras, const TAffine &T )
{
	TAffine			I;
	int				np = gW * gH;
	vector<uint16>	src( np );

	I.InverseOf( T );
	memcpy( &src[0], ras, np * sizeof(uint16) );

	for( int i = 0; i < np; ++i ) {

		int		y = i / gW,
				x = i - gW * y;
		Point	p( x, y );

		I.Transform( p );

		if( p.x >= 0 && p.x < gW - 1 &&
			p.y >= 0 && p.y < gH - 1 ) {

			ras[x+gW*y] =
			(uint16)SafeInterp( p.x, p.y, &src[0], gW, gH );
		}
	}
}

/* --------------------------------------------------------------- */
/* DoChannel ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DoChannel( FILE* fout, int &z, int chan )
{
	char	curwel[8] = {0};

// Reopen rick file for each channel

	FILE		*frick = FileOpenOrDie( rick, "r" );
	CLineScan	LS;

// For each line...
// Get its image-name, x, y
// Do pixel ops on that image and write it to FF folder
// Output new full path to modified image
// Output rescaled x, y and updated z

	while( LS.Get( frick ) > 0 ) {

		char	path[2048], name[64];
		double	x, y;
		int		lname;

		// Get native line data
		sscanf( LS.line, "%s%lf%lf", name, &x, &y );

		lname = strlen( name );

		// Force channel name
		name[lname - 5] = '0' + chan;

		// Get well tag length and test change
		int	len = strchr( name, '_' ) - name;

		if( strncmp( curwel, name, len ) ) {
			// changed
			sprintf( curwel, "%.*s", len, name );
			++z;
		}

		// Read the image
		sprintf( path, "%s/%s", tifdir, name );
		uint16*	ras = Raster16FromTif16( path, gW, gH, flog );

		if( !ras ) {
			fprintf( flog, "Missing image=[%s]\n", path );
			continue;
		}

		// Flat-field the image
		FFChannel( ras, chan );

		// Apply TForm
		if( useT[chan] )
			MagChannel( ras, gT[chan] );

		// Write image file
		sprintf( path, "%s/%.*s.FF.tif", ffdir, lname - 4, name );
		Raster16ToTif16( path, ras, gW, gH, flog );

		// Write output line
		fprintf( fout, "%s\t%.2f\t%.2f\t%d\n",
			path, x * gscale, y * gscale, z );

		// Make mask file
		if( chan == gmchn ) {

			sprintf( path, "%s/%.*s.Mask.tif",
				mskdir, lname - 4, name );

			Mask.MaskFromImage( path, ras, gW, gH );
		}

		RasterFree( ras );
	}

	fclose( frick );
}

/* --------------------------------------------------------------- */
/* WriteMaskLines ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteMaskLines( FILE* fout, int &z, int chan )
{
	char	curwel[8] = {0};

// Reopen rick file for mask channel

	FILE		*frick = FileOpenOrDie( rick, "r" );
	CLineScan	LS;

// For each line...
// Get its image-name, x, y
// Output new full path to mask image
// Output rescaled x, y and updated z

	while( LS.Get( frick ) > 0 ) {

		char	path[2048], name[64];
		double	x, y;
		int		lname;

		// Get native line data
		sscanf( LS.line, "%s%lf%lf", name, &x, &y );

		lname = strlen( name );

		// Force channel name
		name[lname - 5] = '0' + chan;

		// Get well tag length and test change
		int	len = strchr( name, '_' ) - name;

		if( strncmp( curwel, name, len ) ) {
			// changed
			sprintf( curwel, "%.*s", len, name );
			++z;
		}

		// Write output line
		sprintf( path, "%s/%.*s.Mask.tif", 	mskdir, lname - 4, name );

		fprintf( fout, "%s\t%.2f\t%.2f\t%d\n",
			path, x * gscale, y * gscale, z );
	}

	fclose( frick );
}

/* --------------------------------------------------------------- */
/* ChannelLoop --------------------------------------------------- */
/* --------------------------------------------------------------- */

// For ea channel specified in input file:
// Write block of wells/tiles for that channel.
//
// Note that z values advance with changes in well name...
//	and changes in channel.
//
static void ChannelLoop()
{
// Name and open the one output file. The original name
// looks like "...TrackEM2_Chni.txt" and we will change
// to "...TrackEM2_FFall.txt"

	char	buf[2048];
	sprintf( buf, "%.*sFFall.txt", strlen( rick ) - 8, rick );

	FILE	*fout = FileOpenOrDie( buf, "w" );

	int	z = -1;	// changes with chan/well

// Write lines for the FF images

	for( int chan = 0; chan < 4; ++chan ) {

		if( ischn[chan] )
			DoChannel( fout, z, chan );
	}

// Write lines at end for masks

	if( gmchn >= 0 )
		WriteMaskLines( fout, z, gmchn );

	fclose( fout );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	gArgs.SetCmdLine( argc, argv );

	ReadParams();

	ChannelLoop();

	for( int i = 0; i < 4; ++i )
		RasterFree( ffras[i] );

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


