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
static FILE*	flog = NULL;
static char		tifdir[2048],
                ffdir[2048],
                rick[2048];
static uint16*	ffras[4]	= {NULL,NULL,NULL,NULL};
static double	ffave[4]	= {0,0,0,0};
static int		ischn[4]	= {0,0,0,0};
static uint32	gW			= 0,
                gH			= 0;	// universal pic dims






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

// name the FF folder

    sprintf( ffdir, "%s/FF_c2/", tifdir );
    fprintf( flog, "FFpath=%s\n", ffdir );

// scan rick file name

    if( LS.Get( f ) <= 0 || 1 != sscanf( LS.line, "rick=%s", rick ) )
        exit( 42 );
    fprintf( flog, "rick=%s\n", rick );

// now for each channel directive

    while( LS.Get( f ) > 0 ) {

        int		chan, np;

        // scan parameters

        if( 2 != sscanf( LS.line,
            "chan=%d,%s",
            &chan, buf ) ) {

            break;
        }

        fprintf( flog,
        "chan=%d,%s\n",
        chan, buf );

        ischn[chan] = true;

        // external file: compute FF average

        ffras[chan] = Raster16FromTif16( buf, gW, gH, flog );

        if( !ffras[chan] ) {
            fprintf( flog, "Bad ff file, chan %d.\n", chan );
            exit(0);
        }

        np = gW * gH;

        for( int i = 0; i < np; ++i ) {

            if( ffras[chan][i] == 0 )
                ffras[chan][i] = 1;

            ffave[chan] += ffras[chan][i];
        }

        ffave[chan] /= np;
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

            ras[i] = uint16(ras[i]*ffave[chan]/ffras[chan][i]);
        }
    }
}

/* --------------------------------------------------------------- */
/* DoChannel ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DoChannel( int chan )
{
// Reopen rick file for each channel

    FILE		*frick = FileOpenOrDie( rick, "r" );
    CLineScan	LS;

// For each line...
// Get its image-name, x, y
// Do pixel ops on that image and write it to FF folder

    while( LS.Get( frick ) > 0 ) {

        char	path[2048], name[2048];

        // Get native line data
        sscanf( LS.line, "%s", name );

        // This channel?
        char	*c = strrchr( name, 'c' );
        if( atoi( c + 1 ) != chan )
            continue;

        // Read the image
        uint16*	ras = Raster16FromTif16( name, gW, gH, flog );

        if( !ras ) {
            fprintf( flog, "Missing image=[%s]\n", name );
            continue;
        }

        // Flat-field the image
        FFChannel( ras, chan );

        // Write image file
        sprintf( path, "%s%s", ffdir, FileNamePtr( name ) );
        Raster16ToTif16( path, ras, gW, gH, flog );
        RasterFree( ras );
    }

    fclose( frick );
}

/* --------------------------------------------------------------- */
/* ChannelLoop --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ChannelLoop()
{
    for( int chan = 0; chan < 4; ++chan ) {

        if( ischn[chan] )
            DoChannel( chan );
    }
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


