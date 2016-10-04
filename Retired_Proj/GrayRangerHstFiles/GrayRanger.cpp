

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"Maths.h"

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {
public:
    int		z, id;
public:
    Picture( const char* name, int _z );

    bool operator < ( const Picture &rhs ) const
        {return z < rhs.z || (z == rhs.z && id < rhs.id);};
};


class Fix {
public:
    int	z, min, max;
public:
    Fix( int z, int min, int max )
    : z(z), min(min), max(max) {};
};

/* --------------------------------------------------------------- */
/* CArgs_gray ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_gray {

public:
    char	*infile;
    int		zmin,
            zmax,
            chn;

public:
    CArgs_gray()
    {
        infile	= NULL;
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
static uint32		gW = 0,	gH = 0;		// universal pic dims






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

    if( argc < 3 ) {
        printf( "Usage: GrayRanger <xml-file> -chan=c [options].\n" );
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
        else if( GetArg( &chn, "-chan=%d", argv[i] ) )
            ;
        else {
            printf( "Did not understand option '%s'.\n", argv[i] );
            exit( 42 );
        }
    }

    fprintf( flog, "\n\n" );
    fflush( flog );
}

/* --------------------------------------------------------------- */
/* Picture ------------------------------------------------------- */
/* --------------------------------------------------------------- */

Picture::Picture( const char* name, int _z )
{
    const char	*s;

    s = FileNamePtr( name );
    s = strchr( s, '_' ) + 1;

    z		= _z;
    id		= atoi( s );
}

/* --------------------------------------------------------------- */
/* GetTiles ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetTiles(
    vector<Picture>	&vp,
    TiXmlElement*	layer,
    int				z )
{
    TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

    for( ; ptch; ptch = ptch->NextSiblingElement() ) {

        vp.push_back( Picture( ptch->Attribute( "file_path" ), z ) );
    }
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<Picture> &vp )
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

        GetTiles( vp, layer, z );
    }
}

/* --------------------------------------------------------------- */
/* GetLayerLimits ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Given starting index i0, which selects a layer z, set iN to be
// one beyond the highest index of a picture having same z. This
// makes loop limits [i0,iN) exclusive.
//
// If i0 or iN are out of bounds, both are set to -1.
//
static void GetLayerLimits(
    const vector<Picture>	&vp,
    int						&i0,
    int						&iN )
{
    int	np = vp.size();

    if( i0 < 0 || i0 >= np ) {

        i0 = -1;
        iN = -1;
    }
    else {

        int	Z = vp[i0].z;

        for( iN = i0 + 1; iN < np && vp[iN].z == Z; ++iN )
            ;
    }
}

/* --------------------------------------------------------------- */
/* ScaleThisLayer ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ScaleThisLayer(
    vector<Fix>				&fix,
    const vector<Picture>	&vp,
    int						i0,
    int						iN )
{
    const int		nbins = 65536;	// also max val
    vector<double>	bins( nbins, 0.0 );
    double			uflo = 0.0,
                    oflo = 0.0;
    int				T,
                    imin,
                    smin,
                    smax;

// sum layer's histograms

    for( int i = i0; i < iN; ++i ) {

        char	name[2048];
        FILE	*f;

        sprintf( name, "HST/HST_%d_%d_%d.bin",
            vp[i].z, vp[i].id, gArgs.chn );

        if( f = fopen( name, "rb" ) ) {

            vector<double>	binsi( nbins );
            double			ufloi,
                            ofloi;

            fread( &ufloi, sizeof(double), 1, f );
            fread( &ofloi, sizeof(double), 1, f );
            fread( &binsi[0], sizeof(double), nbins, f );

            uflo += ufloi;
            oflo += ofloi;

            for( int j = 0; j < nbins; ++j )
                bins[j] += binsi[j];

            fclose( f );
        }
    }

// smin is between lowest val and 2 sdev below mode
// Omit highest bins to avoid detector saturation

    imin = IndexOfMaxVal( &bins[0], nbins - 100 );
    T	 = int(2.0 * sqrt( imin ));
    smax = imin - T;
    smin = FirstNonzero( &bins[0], nbins );

    if( smax < 0 )
        smax = 0;

    if( smin < 0 )
        smin = 0;

    smin = (smin + smax) / 2;

// Get threshold for bkg-frg segmentation using Otsu method,
// but exclude values lower than imin + 2 sdev.

    imin += T;
    T = int(OtsuThresh( &bins[imin], nbins - imin, imin, nbins ));

// Now, we aren't going to find objects with the Otsu threshold,
// however, we are going to set the scale maximum as 99% of the
// foreground counts.

    smax = PercentileBin( &bins[T], nbins - T, 0.99 );

    fix.push_back( Fix( vp[i0].z, smin, smax ) );
}

/* --------------------------------------------------------------- */
/* ScaleAllLayers ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Loop over layers, process images.
//
static void ScaleAllLayers(
    vector<Fix>				&fix,
    const vector<Picture>	&vp )
{
    int	i0, iN;

    GetLayerLimits( vp, i0 = 0, iN );

    while( iN != -1 ) {

        fprintf( flog, "\n---- Layer %d ----\n", vp[i0].z );

        ScaleThisLayer( fix, vp, i0, iN );

        printf( "Done layer %d\n", vp[i0].z );

        GetLayerLimits( vp, i0 = iN, iN );
    }
}

/* --------------------------------------------------------------- */
/* FixTiles ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void FixTiles( TiXmlElement* layer, const Fix &fix )
{
    TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

    for( ; ptch; ptch = ptch->NextSiblingElement() ) {

        ptch->SetAttribute( "min", fix.min );
        ptch->SetAttribute( "max", fix.max );
    }
}

/* --------------------------------------------------------------- */
/* WriteXML ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteXML( const vector<Fix> &fix )
{
/* ---- */
/* Open */
/* ---- */

    XML_TKEM		xml( gArgs.infile, flog );
    TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---------- */
/* Fix layers */
/* ---------- */

    int	nz = fix.size();

    // for each layer we wish to fix...
    for( int iz = 0; iz < nz; ++iz ) {

        // advance xml to that layer
        while( atoi( layer->Attribute( "z" ) ) < fix[iz].z )
            layer = layer->NextSiblingElement();

        if( !layer )
            break;

        FixTiles( layer, fix[iz] );
    }

/* ---- */
/* Save */
/* ---- */

    xml.Save( "xmltmp.txt", true );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
    vector<Picture>	vp;
    vector<Fix>		fix;

/* ------------------ */
/* Parse command line */
/* ------------------ */

    gArgs.SetCmdLine( argc, argv );

/* ---------------- */
/* Read source file */
/* ---------------- */

    ParseTrakEM2( vp );

    fprintf( flog, "Got %d images.\n", vp.size() );

    if( !vp.size() )
        goto exit;

/* ------------------ */
/* Sort by Z and tile */
/* ------------------ */

    sort( vp.begin(), vp.end() );

/* ----------------- */
/* Calculate scaling */
/* ----------------- */

    ScaleAllLayers( fix, vp );

/* ------------- */
/* Write new xml */
/* ------------- */

    WriteXML( fix );

/* ---- */
/* Done */
/* ---- */

exit:
    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


