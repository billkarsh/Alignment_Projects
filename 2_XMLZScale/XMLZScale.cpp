//
// Multiply each layer's z tag by given -mul= factor.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_xml ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xml {

public:
    double	mul, zmin, zmax;
    char	*infile;

public:
    CArgs_xml()
    {
        mul		= 1.0;
        zmin	= 0.0;
        zmax	= 32768.0;
        infile	= NULL;
    };

    void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

    flog = FileOpenOrDie( "XMLZScale.log", "w" );

// log start time

    time_t	t0 = time( NULL );
    char	atime[32];

    strcpy( atime, ctime( &t0 ) );
    atime[24] = '\0';	// remove the newline

    fprintf( flog, "Start: %s ", atime );

// parse command line args

    if( argc < 5 ) {
        printf( "Usage: XMLZScale <xml-file1> -mul=f -zmin=i -zmax=j.\n" );
        exit( 42 );
    }

    for( int i = 1; i < argc; ++i ) {

        // echo to log
        fprintf( flog, "%s ", argv[i] );

        if( argv[i][0] != '-' )
            infile = argv[i];
        else if( GetArg( &mul, "-mul=%lf", argv[i] ) )
            ;
        else if( GetArg( &zmin, "-zmin=%lf", argv[i] ) )
            ;
        else if( GetArg( &zmax, "-zmax=%lf", argv[i] ) )
            ;
        else {
            printf( "Did not understand option [%s].\n", argv[i] );
            exit( 42 );
        }
    }

    fprintf( flog, "\n\n" );
    fflush( flog );
}

/* --------------------------------------------------------------- */
/* Process ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Process()
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

        double	z = atof( layer->Attribute( "z" ) );

        if( z > gArgs.zmax )
            break;

        if( z < gArgs.zmin )
            continue;

        z *= gArgs.mul;

        layer->SetDoubleAttribute( "z", z );
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
/* ------------------ */
/* Parse command line */
/* ------------------ */

    gArgs.SetCmdLine( argc, argv );

/* ------------- */
/* Write new xml */
/* ------------- */

    Process();

/* ---- */
/* Done */
/* ---- */

    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


