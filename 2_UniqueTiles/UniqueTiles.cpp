//
// List any tiles that occur twice in a layer.
//
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<set>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_heq ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_heq {
public:
    char	*infile;
public:
    CArgs_heq() : infile(NULL) {};

    void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_heq	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_heq::SetCmdLine( int argc, char* argv[] )
{
// start log

    flog = FileOpenOrDie( "UniqueTiles.log", "w" );

// parse command line args

    if( argc < 2 ) {
        printf( "Usage: UniqueTiles <xml-file>.\n" );
        exit( 42 );
    }

    for( int i = 1; i < argc; ++i ) {

        // echo to log
        fprintf( flog, "%s ", argv[i] );

        if( argv[i][0] != '-' ) {

            if( !infile )
                infile = argv[i];
        }
        else {
            printf( "Did not understand option '%s'.\n", argv[i] );
            exit( 42 );
        }
    }

// header

    fprintf( flog, "\n\nDuplicate Tiles---\nZ\tPath\n" );
    fflush( flog );
}

/* --------------------------------------------------------------- */
/* PrintDupTiles ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintDupTiles( TiXmlElement* layer, int z )
{
    set<string>		S;
    TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

    for( ; p; p = p->NextSiblingElement() ) {

        string	s = p->Attribute( "file_path" );

        if( S.find( s ) != S.end() )
            fprintf( flog, "%d\t%s\n", z, s.c_str() );
        else
            S.insert( s );
    }
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2()
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

        int	z = atoi( layer->Attribute( "z" ) );

        PrintDupTiles( layer, z );
    }
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

/* ---------------- */
/* Read source file */
/* ---------------- */

    ParseTrakEM2();

/* ---- */
/* Done */
/* ---- */

    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


