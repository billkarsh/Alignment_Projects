//
// Check how images are distributed among subfolders of
// Plate1_0. That is, usually the unique folder is 'TIF'
// but here we check in each layer if the images are in
// a folder of a different name, or in multiple folders.
//
// Example command:
// >XMLExists xxx.xml -zmin=0 -zmax=10
//


#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<set>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_xex ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xex {

public:
    char	*infile;
    int		zmin, zmax;

public:
    CArgs_xex()
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

static CArgs_xex		gArgs;
static FILE*			flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xex::SetCmdLine( int argc, char* argv[] )
{
// start log

    flog = FileOpenOrDie( "XMLFolders.log", "w" );

// parse command line args

    if( argc < 2 ) {
        printf( "Usage: XMLExists <xml-file> [options].\n" );
        exit( 42 );
    }

    for( int i = 1; i < argc; ++i ) {

        // echo to log
        fprintf( flog, "%s ", argv[i] );

        if( argv[i][0] != '-' ) {

            if( !infile )
                infile = argv[i];
        }
        else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
            ;
        else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
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
/* CollectTileDirs ----------------------------------------------- */
/* --------------------------------------------------------------- */

static int CollectTileDirs( set<string> &dirs, TiXmlElement* layer )
{
    TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

    for( ; ptch; ptch = ptch->NextSiblingElement() ) {

        char	buf[4096], sub[256];
        char	*start, *end;

        sprintf( buf, "%s", ptch->Attribute( "file_path" ) );

        start = strstr( buf, "Plate1_0" ) + 9;
        end = strrchr( buf, '/' );

        sprintf( sub, "%.*s", end - start, start );
        dirs.insert( string( sub ) );
    }

    return dirs.size();
}

/* --------------------------------------------------------------- */
/* ScanXML ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScanXML()
{
/* ---- */
/* Open */
/* ---- */

    XML_TKEM		xml( gArgs.infile, flog );
    TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---- */
/* Scan */
/* ---- */

    for( ; layer; layer = layer->NextSiblingElement() ) {

        int	z = atoi( layer->Attribute( "z" ) );

        if( z > gArgs.zmax )
            break;

        if( z < gArgs.zmin )
            continue;

        if( !(z % 100) )
            printf( "z=%6d\n", z );

        set<string>	dirs;
        int			nd = CollectTileDirs( dirs, layer );

        if( !nd )
            fprintf( flog, "z=%d\t*** No Folder\n", z );
        else if( nd == 1 ) {
            set<string>::iterator	it = dirs.begin();
            if( strcmp( (*it).c_str(), "TIF_HEQ" ) )
                fprintf( flog, "z=%d\t%s\n", z, (*it).c_str() );
        }
        else {
            set<string>::iterator	it = dirs.begin();
            for( int i = 0; i < nd; ++i, ++it ) {
                fprintf( flog, "z=%d\t%s%s",
                z, (!i ? "" : "; "), (*it).c_str() );
            }
            fprintf( flog, "\n" );
        }
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

/*----- */
/* Scan */
/*----- */

    ScanXML();

/* ---- */
/* Done */
/* ---- */

exit:
    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


