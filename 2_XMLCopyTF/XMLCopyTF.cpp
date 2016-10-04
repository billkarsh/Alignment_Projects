//
// Copy tforms from file2=TAffineTable.txt into file1=TrakEM2.xml.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"PipeFiles.h"


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
    char	*xmlfile,
            *tblfile;
public:
    CArgs_xml() : xmlfile(NULL), tblfile(NULL) {};

    void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml			gArgs;
static FILE*				flog = NULL;
static map<MZIDR,TAffine>	M;
static set<int>				Z;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

    flog = FileOpenOrDie( "XMLCopyTF.log", "w" );

// log start time

    time_t	t0 = time( NULL );
    char	atime[32];

    strcpy( atime, ctime( &t0 ) );
    atime[24] = '\0';	// remove the newline

    fprintf( flog, "Start: %s ", atime );

// parse command line args

    if( argc < 3 ) {
        printf(
        "Usage: XMLCopyTF <xml-file> <tbl-file>.\n" );
        exit( 42 );
    }

    for( int i = 1; i < argc; ++i ) {

        // echo to log
        fprintf( flog, "%s ", argv[i] );

        if( argv[i][0] != '-' ) {

            if( !xmlfile )
                xmlfile = argv[i];
            else
                tblfile = argv[i];
        }
        else {
            printf( "Did not understand option [%s].\n", argv[i] );
            exit( 42 );
        }
    }

    fprintf( flog, "\n\n" );
    fflush( flog );
}

/* --------------------------------------------------------------- */
/* CopyMatchingTF ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void CopyMatchingTF( TiXmlElement* layer, int z )
{
    MZIDR			key;
    TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );
    TiXmlElement*	next;

    key.z	= z;
    key.rgn	= 1;

    for( ; p; p = next ) {

        next = p->NextSiblingElement();

        key.id = IDFromPatch( p );

        map<MZIDR,TAffine>::iterator	it = M.find( key );

        if( it != M.end() )
            XMLSetTFVals( p, it->second.t );
        else
            layer->RemoveChild( p );
    }
}

/* --------------------------------------------------------------- */
/* Update -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update()
{
/* ---- */
/* Open */
/* ---- */

    XML_TKEM		xml( gArgs.xmlfile, flog );
    TiXmlElement*	layer	= xml.GetFirstLayer();

/* ------------------------ */
/* Copy matching transforms */
/* ------------------------ */

    for( ; layer; layer = layer->NextSiblingElement() ) {

        int	z = atoi( layer->Attribute( "z" ) );

        if( Z.find( z ) == Z.end() )
            continue;

        CopyMatchingTF( layer, z );
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

/* ------------------------------------- */
/* Load lists of TForms and affected Z's */
/* ------------------------------------- */

    LoadTAffineTbl_AllZ( M, Z, gArgs.tblfile, flog );

/* ------------- */
/* Write new xml */
/* ------------- */

    Update();

/* ---- */
/* Done */
/* ---- */

    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


