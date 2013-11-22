//
// Remove from xml file all layers listed in file2.
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
/* CArgs_xml ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xml {
public:
	char	*xmlfile,
			*listfile;
public:
	CArgs_xml() : xmlfile(NULL), listfile(NULL) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml	gArgs;
static FILE*		flog = NULL;
static set<int>		Z;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "XMLKillZ.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: XMLKillZ <xml-file> <list-file>.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !xmlfile )
				xmlfile = argv[i];
			else
				listfile = argv[i];
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
/* GetZList ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetZList()
{
	FILE		*f = FileOpenOrDie( gArgs.listfile, "r" );
	CLineScan	LS;

	while( LS.Get( f ) > 0 ) {

		int	z;

		if( 1 == sscanf( LS.line, "%d", &z ) )
			Z.insert( z );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Edit ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Edit()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.xmlfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();
	TiXmlNode*		lyrset	= layer->Parent();

/* ---------------------- */
/* Remove matching layers */
/* ---------------------- */

	TiXmlElement*	next;

	for( ; layer; layer = next ) {

		next = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( Z.find( z ) == Z.end() )
			continue;

		lyrset->RemoveChild( layer );
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

/* ---------------- */
/* Load list of Z's */
/* ---------------- */

	GetZList();

/* ------------- */
/* Write new xml */
/* ------------- */

	Edit();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


