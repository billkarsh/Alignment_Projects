//
// > DelTiles file.xml del.txt
//
// Make edited copy of file.xml, deleting each layer/tile pair
// listed in del.txt.
//

#include	"GenDefs.h"
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
			*txtfile;

public:
	CArgs_xml()
	{
		xmlfile	= NULL;
		txtfile	= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml	gArgs;
static FILE*		flog = NULL;
static set<MZID>	M;
static set<int>		Z;
static int			zlo = 32768;
static int			zhi = 0;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "DelTiles.log", "w" );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: DelTiles <xml-file> <txt-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( !xmlfile )
				xmlfile = argv[i];
			else
				txtfile = argv[i];
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
/* LoadList ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void LoadList()
{
	FILE		*f	= FileOpenOrDie( gArgs.txtfile, "r", flog );
	CLineScan	LS;

	for(;;) {

		if( LS.Get( f ) <= 0 )
			break;

		MZID	R;

		sscanf( LS.line, "%d\t%d", &R.z, &R.id );

		if( R.z < zlo )
			zlo = R.z;

		if( R.z > zhi )
			zhi = R.z;

		Z.insert( R.z );
		M.insert( R );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* DeleteTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DeleteTiles( TiXmlElement* layer, int z )
{
	MZID			key;
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );
	TiXmlElement*	nextT;

	key.z = z;

	for( ; p; p = nextT ) {

		nextT = p->NextSiblingElement();

		key.id = IDFromPatch( p );

		if( M.find( key ) != M.end() )
			layer->RemoveChild( p );
	}
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

/* --------- */
/* Do layers */
/* --------- */

	TiXmlNode*		lyrset	= layer->Parent();
	TiXmlElement*	nextL;

	for( ; layer; layer = nextL ) {

		nextL = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > zhi )
			break;

		if( z < zlo || Z.find( z ) == Z.end() )
			continue;

		DeleteTiles( layer, z );

		if( !layer->FirstChildElement( "t2_patch" ) )
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
/* Map (z,id) to sd */
/* ---------------- */

	LoadList();

/* ---- */
/* Edit */
/* ---- */

	Edit();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


