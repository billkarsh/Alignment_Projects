//
// > ViewSD file.xml sdall.txt -sdmin=0 -sdmax=800 -zmin=i -zmax=j
//
// Make edited copy of file.xml having only given layer range and
// only given sd range (sd's from sdall.txt).
//

#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<map>
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
			*sdfile;
	int		zmin,
			zmax,
			sdmin,
			sdmax;

public:
	CArgs_xml()
	{
		xmlfile	= NULL;
		sdfile	= NULL;
		zmin	= 0;
		zmax	= 32768;
		sdmin	= 0;
		sdmax	= 0;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml		gArgs;
static FILE*			flog = NULL;
static map<MZID,int>	M;
static set<int>			Z;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "ViewSD.log", "w" );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: ViewSD <xml-file> <sd-file> -sd= [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( !xmlfile )
				xmlfile = argv[i];
			else
				sdfile = argv[i];
		}
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &sdmin, "-sdmin=%d", argv[i] ) )
			;
		else if( GetArg( &sdmax, "-sdmax=%d", argv[i] ) )
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
/* LoadSD -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LoadSD()
{
	FILE		*f	= FileOpenOrDie( gArgs.sdfile, "r", flog );
	CLineScan	LS;

	for(;;) {

		if( LS.Get( f ) <= 0 )
			break;

		MZID	R;
		int		sd;

		sscanf( LS.line, "%d\t%d\t%d", &R.z, &R.id, &sd );

		Z.insert( R.z );

		M[R] = sd;
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* TrimTiles ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void TrimTiles( FILE* fres, TiXmlElement* layer, int z )
{
	MZID			key;
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );
	TiXmlElement*	nextT;

	key.z = z;

	for( ; p; p = nextT ) {

		nextT = p->NextSiblingElement();

		key.id = IDFromPatch( p );

		map<MZID,int>::iterator	it = M.find( key );

		if( it == M.end() ||
			it->second < gArgs.sdmin ||
			it->second > gArgs.sdmax ) {

			layer->RemoveChild( p );
		}
		else
			fprintf( fres, "%d\t%d\n", key.z, key.id );
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
	TiXmlNode*		lyrset	= xml.GetLayerset();
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* --------- */
/* Do layers */
/* --------- */

	TiXmlElement*	nextL;
	FILE			*fres	= FileOpenOrDie( "Resin.txt", "w" );

	for( ; layer; layer = nextL ) {

		nextL = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax ||
			z < gArgs.zmin ||
			Z.find( z ) == Z.end() ) {

			lyrset->RemoveChild( layer );
			continue;
		}

		TrimTiles( fres, layer, z );

		if( !layer->FirstChildElement( "t2_patch" ) )
			lyrset->RemoveChild( layer );
	}

	fclose( fres );

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

	LoadSD();

/* ---- */
/* Edit */
/* ---- */

	Edit();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


