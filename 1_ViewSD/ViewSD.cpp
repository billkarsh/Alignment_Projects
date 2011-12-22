//
// > ViewSD file.xml sdall.txt -sdmin=0 -sdmax=800 -zmin=i -zmax=j
//
// Make edited copy of file.xml having only given layer range and
// only given sd range (sd's from sdall.txt).
//

#include	"Cmdline.h"
#include	"CRegexID.h"
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

// ----------------------------------------

// Map {z,id} -> sd

class MZID {

public:
	int	z, id;

public:
	MZID() {};

	bool operator < (const MZID &rhs) const
		{
			if( z < rhs.z )
				return true;
			if( z > rhs.z )
				return false;

			return id < rhs.id;
		};

	bool operator == (const MZID &rhs) const
		{return z == rhs.z && id == rhs.id;};
};

/* --------------------------------------------------------------- */
/* CArgs_xml ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xml {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID	re_id;

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

	int IDFromPatch( TiXmlElement *p );
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

	char	*pat;

	re_id.Set( "_Nex.mrc" );

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
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &sdmin, "-sdmin=%d", argv[i] ) )
			;
		else if( GetArg( &sdmax, "-sdmax=%d", argv[i] ) )
			;
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromPatch( TiXmlElement *p )
{
	const char	*name = p->Attribute( "title" );
	int			id;

	if( !re_id.Decode( id, name ) ) {
		printf( "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
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
/* Edit ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Edit()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.xmlfile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog,
		"No trakEM2 node [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog,
		"No t2_layer [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

/* --------- */
/* Do layers */
/* --------- */

	TiXmlNode		*lyrset	= layer->Parent();
	TiXmlElement	*nextL	= NULL;
	FILE			*fres	= FileOpenOrDie( "Resin.txt", "w" );

	for( ; layer; layer = nextL ) {

		nextL = layer->NextSiblingElement();

		MZID	key;

		key.z = atoi( layer->Attribute( "z" ) );

		if( key.z > gArgs.zmax ||
			key.z < gArgs.zmin ||
			Z.find( key.z ) == Z.end() ) {

			lyrset->RemoveChild( layer );
			continue;
		}

		TiXmlElement	*nextT = NULL;

		for(
			TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
			p;
			p = nextT ) {

			nextT = p->NextSiblingElement();

			key.id = gArgs.IDFromPatch( p );

			map<MZID,int>::iterator	it = M.find( key );

			if( it == M.end() ||
				it->second < gArgs.sdmin ||
				it->second > gArgs.sdmax ) {

				layer->RemoveChild( p );
			}
			else
				fprintf( fres, "%d\t%d\n", key.z, key.id );
		}

		if( !layer->FirstChildElement( "t2_patch" ) )
			lyrset->RemoveChild( layer );
	}

	fclose( fres );

/* ---- */
/* Save */
/* ---- */

	doc.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.xmlfile, "xmltmp.txt" );
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

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


