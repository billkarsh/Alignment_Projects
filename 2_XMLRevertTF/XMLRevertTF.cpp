//
// Copy tforms from file2/layer=src to file1/layer=dst-list
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"TAffine.h"

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
	char	*dstfile,
			*srcfile;
	int		zsrc;

public:
	CArgs_xml()
	{
		dstfile	= NULL;
		srcfile	= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml		gArgs;
static FILE*			flog = NULL;
static map<int,TAffine>	M;
static set<int>			Z;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "XMLRevertTF.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: XMLRevertTF <dst-file> <-dst=list>"
		" <src-file> <-src=z>.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		vector<int>	vi;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !dstfile )
				dstfile = argv[i];
			else
				srcfile = argv[i];
		}
		else if( GetArgList( vi, "-dst=", argv[i] ) ) {

			for( int i = 0; i < vi.size(); ++i )
				Z.insert( vi[i] );
		}
		else if( GetArg( &zsrc, "-src=%d", argv[i] ) )
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
/* GetTAffines --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetTAffines( TiXmlElement* layer )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		TAffine	T;
		int		id = IDFromPatch( p );

		T.ScanTrackEM2( p->Attribute( "transform" ) );
		M[id] = T;
	}
}

/* --------------------------------------------------------------- */
/* GetSrc -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetSrc()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.srcfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------------- */
/* Move up to src layer */
/* -------------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z == gArgs.zsrc ) {
			GetTAffines( layer );
			return;
		}
	}

	fprintf( flog, "Src layer not found.\n" );
	exit( 42 );
}

/* --------------------------------------------------------------- */
/* UpdateTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTiles( TiXmlElement* layer )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		int	id = IDFromPatch( p );

		map<int,TAffine>::iterator	it = M.find( id );

		if( it != M.end() )
			XMLSetTFVals( p, it->second.t );
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

	XML_TKEM		xml( gArgs.dstfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ------------------------ */
/* Copy matching transforms */
/* ------------------------ */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( Z.find( z ) == Z.end() )
			continue;

		UpdateTiles( layer );
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

/* -------------- */
/* Get src tforms */
/* -------------- */

	GetSrc();

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


