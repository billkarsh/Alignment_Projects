//
// Scenario:
// Alignment has been run on input xml file 'src.xml' to make
// output file 'dst.xml' (usually the xml file created by lsq).
// We want the oids of the new file to be the same as they were
// on the same objects, but some tiles may have been kicked out,
// and lsq just assigns oids consecutively. Therefore...
//
// Run XMLCopyOids 'dst.xml' 'src.xml' to force like objects
// to get the original oids from 'src.xml'.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<map>
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
	char	*dst,
			*src;

public:
	CArgs_xml()
	{
		dst	= NULL;
		src	= NULL;
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

	flog = FileOpenOrDie( "XMLCopyOids.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: XMLCopyOids <dst.xml> <src.xml> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );
	}

	dst = argv[1];
	src = argv[2];

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* GetTileOids --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetTileOids(
	map<string,int>	&M,
	TiXmlElement*	layer )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		string	s = p->Attribute( "file_path" );

		M[s] = atoi( p->Attribute( "oid" ) );
	}
}

/* --------------------------------------------------------------- */
/* UpdateTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTiles(
	TiXmlElement*	layer,
	map<string,int>	&M )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		string	s = p->Attribute( "file_path" );

		map<string,int>::iterator	mi = M.find( s );

		if( mi == M.end() ) {
			fprintf( flog,
			"Dest tile not found in src '%s'.\n", s.c_str() );
			exit( 42 );
		}

		p->SetAttribute( "oid", mi->second );
	}
}

/* --------------------------------------------------------------- */
/* CopyOids ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void CopyOids()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xmlD( gArgs.dst, flog );
	TiXmlElement*	layerD	= xmlD.GetFirstLayer();

	XML_TKEM		xmlS( gArgs.src, flog );
	TiXmlElement*	layerS	= xmlS.GetFirstLayer();

/* ------- */
/* Process */
/* ------- */

	for( ;
		layerD && layerS;
		layerD = layerD->NextSiblingElement(),
		layerS = layerS->NextSiblingElement() ) {

		int	zD = atoi( layerD->Attribute( "z" ) ),
			zS = atoi( layerS->Attribute( "z" ) );

		if( zD != zS ) {
			fprintf( flog,
			"Unexpected layer mismatch [zD zS] = [%d %d].\n",
			zD, zS );
			exit( 42 );
		}

		/* -------------- */
		/* Copy layer oid */
		/* -------------- */

		int	oidS = atoi( layerS->Attribute( "oid" ) );

		layerD->SetAttribute( "oid", oidS );

		/* ------------------ */
		/* Transfer tile oids */
		/* -------------------*/

		map<string,int>	M;

		GetTileOids( M, layerS );
		UpdateTiles( layerD, M );
	}

/* ---- */
/* Save */
/* ---- */

	xmlD.Save( "xmltmp.txt", true );
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

	CopyOids();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


