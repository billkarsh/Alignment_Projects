//
// Build table of missing items:
// Z	Type	Path
//
// Type is either 'T' for tile or 'D' for directory.
//
// Example command:
// >XMLExists xxx.xml TIF 2 -zmin=0 -zmax=10
//
// Here, 'TIF' is the image folder name following 'Plate1_0'
// and '2' is string for channel 2, but names like 'RGB103'
// could be searched.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"

#include	"tinyxml.h"

#include	<map>
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
	char	*infile,
			*dir,
			*chn;
	int		zmin, zmax;

public:
	CArgs_xex()
	{
		infile	= NULL;
		dir		= NULL;
		chn		= NULL;
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
static map<string,int>	gdirs;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xex::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "XMLExists.log", "w" );

// parse command line args

	if( argc < 4 ) {
		printf( "Usage: XMLExists <xml-file> <dir> <chn> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile )
				infile = argv[i];
			else if( !dir )
				dir = argv[i];
			else
				chn = argv[i];
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

	fprintf( flog, "\n\nZ\tType\tPath\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* IsDir --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool IsDir( string &dir, TiXmlElement* ptch )
{
	char	buf[4096];

	sprintf( buf, "%s", ptch->Attribute( "file_path" ) );
	sprintf( strstr( buf, "Plate1_0" ) + 8, "/%s", gArgs.dir );

	dir = buf;

	map<string,int>::iterator	it = gdirs.find( dir );

	if( it != gdirs.end() )
		return it->second;
	else {

		int	is = DskExists( buf );

		gdirs[dir] = is;
		return is;
	}
}

/* --------------------------------------------------------------- */
/* IsImg --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool IsImg( string &s, TiXmlElement* ptch )
{
	char	til[256], buf[4096];

	sprintf( til, "%s", ptch->Attribute( "title" ) );
	sprintf( strrchr( til, '_' ) + 1, "%s.tif", gArgs.chn );

	sprintf( buf, "%s/%s", s.c_str(), til );
	s = buf;

	return DskBytes( buf ) > 0.0;
}

/* --------------------------------------------------------------- */
/* ScanXML ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScanXML()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	layer = hDoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

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

		// for each tile in this layer...
		for(
			TiXmlElement* ptch =
			layer->FirstChildElement( "t2_patch" );
			ptch;
			ptch = ptch->NextSiblingElement() ) {

			string	dir;

			if( !IsDir( dir, ptch ) )
				fprintf( flog, "%d\tD\t%s\n", z, dir.c_str() );
			else if( !IsImg( dir, ptch ) )
				fprintf( flog, "%d\tT\t%s\n", z, dir.c_str() );
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


