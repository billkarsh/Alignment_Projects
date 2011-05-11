//
// Keep only TrakEM2 layers in given z range.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"


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
	vector<char*>	infile;
	int				zmin, zmax;

public:
	CArgs_xml()
	{
		zmin	= 0;
		zmax	= 32768;
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

	flog = FileOpenOrDie( "XMLExtract.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: XMLExtract <xml-file1> <xml-file2>"
		" [<xml-filej> ...] -zmin=i -zmax=j [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile.push_back( argv[i] );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Extract ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Extract( int i )
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile[i] );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile[i] );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog,
		"No trakEM2 node [%s].\n", gArgs.infile[i] );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog,
		"No t2_layer [%s].\n", gArgs.infile[i] );
		exit( 42 );
	}

/* ------------------------- */
/* Kill layers outside range */
/* ------------------------- */

	TiXmlNode		*lyrset = layer->Parent();
	TiXmlElement	*next;

	for( ; layer; layer = next ) {

		// next layer0 before deleting anything
		next = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( z < gArgs.zmin || z > gArgs.zmax )
			lyrset->RemoveChild( layer );
	}

/* ---- */
/* Save */
/* ---- */

	doc.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile[i], "xmltmp.txt" );
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

	for( int i = 0; i < gArgs.infile.size(); ++i )
		Extract( i );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


