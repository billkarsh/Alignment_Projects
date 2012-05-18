//
// Concatenate TrakEM2 file2 to file1.
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
	char	*infile1,
			*infile2;

public:
	CArgs_xml()
	{
		infile1	= NULL;
		infile2	= NULL;
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

	flog = FileOpenOrDie( "XMLAppend.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: XMLAppend <xml-file1> <xml-file2> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );
	}

	infile1 = argv[1];
	infile2 = argv[2];

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Append -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Append()
{
/* -------------- */
/* Load documents */
/* -------------- */

	TiXmlDocument	doc1( gArgs.infile1 ),
					doc2( gArgs.infile2 );

	if( !doc1.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile1 );
		exit( 42 );
	}

	if( !doc2.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile2 );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc1( &doc1 ),
					hDoc2( &doc2 );
	TiXmlNode*		lyrset1;
	TiXmlElement*	layer1;
	TiXmlElement*	layer2;

	if( !doc1.FirstChild() || !doc2.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	lyrset1 = hDoc1.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.ToNode();

	layer2 = hDoc2.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !lyrset1 || !layer2 ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

/* ------ */
/* Append */
/* ------ */

	layer1 = lyrset1->LastChild( "t2_layer" )->ToElement();

	int	z		= atoi( layer1->Attribute( "z" ) );	// last z
	int	nextoid	= NextOID( hDoc1 );

	do {
		layer2->SetAttribute( "z", ++z );
		nextoid = SetOID( layer2, nextoid );
		lyrset1->InsertEndChild( *layer2 );
	} while( layer2 = layer2->NextSiblingElement() );

/* ---- */
/* Save */
/* ---- */

	doc1.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile1, "xmltmp.txt" );
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

	Append();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


