//
// Insert TrakEM2 file2 into file1 immediately
// after layer 'after' of file1.
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
	int		after;

public:
	CArgs_xml()
	{
		infile1	= NULL;
		infile2	= NULL;
		after	= 0;
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

	flog = FileOpenOrDie( "XMLInsert.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: XMLInsert <xml-file1> <xml-file2> -after=i.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile1 )
				infile1 = argv[i];
			else
				infile2 = argv[i];
		}
		else if( GetArg( &after, "-after=%d", argv[i] ) )
			;
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Insert -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Insert()
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
	TiXmlNode		*lyrset1;
	TiXmlElement	*layer1, *layer2;

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

/* --------------------------- */
/* If 'after' all, then append */
/* --------------------------- */

	layer1 = lyrset1->LastChild( "t2_layer" )->ToElement();

	int	z		= atoi( layer1->Attribute( "z" ) ),	// last z
		nextZ	= gArgs.after + 1,
		nextoid	= NextOID( hDoc1 );

	if( nextZ < 0 )
		nextZ = 0;

	if( z <= gArgs.after ) {

		nextZ = z + 1;

		do {
			layer2->SetAttribute( "z", nextZ++ );
			nextoid = SetOID( layer2, nextoid );
			lyrset1->InsertEndChild( *layer2 );
		} while( layer2 = layer2->NextSiblingElement() );

		goto save;
	}

/* ---------------------------------------- */
/* Advance layer1 to that following 'after' */
/* ---------------------------------------- */

	layer1 = lyrset1->FirstChild( "t2_layer" )->ToElement();

	for(;;) {

		int	z = atoi( layer1->Attribute( "z" ) );

		if( z > gArgs.after )
			break;

		layer1 = layer1->NextSiblingElement();
	}

/* ------------------------ */
/* Insert file2 before this */
/* ------------------------ */

	do {
		layer2->SetAttribute( "z", nextZ++ );
		nextoid = SetOID( layer2, nextoid );
		lyrset1->InsertBeforeChild( layer1, *layer2 );
	} while( layer2 = layer2->NextSiblingElement() );

/* -------------------------------------- */
/* Renumber layer1's following the insert */
/* -------------------------------------- */

	for( ; layer1; layer1 = layer1->NextSiblingElement() )
		layer1->SetAttribute( "z", nextZ++ );

/* ---- */
/* Save */
/* ---- */

save:
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

	Insert();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


