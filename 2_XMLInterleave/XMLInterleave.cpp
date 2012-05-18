//
// Take N TrakEM2 files...
//
// Using z layers of file0, insert matching layers from
// other files into file0.
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

public:
	CArgs_xml()
	{
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

	flog = FileOpenOrDie( "XMLInterleave.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: XMLInterleave <xml-file1> <xml-file2>"
		" [<xml-filej> ...] [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile.push_back( argv[i] );
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Interleave ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Interleave()
{
/* -------------- */
/* Load documents */
/* -------------- */

	int	nf = gArgs.infile.size();

	vector<TiXmlDocument>	doc;

	for( int i = 0; i < nf; ++i ) {

		doc.push_back( TiXmlDocument( gArgs.infile[i] ) );

		if( !doc[i].LoadFile() ) {
			fprintf( flog,
			"Could not open XML file [%s].\n", gArgs.infile[i] );
			exit( 42 );
		}
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	vector<TiXmlHandle>		hdoc;
	vector<TiXmlElement*>	layer( nf );

	for( int i = 0; i < nf; ++i ) {

		hdoc.push_back( TiXmlHandle( &doc[i] ) );

		if( !doc[i].FirstChild() ) {
			fprintf( flog,
			"No trakEM2 node [%s].\n", gArgs.infile[i] );
			exit( 42 );
		}

		layer[i] = hdoc[i].FirstChild( "trakem2" )
					.FirstChild( "t2_layer_set" )
					.FirstChild( "t2_layer" )
					.ToElement();

		if( !layer[i] ) {
			fprintf( flog,
			"No t2_layer [%s].\n", gArgs.infile[i] );
			exit( 42 );
		}
	}

/* --------------------------------------------- */
/* Interleave, adopting layer structure of file0 */
/* --------------------------------------------- */

	TiXmlNode*		lyrset0 = layer[0]->Parent();
	TiXmlElement*	next0;
	int				nextoid	= NextOID( hdoc[0] );

	for( ;; ) {

		// next layer0 before adding anything
		next0 = layer[0]->NextSiblingElement();

		// set z and add with locally reversed order
		for( int i = nf - 1; i > 0; --i ) {

			layer[i]->SetAttribute( "z", layer[0]->Attribute( "z" ) );
			nextoid = SetOID( layer[i], nextoid );
			lyrset0->InsertBeforeChild( layer[0], *layer[i] );
		}

		//advance layers
		if( !(layer[0] = next0) )
			goto save;

		for( int i = 1; i < nf; ++i ) {

			if( !(layer[i] = layer[i]->NextSiblingElement()) )
				goto save;
		}
	}

/* ---- */
/* Save */
/* ---- */

save:
	doc[0].SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile[0], "xmltmp.txt" );
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

	Interleave();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


