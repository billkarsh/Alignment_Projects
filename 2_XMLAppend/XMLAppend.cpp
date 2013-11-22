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
	CArgs_xml() : infile1(NULL), infile2(NULL) {};

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
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml1( gArgs.infile1, flog );
	TiXmlNode*		lyrset1	= xml1.GetLayerset();
	TiXmlElement*	layer1	= xml1.GetLastLayer();

	XML_TKEM		xml2( gArgs.infile2, flog );
	TiXmlElement*	layer2	= xml2.GetFirstLayer();

/* ------ */
/* Append */
/* ------ */

	int	z		= atoi( layer1->Attribute( "z" ) );	// last z
	int	nextoid	= xml1.NextOID();

	do {
		layer2->SetAttribute( "z", ++z );
		nextoid = SetOID( layer2, nextoid );
		lyrset1->InsertEndChild( *layer2 );
	} while( layer2 = layer2->NextSiblingElement() );

/* ---- */
/* Save */
/* ---- */

	xml1.Save( "xmltmp.txt", true );
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


