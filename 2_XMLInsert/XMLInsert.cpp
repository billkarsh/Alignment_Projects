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
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Insert -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Insert()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml1( gArgs.infile1, flog );
	TiXmlNode*		lyrset1	= xml1.GetLayerset();
	TiXmlElement*	layer1	= xml1.GetLastLayer();

	XML_TKEM		xml2( gArgs.infile2, flog );
	TiXmlElement*	layer2	= xml2.GetFirstLayer();

/* --------------------------- */
/* If 'after' all, then append */
/* --------------------------- */

	int	z		= atoi( layer1->Attribute( "z" ) ),	// last z
		nextZ	= gArgs.after + 1,
		nextoid	= xml1.NextOID();

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

	Insert();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


