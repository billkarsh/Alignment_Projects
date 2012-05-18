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
/* ---- */
/* Open */
/* ---- */

	int	nf = gArgs.infile.size();

	vector<XML_TKEM>		xml;
	vector<TiXmlElement*>	layer( nf );

	for( int i = 0; i < nf; ++i ) {

		xml.push_back( XML_TKEM( gArgs.infile[i], flog ) );
		layer[i] = xml[i].GetFirstLayer();
	}

/* --------------------------------------------- */
/* Interleave, adopting layer structure of file0 */
/* --------------------------------------------- */

	TiXmlNode*		lyrset0 = layer[0]->Parent();
	TiXmlElement*	next0;
	int				nextoid	= xml[0].NextOID();

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
	xml[0].Save( "xmltmp.txt", true );
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


