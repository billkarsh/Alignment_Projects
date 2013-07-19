//
// Change all file_path to be of form 'Images/Znn/imagename'.
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
	char	*infile;

public:
	CArgs_xml()
	{
		infile	= NULL;
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

	flog = FileOpenOrDie( "XMLLocalPaths.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 2 ) {
		printf(
		"Usage: XMLLocalPaths <xml-file1>.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* UpdateLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateLayer( TiXmlElement* layer )
{
	int				z = atoi( layer->Attribute( "z" ) );
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		char	name[256], buf[2048];

		sprintf( name, "%s",
			FileNamePtr( ptch->Attribute( "file_path" ) ) );

		sprintf( buf, "Images/Z%d/%s", z, name );
		ptch->SetAttribute( "file_path", buf );
	}
}


/* --------------------------------------------------------------- */
/* Update -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update()
{
	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

	for( ; layer; layer = layer->NextSiblingElement() )
		UpdateLayer( layer );

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

/* ------------- */
/* Write new xml */
/* ------------- */

	Update();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


