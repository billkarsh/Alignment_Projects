//
// Reformat either:
// {-x=xml file, -r=Rick file, -d-idb}.
//
// New xml files have 'title' attributes like this:
// "z.id:rgn" or,
// "z.id:rgn_col.row.cam" (if data present).
//

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"PipeFiles.h"


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

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID	re_id;

public:
	char	*inpath;
	int		cmd,		// {'x','r','d'}
			zmin,
			zmax;

public:
	CArgs_xml()
	{
		inpath	= NULL;
		cmd		= 0;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );

	int IDFromPatch( TiXmlElement* p );
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

	flog = FileOpenOrDie( "Reformat.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	char	*pat;

	re_id.Set( "_N_" );

	if( argc < 6 ) {
		printf(
		"Usage: reformat path <-x,-r,-d> -p=_Nex.mrc -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			inpath = argv[i];
		else if( IsArg( "-x", argv[i] ) )
			cmd = 'x';
		else if( IsArg( "-r", argv[i] ) )
			cmd = 'r';
		else if( IsArg( "-d", argv[i] ) )
			cmd = 'd';
		else if( GetArgStr( pat, "-p=", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromPatch( TiXmlElement* p )
{
	const char	*name = p->Attribute( "title" );
	int			id;

	if( !re_id.Decode( id, name ) ) {
		fprintf( flog, "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* UpdateXMLLayer ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void UpdateXMLLayer( TiXmlElement* layer, int z )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		char	title[128];
		int		id = gArgs.IDFromPatch( p );

		const char	*c, *name = p->Attribute( "title" );

		if( c = strstr( name, "col" ) ) {

			int	col = -1, row = -1, cam = -1;
			sscanf( c, "col%d_row%d_cam%d", &col, &row, &cam );

			sprintf( title, "%d.%d:1_%d.%d.%d",
				z, id, col, row, cam );
		}
		else
			sprintf( title, "%d.%d:1", z, id );

		p->SetAttribute( "title", title );
	}
}

/* --------------------------------------------------------------- */
/* UpdateXML ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateXML()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.inpath, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		UpdateXMLLayer( layer, z );
	}

/* ---- */
/* Save */
/* ---- */

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

	switch( gArgs.cmd ) {

		case 'x':
			UpdateXML();
		break;
	}

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


