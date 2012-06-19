//
// Given xml file pointing at apig-named images of form:
// A7_id_chn.tif, set all image names to given channel.
//
// Channel is any string having no white space or '_'.
//
// Accepts layer range [zmin, zmax].
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
/* CArgs_rgbm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_rgbm {

public:
	char	*infile,
			*chn;
	int		zmin, zmax;

public:
	CArgs_rgbm()
	{
		infile	= NULL;
		chn		= NULL;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_rgbm	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_rgbm::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "SetApigChan.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: SetApigChan <xml-file> <chn> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile )
				infile = argv[i];
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

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* EditTag ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// The channel part of path/A7_id_chn.xxx.tif is between
// the first dot and the previous underscore.
//
static void EditTag( TiXmlElement* ptch, const char *tag )
{
	char		buf[2048];
	const char	*n = ptch->Attribute( tag );
	const char	*d = strchr( n, '.' );

	if( !d )
		d = n + strlen( n );

	const char	*u = d;

	while( *--u != '_' )
		;

	sprintf( buf, "%.*s_%s%s", u - n, n, gArgs.chn, d );

	ptch->SetAttribute( tag, buf );
}

/* --------------------------------------------------------------- */
/* UpdateTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTiles( TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		EditTag( ptch, "title" );
		EditTag( ptch, "file_path" );
	}
}

/* --------------------------------------------------------------- */
/* WriteXML ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteXML()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---------- */
/* Fix layers */
/* ---------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		UpdateTiles( layer );
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

	WriteXML();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


