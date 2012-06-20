//
// Using same tag and z range as for HEQLayers, update
// the corresponding xml to point to the HEQ images.
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
	char	dtag[32],
			utag[32];
	char	*infile,
			*tag;
	int		zmin, zmax,
			ltag;

public:
	CArgs_rgbm()
	{
		infile	= NULL;
		tag		= NULL;
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

	flog = FileOpenOrDie( "HEQXML.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
usage:
		printf( "Usage: HEQXML <xml-file> <tag> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile )
				infile = argv[i];
			else
				tag = argv[i];
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

	if( tag ) {

		ltag = sprintf( dtag, ".%s", tag );

		utag[0] = '_';
		strcpy( utag + 1, tag );
	}
	else
		goto usage;

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* EditTitleAndPath ---------------------------------------------- */
/* --------------------------------------------------------------- */

// If native image path has form par/dir/name.tif,
// and user supplied tag is HEQ, the new path is
// par/dir_HEQ/name.HEQ.tif.
//
static void EditTitleAndPath( TiXmlElement* ptch )
{
	char		buf[2048];
	const char	*n, *e;

// title first

	n = ptch->Attribute( "title" );
	e = FileDotPtr( n );

	sprintf( buf, "%.*s%s", e - n, n, gArgs.dtag );

	ptch->SetAttribute( "title", buf );

// now path

	n = ptch->Attribute( "file_path" );

	char	*s = strrchr( n, '/' );

	sprintf( buf,
		"%.*s%s"	// path excl / + _HEQ
		"/%s.tif",	// / + name + .tif
		s - n, n, gArgs.utag,
		ptch->Attribute( "title" ) );

	ptch->SetAttribute( "file_path", buf );
}

/* --------------------------------------------------------------- */
/* UpdateTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTiles( TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		// edit attributes
		EditTitleAndPath( ptch );
		ptch->SetAttribute( "type", 0 );
		ptch->SetAttribute( "min", 0 );
		ptch->SetAttribute( "max", 255 );
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


