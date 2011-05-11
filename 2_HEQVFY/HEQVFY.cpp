//
// List missing HEQ files
//
//

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"Disk.h"
#include	"File.h"

#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_heq ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_heq {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images.
	CRegexID	re_id;

public:
	char	*infile,
			*tag;
	int		zmin, zmax, chn;

public:
	CArgs_heq()
	{
		infile	= NULL;
		tag		= NULL;
		zmin	= 0;
		zmax	= 32768;
		chn		= -1;
	};

	void SetCmdLine( int argc, char* argv[] );

	int IDFromPatch( TiXmlElement *p );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_heq	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_heq::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "HEQVFY.log", "w" );

// parse command line args

	char	*pat;

	re_id.Set( "_N_" );

	if( argc < 6 ) {
		printf( "Usage: HEQVFY <xml-file> <tag>"
		" -zmin=i -zmax=j -chn=c.\n" );
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
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &chn, "-chn=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

// header

	fprintf( flog, "\n\nMissing HEQ---\nZ\tID\n" );
	fflush( flog );
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_heq::IDFromPatch( TiXmlElement *p )
{
	const char	*name = p->Attribute( "title" );
	int			id;

	if( !re_id.Decode( id, name ) ) {
		printf( "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* OutName ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char *OutName( char *buf, TiXmlElement *p )
{
	char		tag[32];
	const char	*n = p->Attribute( "file_path" );
	char		*s = strrchr( n, '/' ),
				*c = strrchr( s, '_' );
	int			len = sprintf( tag, "_%s", gArgs.tag );

	if( !strncmp( tag, s - len, len ) ) {

		// dir already has tag - look there

		sprintf( buf, "%.*s_%d.tif", c - n, n, gArgs.chn );
	}
	else {

		// compose tagged name - look there

		sprintf( buf, "%.*s%s%.*s_%d.tif",
			s - n, n, tag, c - s, s, gArgs.chn );
	}

	return buf;
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	layer = hDoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

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

		/* ------------------------------ */
		/* For each patch (tile) in layer */
		/* ------------------------------ */

		for(
			TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
			p;
			p = p->NextSiblingElement() ) {

			char	buf[2048];
			int		id = gArgs.IDFromPatch( p );

			if( !DskExists( OutName( buf, p ) ) )
				fprintf( flog, "%d\t%d\n", z, id );
		}
	}
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

/* ---------------- */
/* Read source file */
/* ---------------- */

	ParseTrakEM2();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


