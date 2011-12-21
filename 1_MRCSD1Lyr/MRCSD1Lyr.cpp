//
// > MRCSD1Lyr file.xml -z=121 -p_Nex.mrc
//
// Given xml file of 16-bit mrc images, create text file named
// sd_zzz.txt, each line of which lists: z, tileID, image stddev.
//
// Memory usage profiling can be enabled (last line) for testing.
//

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"mrc.h"
#include	"Memory.h"

#include	"tinyxml.h"


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
	char	*xmlfile;
	int		z;

public:
	CArgs_xml()
	{
		xmlfile	= NULL;
		z		= 0;
	};

	void SetCmdLine( int argc, char* argv[] );

	int IDFromPatch( TiXmlElement *p );
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
// parse command line args

	char	*pat;

	re_id.Set( "_Nex.mrc" );

	if( argc < 3 ) {
		printf(
		"Usage: MRCSD1Lyr <xml-file> -z= [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			xmlfile = argv[i];
		else if( GetArg( &z, "-z=%d", argv[i] ) )
			;
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
	}

	re_id.Compile( stdout );

// start log

	char	buf[256];

	sprintf( buf, "sd_%d.txt", z );
	flog = FileOpenOrDie( buf, "w" );
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromPatch( TiXmlElement *p )
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
/* GetSD --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static int GetSD( TiXmlElement *p )
{
	vector<uint16*>	vras;
	uint32			w, h;

	if( 1 > ReadRawMRCFile( vras, p->Attribute( "file_path" ),
				w, h, NULL ) ) {

		return 0;
	}

	const uint16*	V = &vras[0][0];
	double			sd, sm = 0.0, sm2 = 0.0;
	int				n = w * h;

	for( int i = 0; i < n; ++i ) {

		double	d = V[i];

		sm  += d;
		sm2 += d * d;
	}

	sd = sqrt( (sm2 - sm*sm/n) / (n - 1.0) );

	FreeMRC( vras );

	return (int)sd;
}

/* --------------------------------------------------------------- */
/* Process ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Process()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.xmlfile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog,
		"No trakEM2 node [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog,
		"No t2_layer [%s].\n", gArgs.xmlfile );
		exit( 42 );
	}

/* -------- */
/* Do layer */
/* -------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.z )
			break;

		if( z < gArgs.z )
			continue;

		for(
			TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
			p;
			p = p->NextSiblingElement() ) {

			fprintf( flog, "%d\t%d\t%d\n",
			z, gArgs.IDFromPatch( p ), GetSD( p ) );
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

/* ------- */
/* Process */
/* ------- */

	Process();

/* ---- */
/* Done */
/* ---- */

exit:
//	VMStats( flog );
	fclose( flog );

	return 0;
}


