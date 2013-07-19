//
// > MRCSD1Lyr file.xml -z=121
//
// Given xml file of 16-bit mrc images, create text file named
// sd_zzz.txt, each line of which lists: z, tileID, image stddev.
//
// Memory usage profiling can be enabled (last line) for testing.
//

#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"mrc.h"
#include	"Memory.h"


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
	char	*xmlfile;
	int		z;

public:
	CArgs_xml()
	{
		xmlfile	= NULL;
		z		= 0;
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
// parse command line args

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
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

// start log

	char	buf[256];

	sprintf( buf, "sd_%d.txt", z );
	flog = FileOpenOrDie( buf, "w" );
}

/* --------------------------------------------------------------- */
/* GetSD --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static int GetSD( TiXmlElement* p )
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
/* GetTileSDs ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetTileSDs( TiXmlElement* layer, int z )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		fprintf( flog, "%d\t%d\t%d\n",
		z, IDFromPatch( p ), GetSD( p ) );
	}
}

/* --------------------------------------------------------------- */
/* Process ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Process()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.xmlfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------- */
/* Do layer */
/* -------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.z )
			break;

		if( z < gArgs.z )
			continue;

		GetTileSDs( layer, z );
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


