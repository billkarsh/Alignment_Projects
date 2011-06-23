//
// Using GRScales.txt, update the corresponding xml layers
// to have those scales and given channel.
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

class Scale {

public:
	int	z, smin, smax;

public:
	Scale( int _z, int _min, int _max )
		{z = _z; smin = _min; smax = _max;};

	bool operator < (const Scale &rhs) const
		{return z < rhs.z;};
};

/* --------------------------------------------------------------- */
/* CArgs_gray ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_gray {

public:
	char	*xmlfile;
	int		chn;

public:
	CArgs_gray()
	{
		xmlfile	= NULL;
		chn		= -1;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_gray	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_gray::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "GraRanXML.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: GraRanXML <xml-file> <scale-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			xmlfile = argv[i];
		else if( GetArg( &chn, "-chn=%d", argv[i] ) )
			;
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* ReadScales ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadScales( vector<Scale> &vs )
{
	int	z0, nz;

	FILE	*f = FileOpenOrDie( "GRTemp/range.txt", "r", flog );
	fscanf( f, "%d\t%d", &z0, &nz );
	fclose( f );

	for( int i = 0; i < nz; ++i ) {

		char	buf[2048];
		sprintf( buf, "GRTemp/z_%d.txt", z0+i );
		FILE	*f = fopen( buf, "r" );

		if( f ) {

			int	smin, smax;

			fscanf( f, "%d\t%d", &smin, &smax );
			fclose( f );

			vs.push_back( Scale( z0+i, smin, smax ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* EditTitleAndPath ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void EditTitleAndPath( TiXmlElement* ptch )
{
	char	buf[2048];
	int		len;

// title first

	len = sprintf( buf, "%s", ptch->Attribute( "title" ) );
	buf[len - 5] = '0' + gArgs.chn;
	ptch->SetAttribute( "title", buf );

// now path

	len = sprintf( buf, "%s", ptch->Attribute( "file_path" ) );
	buf[len - 5] = '0' + gArgs.chn;
	ptch->SetAttribute( "file_path", buf );
}

/* --------------------------------------------------------------- */
/* WriteXML ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteXML( const vector<Scale> &vs )
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.xmlfile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.xmlfile );
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

/* ---------- */
/* Fix layers */
/* ---------- */

	int	nz = vs.size();

	// for each layer we wish to fix...
	for( int iz = 0; iz < nz; ++iz ) {

		// advance xml to that layer
		while( atoi( layer->Attribute( "z" ) ) < vs[iz].z )
			layer = layer->NextSiblingElement();

		if( !layer )
			break;

		// for each tile in this layer...
		for(
			TiXmlElement* ptch =
			layer->FirstChildElement( "t2_patch" );
			ptch;
			ptch = ptch->NextSiblingElement() ) {

			// edit names
			if( gArgs.chn >= 0 )
				EditTitleAndPath( ptch );

			// edit min and max
			ptch->SetAttribute( "min", vs[iz].smin );
			ptch->SetAttribute( "max", vs[iz].smax );
		}
	}

/* ---- */
/* Save */
/* ---- */

	doc.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.xmlfile, "xmltmp.txt" );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Scale>	vs;

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ------------------------ */
/* Read and sort scale file */
/* ------------------------ */

	ReadScales( vs );

	if( !vs.size() )
		goto exit;

//	sort( vs.begin(), vs.end() );

/* ------------- */
/* Write new xml */
/* ------------- */

	WriteXML( vs );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


