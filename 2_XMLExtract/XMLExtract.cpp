//
// Keep only TrakEM2 layers in given z range and, if using
// -lrbt= or -xyr= options, in given XY-box.
//
// -lrbt= specifies left,right,bottom,top of a bbox.
// -xyr=  specifies a box using Xcenter,Ycenter,radius.
//

#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"TAffine.h"


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
	vector<double>	lrbt;		// if not empty, forced bbox
	char			*infile;
	int				zmin, zmax;

public:
	CArgs_xml()
	{
		infile	= NULL;
		zmin	= 0;
		zmax	= 32768;
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

	flog = FileOpenOrDie( "XMLExtract.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf( "Usage: XMLExtract <xml-file1> -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	vector<double>	xyr;

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArgList( lrbt, "-lrbt=", argv[i] ) ) {

			if( 4 != lrbt.size() ) {
				printf( "Bad format in -lrbt [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else if( GetArgList( xyr, "-xyr=", argv[i] ) ) {

			if( 3 != xyr.size() ) {
				printf( "Bad format in -xyr [%s].\n", argv[i] );
				exit( 42 );
			}

			lrbt.resize( 4 );
			lrbt[0] = xyr[0] - xyr[2];
			lrbt[1] = xyr[0] + xyr[2];
			lrbt[2] = xyr[1] - xyr[2];
			lrbt[3] = xyr[1] + xyr[2];
		}
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* TrimTiles ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void TrimTiles( TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );
	TiXmlElement*	pnext;

	for( ; ptch; ptch = pnext ) {

		pnext = ptch->NextSiblingElement();

		vector<Point>	cnr( 4 );
		TAffine			T;
		DBox			B;
		int				w, h;

		w = atoi( ptch->Attribute( "width" ) );
		h = atoi( ptch->Attribute( "height" ) );
		T.ScanTrackEM2( ptch->Attribute( "transform" ) );

		B.L = BIGD, B.R = -BIGD,
		B.B = BIGD, B.T = -BIGD;

		cnr[0] = Point( 0.0, 0.0 );
		cnr[1] = Point( w-1, 0.0 );
		cnr[2] = Point( w-1, h-1 );
		cnr[3] = Point( 0.0, h-1 );

		T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			B.L = fmin( B.L, cnr[k].x );
			B.R = fmax( B.R, cnr[k].x );
			B.B = fmin( B.B, cnr[k].y );
			B.T = fmax( B.T, cnr[k].y );
		}

		if( B.R <= gArgs.lrbt[0] || B.L >= gArgs.lrbt[1] ||
			B.T <= gArgs.lrbt[2] || B.B >= gArgs.lrbt[3] ) {

			layer->RemoveChild( ptch );
		}
	}
}

/* --------------------------------------------------------------- */
/* Extract ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Extract()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlNode*		lyrset	= xml.GetLayerset();
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ------------------------- */
/* Kill layers outside range */
/* ------------------------- */

	TiXmlElement*	next;

	for( ; layer; layer = next ) {

		/* ---------------- */
		/* Layer in bounds? */
		/* ---------------- */

		// next layer0 before deleting anything
		next = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( z < gArgs.zmin || z > gArgs.zmax )
			lyrset->RemoveChild( layer );

		/* --------------- */
		/* Tile in bounds? */
		/* --------------- */

		if( !gArgs.lrbt.size() )
			continue;

		TrimTiles( layer );
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

	Extract();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


