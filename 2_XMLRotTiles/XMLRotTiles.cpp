//
// Three modes of operation:
//
// (0)
// XMLRotTiles file.xml -mode=0 -id1=i -id2=j -zmin=z1 -zmax=z2
// Using spec tiles (id1, id2) which are left-to-right across
// horz line, define global layer angle, and, generate table:
//
// Z  Layer-ang  id-tile-ang  delta-ang  *-if-bigger-than-2deg
//
// (1)
// XMLRotTiles file.xml -mode=1 -x=4096 -y=4096 -z=i -degcw=10
// Rotate tiles in spec layer by spec clockwise angle. Each tile
// is dim (x,y).
//
// (2)
// XMLRotTiles file.xml -mode=2 -id1=i -id2=j -x=4096 -y=4096 \
//	-zmin=z1 -zmax=z2 -tdeg=0
// Using same calcs as mode (0), actually rotate the layer if delta
// exceeds spec thresh tdeg. Each tile is dim (x,y).
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
	double		degcw, tdeg;
	char		*infile;
	int			mode, id1, id2, x, y, z, zmin, zmax;

public:
	CArgs_xml()
	{
		degcw	= 0.0;
		tdeg	= 0.0;
		infile	= NULL;
		mode	= 0;
		x		= 1376;
		y		= 1040;
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

	flog = FileOpenOrDie( "XMLRotTiles.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf( "Usage: XMLRotTiles <xml-file> -z=i -degcw=f\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &mode, "-mode=%d", argv[i] ) )
			;
		else if( GetArg( &id1, "-id1=%d", argv[i] ) )
			;
		else if( GetArg( &id2, "-id2=%d", argv[i] ) )
			;
		else if( GetArg( &x, "-x=%d", argv[i] ) )
			;
		else if( GetArg( &y, "-y=%d", argv[i] ) )
			;
		else if( GetArg( &z, "-z=%d", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &degcw, "-degcw=%lf", argv[i] ) )
			;
		else if( GetArg( &tdeg, "-tdeg=%lf", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* RotateLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RotateLayer( TiXmlElement* layer, double degcw )
{
/* --------------- */
/* Create rotation */
/* --------------- */

	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );
	TAffine			TR;

	TR.SetCWRot( degcw, Point( gArgs.x/2, gArgs.y/2 ) );

	for( ; p; p = p->NextSiblingElement() ) {

		TAffine	T0, T;

		T0.ScanTrackEM2( p->Attribute( "transform" ) );
		T = T0 * TR;
		XMLSetTFVals( p, T.t );
	}
}

/* --------------------------------------------------------------- */
/* GetTheTwoTAffines --------------------------------------------- */
/* --------------------------------------------------------------- */

static bool GetTheTwoTAffines(
	TAffine			&T1,
	TAffine			&T2,
	TiXmlElement*	layer )
{
	TiXmlElement*	p		= layer->FirstChildElement( "t2_patch" );
	int				ngot	= 0;

	for( ; ngot < 2 && p; p = p->NextSiblingElement() ) {

		int	id = IDFromPatch( p );

		if( id == gArgs.id1 ) {
			T1.ScanTrackEM2( p->Attribute( "transform" ) );
			++ngot;
		}
		else if( id == gArgs.id2 ) {
			T2.ScanTrackEM2( p->Attribute( "transform" ) );
			++ngot;
		}
	}

	return (ngot == 2);
}

/* --------------------------------------------------------------- */
/* Report -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Report()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---------- */
/* All layers */
/* ---------- */

	fprintf( flog, "Z\tGlobal\tTileAve\tdegCW\tBig\n" );

	for( ; layer; layer = layer->NextSiblingElement() ) {

		TAffine	T1, T2;
		int		z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		if( !GetTheTwoTAffines( T1, T2, layer ) ) {
			fprintf( flog, "%d\tMissing ref tile\n", z );
			continue;
		}

		double	base, tile, cw;
		char	big = ' ';

		base = -(T2.t[5] - T1.t[5]) / (T2.t[2] - T1.t[2]);
		base = 180/PI * atan( base );

		tile = T1.GetRadians() + T2.GetRadians();
		tile = 180/PI * tile / 2;

		cw = -(base + tile);

		if( fabs( cw ) > 2 )
			big = '*';

		fprintf( flog, "%d\t%.2f\t%.2f\t%.2f\t%c\n",
			z, base, tile, cw, big );
	}
}

/* --------------------------------------------------------------- */
/* RotateAuto ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RotateAuto()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---------- */
/* All layers */
/* ---------- */

	fprintf( flog, "Z\tGlobal\tTileAve\tdegCW\tBig\n" );

	for( ; layer; layer = layer->NextSiblingElement() ) {

		TAffine	T1, T2;
		int		z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		if( !GetTheTwoTAffines( T1, T2, layer ) ) {
			fprintf( flog, "%d\tMissing ref tile\n", z );
			continue;
		}

		double	base, tile, cw;

		base = -(T2.t[5] - T1.t[5]) / (T2.t[2] - T1.t[2]);
		base = 180/PI * atan( base );

		tile = T1.GetRadians() + T2.GetRadians();
		tile = 180/PI * tile / 2;

		cw = -(base + tile);

		if( fabs( cw ) >= gArgs.tdeg )
			RotateLayer( layer, cw );
	}

/* ---- */
/* Save */
/* ---- */

	xml.Save( "xmltmp.txt", true );
}

/* --------------------------------------------------------------- */
/* Rotate1Layer -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Rotate1Layer()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* Find our layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z != gArgs.z )
			continue;

		RotateLayer( layer, gArgs.degcw );
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

	if( gArgs.mode == 1 )
		Rotate1Layer();
	else if( gArgs.mode == 2 )
		RotateAuto();
	else
		Report();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


