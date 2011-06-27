//
// Three modes of operation:
//
// (0)
// XMLRotTiles file.xml -mode=0 -p_N_ -id1=i -id2=j -zmin=z1 -zmax=z2
// Using spec tiles (id1, id2) which are left-to-right across
// horz line, define global layer angle, and, generate table:
//
// Z  Layer-ang  id-tile-ang  delta-ang  *-if-bigger-than-2deg
//
// (1)
// XMLRotTiles file.xml -mode=1 -z=i -degcw=10
// Rotate tiles in spec layer by spec clockwise angle.
//
// (2)
// XMLRotTiles file.xml -mode=2 -p_N_ -id1=i -id2=j -zmin=z1 -zmax=z2 -tdeg=0
// Using same calcs as mode (0), actually rotate the layer if delta
// exceeds spec thresh tdeg.
//
// Uses fixed image dims (x,y) = (1376,1040).
//

#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"CTForm.h"


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
	// "/N" used for EM projects, "_N_" for APIG images.
	CRegexID	re_id;

public:
	double		degcw, tdeg;
	char		*infile;
	int			mode, id1, id2, z, zmin, zmax;

public:
	CArgs_xml()
	{
		degcw	= 0.0;
		tdeg	= 0.0;
		infile	= NULL;
		mode	= 0;
		zmin	= 0;
		zmax	= 32768;
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
// start log

	flog = FileOpenOrDie( "XMLRotTiles.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	char	*pat;

	re_id.Set( "_N_" );

	if( argc < 4 ) {
		printf( "Usage: XMLRotTiles <xml-file> -z=i -degcw=f.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &mode, "-mode=%d", argv[i] ) )
			;
		else if( GetArg( &id1, "-id1=%d", argv[i] ) )
			;
		else if( GetArg( &id2, "-id2=%d", argv[i] ) )
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
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
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
/* RotateLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RotateLayer( TiXmlElement* layer, double degcw )
{
/* --------------- */
/* Create rotation */
/* --------------- */

	TForm	TR;

	CreateCWRot( TR, degcw, Point( 1376/2, 1040/2 ) );

	for(
		TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
		p;
		p = p->NextSiblingElement() ) {

		char	buf[256];
		TForm	T0, T;

		T0.ScanTrackEM2( p->Attribute( "transform" ) );
		MultiplyTrans( T, T0, TR );

		sprintf( buf, "matrix(%f,%f,%f,%f,%f,%f)",
		T.t[0], T.t[3], T.t[1], T.t[4], T.t[2], T.t[5] );

		p->SetAttribute( "transform", buf );
	}
}

/* --------------------------------------------------------------- */
/* Report -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Report()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node [%s].\n", gArgs.infile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No t2_layer [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------- */
/* All layers */
/* ---------- */

	fprintf( flog, "Z\tGlobal\tTileAve\tdegCW\tBig\n" );

	for( ; layer; layer = layer->NextSiblingElement() ) {

		TForm	T1, T2;
		int		ngot = 0, z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		for(
			TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
			p;
			p = p->NextSiblingElement() ) {

			int	id = gArgs.IDFromPatch( p );

			if( id == gArgs.id1 ) {
				T1.ScanTrackEM2( p->Attribute( "transform" ) );
				++ngot;
			}
			else if( id == gArgs.id2 ) {
				T2.ScanTrackEM2( p->Attribute( "transform" ) );
				++ngot;
			}
		}

		if( ngot < 2 ) {
			fprintf( flog, "%d\tMissing ref tile\n", z );
			continue;
		}

		double	base, tile, cw;
		char	big = ' ';

		base = -(T2.t[5] - T1.t[5]) / (T2.t[2] - T1.t[2]);
		base = 180/PI * atan( base );

		tile = atan2( T1.t[3], T1.t[0] ) + atan2( T2.t[3], T2.t[0] );
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
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node [%s].\n", gArgs.infile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No t2_layer [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------- */
/* All layers */
/* ---------- */

	fprintf( flog, "Z\tGlobal\tTileAve\tdegCW\tBig\n" );

	for( ; layer; layer = layer->NextSiblingElement() ) {

		TForm	T1, T2;
		int		ngot = 0, z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		for(
			TiXmlElement *p = layer->FirstChildElement( "t2_patch" );
			p;
			p = p->NextSiblingElement() ) {

			int	id = gArgs.IDFromPatch( p );

			if( id == gArgs.id1 ) {
				T1.ScanTrackEM2( p->Attribute( "transform" ) );
				++ngot;
			}
			else if( id == gArgs.id2 ) {
				T2.ScanTrackEM2( p->Attribute( "transform" ) );
				++ngot;
			}
		}

		if( ngot < 2 ) {
			fprintf( flog, "%d\tMissing ref tile\n", z );
			continue;
		}

		double	base, tile, cw;

		base = -(T2.t[5] - T1.t[5]) / (T2.t[2] - T1.t[2]);
		base = 180/PI * atan( base );

		tile = atan2( T1.t[3], T1.t[0] ) + atan2( T2.t[3], T2.t[0] );
		tile = 180/PI * tile / 2;

		cw = -(base + tile);

		if( fabs( cw ) >= gArgs.tdeg )
			RotateLayer( layer, cw );
	}

/* ---- */
/* Save */
/* ---- */

	doc.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile, "xmltmp.txt" );
}

/* --------------------------------------------------------------- */
/* Rotate1Layer -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Rotate1Layer()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement	*layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node [%s].\n", gArgs.infile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No t2_layer [%s].\n", gArgs.infile );
		exit( 42 );
	}

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

	doc.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile, "xmltmp.txt" );
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

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


