//
// Use montage/TAffineTable files from given working alignment dir
// to update given xml file.
//
// Default is to keep layer orientation same as xml file,
// effectively just improving the accuracy of the xml file.
//
// However, -force option replaces xml transforms with the
// lsq montage solutions.
//

#include	"Cmdline.h"
#include	"Disk.h"
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

public:
	char	*xmlfile,
			*tempdir;
	int		zmin,
			zmax;
	bool	force;

public:
	CArgs_xml()
	{
		xmlfile	= NULL;
		tempdir	= NULL;
		zmin	= 0;
		zmax	= 32768;
		force	= false;
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

	flog = FileOpenOrDie( "updatemontages.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: updatemontages <xml-file> <aln-dir> -zmin=i, -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !xmlfile )
				xmlfile = argv[i];
			else
				tempdir = argv[i];
		}
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( IsArg( "-force", argv[i] ) )
			force = true;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* GetTable ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool GetTable( map<MZIDR,TAffine> &M, int z )
{
	char	file[2048];

	sprintf( file, "%s/%d/montage/TAffineTable.txt",
		gArgs.tempdir, z );

	if( !DskExists( file ) ) {
		fprintf( flog, "Layer %d: No TAffineTable.\n", z );
		return false;
	}

	LoadTAffineTbl_RngZ( M, z, z, file, flog );
	return true;
}

/* --------------------------------------------------------------- */
/* CalcTupdt ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool CalcTupdt(
	TAffine				&Tupdt,
	TiXmlElement*		layer,
	int					z,
	map<MZIDR,TAffine>	&M )
{
	MZIDR			key;
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	key.z	= z;
	key.rgn	= 1;

	for( ; p; p = p->NextSiblingElement() ) {

		key.id = IDFromPatch( p );

		map<MZIDR,TAffine>::iterator	it = M.find( key );

		if( it != M.end() ) {

			TAffine	T0, TR;

			// get previous Transform
			T0.ScanTrackEM2( p->Attribute( "transform" ) );

			// set delta rotation part of Tupdt
			Tupdt.NUSetRot(
				T0.GetRadians() - it->second.GetRadians() );

			// set delta translation part of Tupdt
			TR = Tupdt * it->second;
			Tupdt.SetXY(
				T0.t[2] - TR.t[2],
				T0.t[5] - TR.t[5] );

			return true;
		}
	}

	fprintf( flog, "Layer %d: No tile matched.\n", z );
	return false;
}

/* --------------------------------------------------------------- */
/* Apply --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Apply(
	TiXmlElement*		layer,
	int					z,
	const TAffine		&Tupdt,
	map<MZIDR,TAffine>	&M )
{
	MZIDR			key;
	TiXmlElement*	next;
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	key.z	= z;
	key.rgn	= 1;

	for( ; p; p = next ) {

		next = p->NextSiblingElement();

		key.id = IDFromPatch( p );

		map<MZIDR,TAffine>::iterator	it = M.find( key );

		if( it != M.end() ) {

			TAffine	T = Tupdt * it->second;
			XMLSetTFVals( p, T.t );
		}
		else {
			layer->RemoveChild( p );
			fprintf( flog, "Layer %d: Tile %d dropped.\n", z, key.id );
		}
	}
}

/* --------------------------------------------------------------- */
/* Force --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Force(
	TiXmlElement*		layer,
	int					z,
	map<MZIDR,TAffine>	&M )
{
	MZIDR			key;
	TiXmlElement*	next;
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	key.z	= z;
	key.rgn	= 1;

	for( ; p; p = next ) {

		next = p->NextSiblingElement();

		key.id = IDFromPatch( p );

		map<MZIDR,TAffine>::iterator	it = M.find( key );

		if( it != M.end() )
			XMLSetTFVals( p, it->second.t );
		else {
			layer->RemoveChild( p );
			fprintf( flog, "Layer %d: Tile %d dropped.\n", z, key.id );
		}
	}
}

/* --------------------------------------------------------------- */
/* DoLayer ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DoLayer( TiXmlElement* layer, int z )
{
	map<MZIDR,TAffine>	M;
	TAffine				Tupdt;

	if( GetTable( M, z ) && CalcTupdt( Tupdt, layer, z, M ) )
		Apply( layer, z, Tupdt, M );
}

/* --------------------------------------------------------------- */
/* DoLayerForce -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void DoLayerForce( TiXmlElement* layer, int z )
{
	map<MZIDR,TAffine>	M;

	if( GetTable( M, z ) )
		Force( layer, z, M );
}

/* --------------------------------------------------------------- */
/* Update -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.xmlfile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ------------------------ */
/* Copy matching transforms */
/* ------------------------ */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		if( gArgs.force )
			DoLayerForce( layer, z );
		else
			DoLayer( layer, z );
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

	Update();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


