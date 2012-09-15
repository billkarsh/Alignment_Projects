//
// Transfer edited rough stack back to all tiles.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"CTileSet.h"
#include	"TrakEM2_UTL.h"
#include	"../1_AlignMontages2/ScapeMeta.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_alnmon -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_alnmon {

public:
	char	*xmlfile,
			*pat;
	char	outdir[2048];
	int		zmin,
			zmax,
			xml_type,
			xml_min,
			xml_max;

public:
	CArgs_alnmon()
	{
		xmlfile		= NULL;
		pat			= "/N";
		zmin		= 0;
		zmax		= 32768;
		xml_type	= 0;
		xml_min		= 0;
		xml_max		= 0;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_alnmon	gArgs;
static CTileSet		TS;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_alnmon::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "alignmontages3.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Assemble stack: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: alignmontages3 <xml-file> -d. -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		char	*_outdir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			xmlfile = argv[i];
		else if( GetArgStr( pat, "-p", argv[i] ) )
			;
		else if( GetArgStr( _outdir, "-d", argv[i] ) )
			DskAbsPath( outdir, sizeof(outdir), _outdir, flog );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			;
		else if( GetArg( &xml_min, "-xmlmin=%d", argv[i] ) )
			;
		else if( GetArg( &xml_max, "-xmlmax=%d", argv[i] ) )
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
/* LoadTForms ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LoadTForms( vector<TForm> &vT )
{
	char	path[2048];
	sprintf( path, "%s/LowResMons.xml", gArgs.outdir );

	XML_TKEM		xml( path, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();
	int				nz = 0;

	for( ; layer; layer = layer->NextSiblingElement() ) {

		TiXmlElement* p = layer->FirstChildElement( "t2_patch" );

		vT[nz++].ScanTrackEM2( p->Attribute( "transform" ) );
	}
}

/* --------------------------------------------------------------- */
/* UpdateTForms -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTForms()
{
// Get log data

	vector<CLog>	vL;
	int				nL = ReadLogs( vL, gArgs.outdir,
							gArgs.zmin, gArgs.zmax );

	if( !nL )
		return;

// Get scaled-down image -> image TForms

	vector<TForm>	vTs( nL );

	LoadTForms( vTs );

// Build whole montage TForms

	vector<TForm>	vTm( nL );

	for( int ia = 1; ia < nL; ++ia ) {

		const CScapeMeta	&Ma = vL[ia].M;
		const CScapeMeta	&M0 = vL[0].M;
		TForm				R0i, Ra, t, s;

		// A-montage -> A-oriented
		Ra.NUSetRot( Ma.deg*PI/180 );

		// A-oriented -> A-image (like Scape.cpp)
		s.NUSetScl( 1.0/Ma.scl );
		MultiplyTrans( t, s, Ra );
		t.AddXY( -Ma.x0, -Ma.y0 );

		// A-image -> 0-image
		MultiplyTrans( t, vTs[ia], TForm( t ) );

		// 0-image -> 0-image oriented
		InvertTrans( R0i, vTs[0] );
		MultiplyTrans( t, R0i, TForm( t ) );

		// 0-image oriented -> 0-oriented (reverse Scape.cpp)
		t.AddXY( M0.x0, M0.y0 );
		s.NUSetScl( M0.scl );
		MultiplyTrans( t, s, TForm( t ) );

		// 0-oriented -> 0-montage
		R0i.NUSetRot( -M0.deg*PI/180 );
		MultiplyTrans( vTm[ia], R0i, t );
	}

// Apply

	int	is0, isN;

	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 ) {

		// find TForm for this z

		int	z = TS.vtil[is0].z, ib;

		for( ib = 0; ib < nL; ++ib ) {

			if( vL[ib].M.z == z )
				break;
		}

		// apply

		if( ib < nL ) {

			for( int i = is0; i < isN; ++i ) {
				TForm	&T = TS.vtil[i].T;
				MultiplyTrans( T, vTm[ib], TForm( T ) );
			}
		}

		TS.GetLayerLimits( is0 = isN, isN );
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

	TS.SetLogFile( flog );
	TS.SetDecoderPat( gArgs.pat );

/* -------------- */
/* Read dest file */
/* -------------- */

	TS.FillFromTrakEM2( gArgs.xmlfile, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SortAll_z_id();

/* ------------- */
/* Update TForms */
/* ------------- */

	UpdateTForms();

	TS.WriteTrakEM2_EZ( "AlignedTiles.xml",
		gArgs.xml_type, gArgs.xml_min, gArgs.xml_max );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


