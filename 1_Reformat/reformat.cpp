//
// Reformat either:
// {-x=xml file, -r=Rick file, -d-idb}.
//
// New xml files have 'title' attributes like this:
// "z.id-rgn", or,
// "z.id-rgn_col.row.cam" (if data present).
//
// New Rick files have lines like this:
// "z id 6_aff_val_by_row -999 -999 -999 path", or,
// "z id 6_aff_val_by_row col row cam path" (if data present).
//
// New idb TileToImage lines look like this:
// "id 6_aff_val_by_row -999 -999 -999 path", or,
// "id 6_aff_val_by_row col row cam path" (if data present).
//

#include	"Cmdline.h"
#include	"CRegexID.h"
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

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID	re_id;

public:
	char	*inpath;
	int		cmd,		// {'x','r','d'}
			zmin,
			zmax;

public:
	CArgs_xml()
	{
		inpath	= NULL;
		cmd		= 0;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );

	int IDFromName( const char *name );
	int IDFromPatch( TiXmlElement* p );
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

	flog = FileOpenOrDie( "Reformat.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	const char	*pat;

	re_id.Set( "_N_" );

	if( argc < 4 ) {
		printf(
		"Usage: reformat path <-x,-r,-d> -p=_Nex.mrc -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			inpath = argv[i];
		else if( IsArg( "-x", argv[i] ) )
			cmd = 'x';
		else if( IsArg( "-r", argv[i] ) )
			cmd = 'r';
		else if( IsArg( "-d", argv[i] ) )
			cmd = 'd';
		else if( GetArgStr( pat, "-p=", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
}

/* -------------------------------------------------------------- */
/* IDFromName --------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromName( const char *name )
{
	int	id;

	if( !re_id.Decode( id, FileNamePtr( name ) ) ) {
		fprintf( flog, "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromPatch( TiXmlElement* p )
{
	return IDFromName( p->Attribute( "title" ) );
}

/* --------------------------------------------------------------- */
/* UpdateXMLLayer ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void UpdateXMLLayer( TiXmlElement* layer, int z )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		char	title[128];
		int		id = gArgs.IDFromPatch( p );

		const char	*c, *n = p->Attribute( "title" );

		if( c = strstr( n, "col" ) ) {

			int	col = -1, row = -1, cam = 0;
			sscanf( c, "col%d_row%d_cam%d", &col, &row, &cam );

			sprintf( title, "%d.%d-1_%d.%d.%d",
				z, id, col, row, cam );
		}
		else
			sprintf( title, "%d.%d-1", z, id );

		p->SetAttribute( "title", title );
	}
}

/* --------------------------------------------------------------- */
/* UpdateXML ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateXML()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.inpath, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

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

		UpdateXMLLayer( layer, z );
	}

/* ---- */
/* Save */
/* ---- */

	xml.Save( "xmltmp.txt", true );
}

/* --------------------------------------------------------------- */
/* UpdateRick ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateRick()
{
	FILE	*in  = FileOpenOrDie( gArgs.inpath, "r" );

// name and open out file: path/name.ext -> path/name_v2.ext

	char	buf[2048];
	int		len = sprintf( buf, "%s", gArgs.inpath ) - 4;
	sprintf( buf + len, "_v2%s", gArgs.inpath + len );

	FILE	*out = FileOpenOrDie( buf, "w" );

// process line by line

	for( ;; ) {

		TAffine	T;
		int		z, col = -999, row = -999, cam = 0;

		if( fscanf( in, "%s%lf%lf%d\n",
			buf, &T.t[2], &T.t[5], &z ) != 4 ) {

			break;
		}

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		const char *c, *n = FileNamePtr( buf );

		if( c = strstr( n, "col" ) )
			sscanf( c, "col%d_row%d_cam%d", &col, &row, &cam );

		fprintf( out, "%d\t%d"
		"\t%f\t%f\t%f\t%f\t%f\t%f"
		"\t%d\t%d\t%d\t%s\n",
		z, gArgs.IDFromName( n ),
		T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
		col, row, cam, buf );
	}

	fclose( out );
	fclose( in );
}

/* --------------------------------------------------------------- */
/* UpdateIDB ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateIDB()
{
// loop over dirs

	for( int z = gArgs.zmin; z <= gArgs.zmax; ++z ) {

		// z exists?

		char	buf[4096];
		sprintf( buf, "%s/%d/TileToImage.txt", gArgs.inpath, z );

		if( !DskExists( buf ) )
			continue;

		// open input and skip header

		FILE		*in  = FileOpenOrDie( buf, "r" );
		CLineScan	LS;

		if( LS.Get( in ) <= 0 ) {
			fprintf( flog, "UpdateIDB: Empty file [%s].\n", buf );
			fclose( in );
			continue;
		}

		// name and open out file: same.tmp

		sprintf( buf, "%s/%d/TileToImage.tmp", gArgs.inpath, z );

		FILE	*out = FileOpenOrDie( buf, "w" );

		// write header

		fprintf( out,
		"ID\tT0\tT1\tX\tT3\tT4\tY\tCol\tRow\tCam\tPath\n" );

		// process line by line

		while( LS.Get( in ) > 0 ) {

			TAffine	T;
			int		id, col = -999, row = -999, cam = 0;

			sscanf( LS.line,
			"%d\t%lf\t%lf\t%lf"
			"\t%lf\t%lf\t%lf\t%[^\t\n]",
			&id,
			&T.t[0], &T.t[1], &T.t[2],
			&T.t[3], &T.t[4], &T.t[5],
			buf );

			const char *c, *n = FileNamePtr( buf );

			if( c = strstr( n, "col" ) )
				sscanf( c, "col%d_row%d_cam%d", &col, &row, &cam );

			fprintf( out, "%d"
			"\t%f\t%f\t%f\t%f\t%f\t%f"
			"\t%d\t%d\t%d\t%s\n",
			id,
			T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
			col, row, cam, buf );
		}

		fclose( out );
		fclose( in );

		// replace file

		sprintf( buf,
		"mv %s/%d/TileToImage.tmp %s/%d/TileToImage.txt",
		gArgs.inpath, z, gArgs.inpath, z );

		system( buf );
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

/* ------------- */
/* Write new xml */
/* ------------- */

	switch( gArgs.cmd ) {

		case 'x':
			UpdateXML();
		break;

		case 'r':
			UpdateRick();
		break;

		case 'd':
			UpdateIDB();
		break;
	}

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


