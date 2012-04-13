

#include	"PipeFiles.h"

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"CTileSet.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	minolap	0.025

// Special override for Davi: allow tiny overlap
//#define	minolap	0.0003

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Pair {

public:
	int	a, b;

public:
	Pair( int _a, int _b )	{a = _a; b = _b;};
};

/* --------------------------------------------------------------- */
/* CArgs_scr ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scr {

public:
	char	*infile,
			*outdir,
			*pat;
	int		zmin,
			zmax;
	bool	Connect,		// just connect two layers
			Simple,
			NoFolds,
			NoDirs;

public:
	CArgs_scr()
	{
		infile		=
		outdir		= "NoSuch";	// prevent overwriting real dir
		pat			= "/N";
		zmin		= 0;
		zmax		= 32768;
		Connect		= false;
		Simple		= false;
		NoFolds		= false;
		NoDirs		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char			gtopdir[2048];
static CArgs_scr	gArgs;
static CTileSet		TS;
static FILE*		flog	= NULL;
static uint32		gW		= 0,	// universal pic dims
					gH		= 0;
static int			gZMax	= 0;
static int			ismrc	= false;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "scr.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make scripts: %s ", atime );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: scr <source-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( IsArg( "-connect", argv[i] ) )
			Connect = true;
		else if( GetArgStr( pat, "-p", argv[i] ) )
			;
		else if( IsArg( "-simple", argv[i] ) )
			Simple = true;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else if( IsArg( "-nd", argv[i] ) )
			NoFolds = NoDirs = true;
		else if( GetArgStr( outdir, "-d", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
	char	name[2048];

// gtopdir gets the full path to the top directory
	DskAbsPath( gtopdir, sizeof(gtopdir), gArgs.outdir, flog );

// create the top dir
	DskCreateDir( gArgs.outdir, flog );

// create stack subdir
	sprintf( name, "%s/stack", gArgs.outdir );
	DskCreateDir( name, flog );
}

/* --------------------------------------------------------------- */
/* WriteImageparamsFile ------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteImageparamsFile()
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/imageparams.txt", gArgs.outdir );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "IMAGESIZE %d %d\n", gW, gH );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* CreateLayerDir ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Each layer gets a directory named by its z-index. All content
// pertains to this layer, or to this layer acting as a source
// onto itself or other layers.
//
// For example, make.down will contain ptest jobs aligning this
// layer onto that below (z-1).
//
static void CreateLayerDir( char *lyrdir, int L )
{
	fprintf( flog, "\n\nCreateLayerDir: layer %d\n", L );

	sprintf( lyrdir, "%s/%d", gArgs.outdir, L );
	DskCreateDir( lyrdir, flog );
}

/* --------------------------------------------------------------- */
/* CreateTileSubdirs --------------------------------------------- */
/* --------------------------------------------------------------- */

// Each tile gets a directory named by its picture id. All content
// pertains to this tile, or this tile acting as a source onto
// other tiles.
//
// For example, folder 8/10 contains the foldmask fm.png for tile
// 10 in layer 8. If this folder contains file 7.11.tf.txt it
// lists the transforms mapping tile 8/10 onto tile 7/11.
//
static void CreateTileSubdirs( const char *lyrdir, int is0, int isN )
{
	fprintf( flog, "--CreateTileSubdirs: layer %d\n",
		TS.vtil[is0].z );

	for( int i = is0; i < isN; ++i ) {

		char	subdir[2048];

		sprintf( subdir, "%s/%d", lyrdir, TS.vtil[i].id );
		DskCreateDir( subdir, flog );
	}
}

/* --------------------------------------------------------------- */
/* Make_TileToImage ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Write TileToImage.txt for this layer.
// Each row gives {tile, global-Tr, image path}.
//
static void Make_TileToImage( const char *lyrdir, int is0, int isN )
{
// Locally sort entries in this file by tile-id

	vector<int>	order;

	isN = TS.GetOrder_id( order, is0, isN );

// Open file

	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/TileToImage.txt", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Header

	fprintf( f, "Tile\tT0\tT1\tX\tT3\tT4\tY\tPath\n" );

// Write sorted entries
// Use given name, unless mrc images.

	if( !ismrc ) {

		for( int i = 0; i < isN; ++i ) {

			const CUTile&	U = TS.vtil[order[i]];
			const double*	T = U.T.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
				U.id, T[0], T[1], T[2], T[3], T[4], T[5],
				U.name.c_str() );
		}
	}
	else {

		char	basepath[2048];

		sprintf( basepath, "%s/%d/", gtopdir, TS.vtil[order[0]].z );

		for( int i = 0; i < isN; ++i ) {

			const CUTile&	U = TS.vtil[order[i]];
			const double*	T = U.T.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s%d/nmrc_%d_%d.png\n",
				U.id, T[0], T[1], T[2], T[3], T[4], T[5],
				basepath, U.id, U.z, U.id );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_ThmPairFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void OneTprFile( const char *lyrdir, int a, int b )
{
    char	name[2048];
    FILE	*f;

	sprintf( name, "%s/ThmPair_%d_@_%d.txt", lyrdir, a, b );
	f = FileOpenOrDie( name, "w", flog );
	WriteThmPairHdr( f );
	fclose( f );
}


// For a layer, make ThmPair files for this and adjacent layers.
//
static void Make_ThmPairFile(
	const char				*lyrdir,
	int						is0,
	int						id0,
	int						iu0 )
{
	fprintf( flog, "--Make_ThmPairFile: layer %d\n", TS.vtil[is0].z );

	OneTprFile( lyrdir, TS.vtil[is0].z, TS.vtil[is0].z );

	if( id0 != -1 )
		OneTprFile( lyrdir, TS.vtil[is0].z, TS.vtil[id0].z );

	//if( iu0 != -1 )
	//	OneTprFile( lyrdir, TS.vtil[is0].z, TS.vtil[iu0].z );
}

/* --------------------------------------------------------------- */
/* ABOlap -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return area(intersection) / (gW*gH).
//
// We construct the polygon enclosing the intersection in
// b-space. Its vertices comprise the set of rectangle edge
// crossings + any rectangle vertices that lie interior to
// the other.
//
// Once we've collected the vertices we order them by angle
// about their common centroid so they can be assembled into
// directed and ordered line segments for the area calculator.
//
static bool ABOlap( int a, int b )
{
	double	A = TS.ABOlap( a, b );

	fprintf( flog, "----ABOlap: Tile %3d - %3d; area frac %f.\n",
	TS.vtil[a].id, TS.vtil[b].id, A );

	return A > minolap;
}

/* --------------------------------------------------------------- */
/* ConvertSpaces ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Make files have targets, dependencies and rules. Make can not
// scan dependency strings that contain embedded spaces. However,
// it works if " " is substituted by "\ ".
//
// Note too, that make does not like dependency strings to be
// enclosed in any quotes, so that will not solve this issue.
//
static void ConvertSpaces( char *out, const char *in )
{
	while( *out++ = *in++ ) {

		if( in[-1] == ' ' ) {

			out[-1]	= '\\';
			*out++	= ' ';
		}
	}
}

/* --------------------------------------------------------------- */
/* WriteThumbMakeFile -------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell thumbs to process the pairs
// of images described by (P). Argument mkname is a string from
// {"same", "up", "down"}.
//
static void WriteThumbMakeFile(
	const char				*lyrdir,
	const char				*mkname,
	const vector<Pair>		&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/thumbs.%s", lyrdir, mkname );

	f = FileOpenOrDie( name, "w", flog );

// write 'all' targets line

	fprintf( f, "all:\n" );

// rule lines

	for( int i = 0; i < np; ++i ) {

		const CUTile&	A = TS.vtil[P[i].a];
		const CUTile&	B = TS.vtil[P[i].b];

		fprintf( f,
		"\tthumbs %d/%d@%d/%d ${EXTRA}\n",
		A.z, A.id, B.z, B.id );
	}

	fprintf( f, "\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_ThumbsSame ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting thumbs jobs for pairs
// of intersecting images within this layer.
//
// (a, b) = (source, target).
//
static void Make_ThumbsSame( const char *lyrdir, int is0, int isN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_ThumbsSame: layer %d\n", TS.vtil[is0].z );

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = a + 1; b < isN; ++b ) {

			if( ABOlap( a, b ) )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

	WriteThumbMakeFile( lyrdir, "same", P );
}

/* --------------------------------------------------------------- */
/* Make_ThumbsDown ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting thumbs jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[below layer]).
//
static void Make_ThumbsDown(
	const char				*lyrdir,
	int						is0,
	int						isN,
	int						id0,
	int						idN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_ThumbsDown: layer %d @ %d\n",
		TS.vtil[is0].z, (id0 != -1 ? TS.vtil[id0].z : -1) );

// write dummy file even if no targets

	if( id0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = id0; b < idN; ++b ) {

			if( ABOlap( a, b ) )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteThumbMakeFile( lyrdir, "down", P );
}

/* --------------------------------------------------------------- */
/* Make_MakeFM --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write script to tell program tiny to calculate a foldmask for
// each tile.
//
static void Make_MakeFM( const char *lyrdir, int is0, int isN )
{
	char	name[2048];
	FILE	*f;

	fprintf( flog, "--Make_MakeFM: layer %d\n", TS.vtil[is0].z );

	sprintf( name, "%s/make.fm", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// master target depends on all others

	fprintf( f, "all: " );

	for( int i = is0; i < isN; ++i )
		fprintf( f, "%d/fm.png ", TS.vtil[i].id );

	fprintf( f, "\n\n" );

// subtargets and rules

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];
		char dep[2048];

		ConvertSpaces( dep, U.name.c_str() );

		fprintf( f, "%d/fm.png: %s\n", U.id, dep );

		if( gArgs.NoFolds ) {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=%d/nmrc_%d_%d.png'"
				" ${EXTRA}\n",
				U.z, U.id, U.name.c_str(),
				U.id, U.z, U.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" ${EXTRA}\n",
				U.z, U.id, U.name.c_str() );
			}
		}
		else {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=%d/nmrc_%d_%d.png'"
				" '-fm=%d/fm.png'"
				" '-fmd=%d/fmd.png'"
				" ${EXTRA}\n",
				U.z, U.id, U.name.c_str(),
				U.id, U.z, U.id,
				U.id,
				U.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-fm=%d/fm.png'"
				" '-fmd=%d/fmd.png'"
				" ${EXTRA}\n",
				U.z, U.id, U.name.c_str(),
				U.id,
				U.id );
			}
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_fmsame --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write fm.same file with FOLDMAP2 entry for each tile.
//
// This is used only for '-nf' option because these entries
// all get connected-region count = 1.
//
static void Make_fmsame( const char *lyrdir, int is0, int isN )
{
	char	name[2048];
	FILE	*f;

	fprintf( flog, "--Make_fmsame: layer %d\n", TS.vtil[is0].z );

	sprintf( name, "%s/fm.same", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// FOLDMAP2 entries

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		fprintf( f, "FOLDMAP2 %d.%d 1\n", U.z, U.id );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P). Argument mkname is a string from
// {"same", "up", "down"}.
//
static void WriteMakeFile(
	const char				*lyrdir,
	const char				*mkname,
	const vector<Pair>		&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/make.%s", lyrdir, mkname );

	f = FileOpenOrDie( name, "w", flog );

// write 'all' targets line

	fprintf( f, "all: " );

	for( int i = 0; i < np; ++i ) {

		const CUTile&	A = TS.vtil[P[i].a];
		const CUTile&	B = TS.vtil[P[i].b];

		fprintf( f, "%d/%d.%d.map.tif ", A.id, B.z, B.id );
	}

	fprintf( f, "\n\n" );

// Write each 'target: dependencies' line
//		and each 'rule' line

	if( gArgs.NoFolds ) {

		for( int i = 0; i < np; ++i ) {

			const CUTile&	A = TS.vtil[P[i].a];
			const CUTile&	B = TS.vtil[P[i].b];

			fprintf( f,
			"%d/%d.%d.map.tif:\n",
			A.id, B.z, B.id );

			fprintf( f,
			"\tptest %d/%d@%d/%d -nf ${EXTRA}\n\n",
			A.z, A.id, B.z, B.id );
		}
	}
	else {

		for( int i = 0; i < np; ++i ) {

			const CUTile&	A = TS.vtil[P[i].a];
			const CUTile&	B = TS.vtil[P[i].b];
			char depA[2048], depB[2048];

			ConvertSpaces( depA, A.name.c_str() );
			ConvertSpaces( depB, B.name.c_str() );

			fprintf( f,
			"%d/%d.%d.map.tif:"
			" %s %s"
			" ../%d/%d/fm.png"
			" ../%d/%d/fm.png\n",
			A.id, B.z, B.id,
			depA, depB,
			A.z, A.id,
			B.z, B.id );

			fprintf( f,
			"\tptest %d/%d@%d/%d ${EXTRA}\n\n",
			A.z, A.id, B.z, B.id );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_MakeSame ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images within this layer.
//
// (a, b) = (source, target).
//
static void Make_MakeSame( const char *lyrdir, int is0, int isN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeSame: layer %d\n", TS.vtil[is0].z );

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = a + 1; b < isN; ++b ) {

			if( ABOlap( a, b ) )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

	WriteMakeFile( lyrdir, "same", P );
}

/* --------------------------------------------------------------- */
/* Make_MakeDown ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[below layer]).
//
static void Make_MakeDown(
	const char				*lyrdir,
	int						is0,
	int						isN,
	int						id0,
	int						idN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeDown: layer %d @ %d\n",
		TS.vtil[is0].z, (id0 != -1 ? TS.vtil[id0].z : -1) );

// write dummy file even if no targets

	if( id0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = id0; b < idN; ++b ) {

			if( ABOlap( a, b ) )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteMakeFile( lyrdir, "down", P );
}

/* --------------------------------------------------------------- */
/* Make_MakeUp --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[above layer]).
//
static void Make_MakeUp(
	const char				*lyrdir,
	int						is0,
	int						isN,
	int						iu0,
	int						iuN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeUp: layer %d @ %d\n",
		TS.vtil[is0].z, (iu0 != -1 ? TS.vtil[iu0].z : -1) );

// write dummy file even if no targets

	if( iu0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = iu0; b < iuN; ++b ) {

			if( ABOlap( a, b ) )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteMakeFile( lyrdir, "up", P );
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer()
{
	int		id0, idN, is0, isN, iu0, iuN;

	id0 = -1;
	idN = -1;
	TS.GetLayerLimits( is0 = 0, isN );
	TS.GetLayerLimits( iu0 = isN, iuN );

	while( isN != -1 ) {

		char	lyrdir[2048];

		CreateLayerDir( lyrdir, TS.vtil[is0].z );

		if( !gArgs.NoDirs )
			CreateTileSubdirs( lyrdir, is0, isN );

		Make_TileToImage( lyrdir, is0, isN );

		Make_ThmPairFile( lyrdir, is0, id0, iu0 );

		//Make_ThumbsSame( lyrdir, is0, isN );
		//Make_ThumbsDown( lyrdir, is0, isN, id0, idN );

		if( gArgs.NoFolds ) {

			if( ismrc )
				Make_MakeFM( lyrdir, is0, isN );
			else
				Make_fmsame( lyrdir, is0, isN );
		}
		else
			Make_MakeFM( lyrdir, is0, isN );

		Make_MakeSame( lyrdir, is0, isN );
		Make_MakeDown( lyrdir, is0, isN, id0, idN );
		//Make_MakeUp( lyrdir, is0, isN, iu0, iuN );

		id0 = is0;
		idN = isN;
		is0 = iu0;
		isN = iuN;
		TS.GetLayerLimits(  iu0 = iuN, iuN );
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

/* ---------------- */
/* Read source file */
/* ---------------- */

	if( gArgs.Simple )
		TS.FillFromRickFile( gArgs.infile, gArgs.zmin, gArgs.zmax );
	else
		TS.FillFromTrakEM2( gArgs.infile, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SetTileDimsFromImageFile();

	ismrc = strstr( TS.vtil[0].name.c_str(), ".mrc" ) != NULL;

/* ------------------------------------------------- */
/* Within each layer, sort tiles by dist from center */
/* ------------------------------------------------- */

	TS.SortAll_z_r();

/* ------------------- */
/* Handle connect mode */
/* ------------------- */

// Connect mode expects exactly two adjacent layers
// and creates only make.up and make.down for them.

	if( gArgs.Connect ) {

		int	nt = TS.vtil.size();

		if( nt < 2 || TS.vtil[0].z != TS.vtil[nt-1].z - 1 ) {

			fprintf( flog,
			"Bogons! Expected two consecutive layers for"
			" connect mode.\n" );
			exit( 42 );
		}

		char	lyrdir[2048];
		int		is0, isN, iu0, iuN;

		// makeUp for layer 0 onto 1

		TS.GetLayerLimits( is0 = 0, isN );
		iu0 = isN;
		iuN = nt;

		//sprintf( lyrdir, "%s/0", gArgs.outdir );
		//Make_MakeUp( lyrdir, is0, isN, iu0, iuN );

		// makeDown for layer 1 onto 0

		sprintf( lyrdir, "%s/1", gArgs.outdir );
		Make_MakeDown( lyrdir, iu0, iuN, is0, isN );

		goto exit;
	}

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteImageparamsFile();

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


