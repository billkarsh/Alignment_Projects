//
// makealn reads an IDB 'image database' and creates an alignment
// workspace with this structure:
//
//	folder 'alnname'			// top folder
//		imageparams.txt			// IDBPATH, IMAGESIZE tags
//		folder '0'				// folder per layer, here, '0'
//			ThmPair_0_@_j.txt	// table of thumbnail results
//			make.down			// make file for cross layers
//			make.same			// make file for same layer
//			folder '0'			// output folder per tile, here '0'
//


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
	string	idbpath;
	char	*outdir,
			*exenam;
	int		zmin,
			zmax;
	bool	NoFolds,
			NoDirs;

public:
	CArgs_scr()
	{
		outdir		= "NoSuch";	// prevent overwriting real dir
		exenam		= "ptest";
		zmin		= 0;
		zmax		= 32768;
		NoFolds		= false;
		NoDirs		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scr	gArgs;
static CTileSet		TS;
static FILE*		flog	= NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "makealn.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make alignment workspace: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: makealn <idbpath> -dtemp -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			idbpath = argv[i];
		else if( GetArgStr( outdir, "-d", argv[i] ) )
			;
		else if( GetArgStr( exenam, "-exe=", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else if( IsArg( "-nd", argv[i] ) )
			NoDirs = true;
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

// create the top dir
	DskCreateDir( gArgs.outdir, flog );

// copy imageparams here
	sprintf( name, "cp %s/imageparams.txt %s",
		gArgs.idbpath.c_str(), gArgs.outdir );
	system( name );

// create stack subdir
	sprintf( name, "%s/stack", gArgs.outdir );
	DskCreateDir( name, flog );

// create mosaic subdir
	sprintf( name, "%s/mosaic", gArgs.outdir );
	DskCreateDir( name, flog );
}

/* --------------------------------------------------------------- */
/* ScriptPerms --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScriptPerms( const char *path )
{
	char	buf[2048];

	sprintf( buf, "chmod ug=rwx,o=rx %s", path );
	system( buf );
}

/* --------------------------------------------------------------- */
/* WriteRunlsqFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteRunlsqFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/stack/runlsq", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "lsq pts.all -scale=.1 -square=.1 > lsq.txt\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmosFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmosFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/mosaic/submos", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "setenv MRC_TRIM 12\n\n" );

	fprintf( f, "foreach i (`seq $1 $2`)\n" );

	fprintf( f,
	"\tqsub -N mos-$i -cwd -V -b y -pe batch 8"
	" \"mos ../stack/simple 0,0,-1,-1 $i,$i -warp%s"
	" > mos_$i.txt\"\n",
	(gArgs.NoFolds ? " -nf" : "") );

	fprintf( f, "end\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSub8File ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSub8File()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sub8", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "setenv MRC_TRIM 12\n\n" );

	fprintf( f, "if ($#argv == 1) then\n" );
	fprintf( f, "\tset last = $1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tset last = $2\n" );
	fprintf( f, "endif\n\n" );

	fprintf( f, "foreach i (`seq $1 $last`)\n" );
	fprintf( f, "\techo $i\n" );
	fprintf( f, "\tif (-d $i) then\n" );
	fprintf( f, "\t\tcd $i\n" );
	fprintf( f, "\t\tqsub -N lou-s-$i -cwd -V -b y -pe batch 8 make -f make.same -j 8 EXTRA='\"\"'\n" );
	fprintf( f, "\t\tqsub -N lou-d-$i -cwd -V -b y -pe batch 8 make -f make.down -j 8 EXTRA='\"\"'\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tendif\n" );
	fprintf( f, "end\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSub4File ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSub4File()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sub4", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "setenv MRC_TRIM 12\n\n" );

	fprintf( f, "if ($#argv == 1) then\n" );
	fprintf( f, "\tset last = $1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tset last = $2\n" );
	fprintf( f, "endif\n\n" );

	fprintf( f, "foreach i (`seq $1 $last`)\n" );
	fprintf( f, "\techo $i\n" );
	fprintf( f, "\tif (-d $i) then\n" );
	fprintf( f, "\t\tcd $i\n" );
	fprintf( f, "\t\tqsub -N lou-s-$i -cwd -V -b y -pe batch 8 make -f make.same -j 4 EXTRA='\"\"'\n" );
	fprintf( f, "\t\tqsub -N lou-d-$i -cwd -V -b y -pe batch 8 make -f make.down -j 4 EXTRA='\"\"'\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tendif\n" );
	fprintf( f, "end\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/report", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "ls -l */lou-s*.e* > SamErrs.txt\n" );
	fprintf( f, "ls -l */lou-d*.e* > DwnErrs.txt\n\n" );

	fprintf( f, "ls -l */pts.same > SamPts.txt\n" );
	fprintf( f, "ls -l */pts.down > DwnPts.txt\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteCombineFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCombineFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/combine", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "rm -f pts.all\n\n" );

	fprintf( f, "#get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
	fprintf( f, "set idb = `sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' <imageparams.txt`\n\n" );

	fprintf( f, "cp imageparams.txt pts.all\n\n" );

	fprintf( f, "foreach i (`seq $1 $2`)\n" );
	fprintf( f, "\tcat $idb/$i/fm.same >> pts.all\n" );
	fprintf( f, "end\n\n" );

	fprintf( f, "foreach i (`seq $1 $2`)\n" );
	fprintf( f, "\techo $i\n" );
	fprintf( f, "\tif ( $i == $1 ) then\n" );
	fprintf( f, "\t\tcat $i/pts.{same} >> pts.all\n" );
	fprintf( f, "\telse\n" );
	fprintf( f, "\t\tcat $i/pts.{down,same} >> pts.all\n" );
	fprintf( f, "\tendif\n" );
	fprintf( f, "end\n\n" );

	fprintf( f, "mv pts.all stack\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteFinishFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFinishFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/finish", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "./combine %d %d\n", gArgs.zmin, gArgs.zmax );
	fprintf( f, "cd stack\n" );
	fprintf( f, "./runlsq\n\n" );

	fclose( f );
	ScriptPerms( buf );
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

	const char	*option_nf = (gArgs.NoFolds ? " -nf" : "");

	for( int i = 0; i < np; ++i ) {

		const CUTile&	A = TS.vtil[P[i].a];
		const CUTile&	B = TS.vtil[P[i].b];

		fprintf( f,
		"%d/%d.%d.map.tif:\n",
		A.id, B.z, B.id );

		fprintf( f,
		"\t%s %d/%d@%d/%d%s ${EXTRA}\n\n",
		gArgs.exenam, A.z, A.id, B.z, B.id, option_nf );
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

		Make_ThmPairFile( lyrdir, is0, id0, iu0 );

		//Make_ThumbsSame( lyrdir, is0, isN );
		//Make_ThumbsDown( lyrdir, is0, isN, id0, idN );

		Make_MakeSame( lyrdir, is0, isN );
		Make_MakeDown( lyrdir, is0, isN, id0, idN );
		//Make_MakeUp( lyrdir, is0, isN, iu0, iuN );

		id0 = is0;
		idN = isN;
		is0 = iu0;
		isN = iuN;
		TS.GetLayerLimits( iu0 = iuN, iuN );
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

/* ---------------- */
/* Read source data */
/* ---------------- */

	TS.FillFromIDB( gArgs.idbpath, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SetTileDimsFromIDB( gArgs.idbpath );

/* ------------------------------------------------- */
/* Within each layer, sort tiles by dist from center */
/* ------------------------------------------------- */

	TS.SortAll_z_r();

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteRunlsqFile();
	WriteSubmosFile();

	WriteSub8File();
	WriteSub4File();
	WriteReportFile();
	WriteCombineFile();
	WriteFinishFile();

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


