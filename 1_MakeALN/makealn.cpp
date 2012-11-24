//
// makealn reads an IDB 'image database' and creates an alignment
// workspace with this structure:
//
//	folder 'alnname'				// top folder
//		imageparams.txt				// IDBPATH, IMAGESIZE tags
//		folder '0'					// folder per layer, here, '0'
//			folder '0'				// output folder per tile, here '0'
//			S0_0					// same layer jobs
//				make.same			// make file for same layer
//				ThmPair_0_@_0.txt	// table of thumbnail results
//			D0_0					// down layer jobs
//				make.down			// make file for cross layers
//				ThmPair_0_@_j.txt	// table of thumbnail results
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

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
	// minareafrac used for same and cross layer.
	// Typical vals:
	// 0.025	historic typical
	// 0.020	Davi and new typical
	// 0.0003	Davi older Harvard data
	//
	// downradiuspix = 0, or, sloppy match radius for downs
	//
	double	minareafrac;
	string	idbpath;
	char	*outdir,
			*exenam;
	int		downradiuspix,
			zmin,
			zmax,
			xml_type,
			xml_min,
			xml_max;
	bool	NoFolds,
			NoDirs;

public:
	CArgs_scr()
	{
		minareafrac		= 0.020;
		outdir			= "NoSuch";	// prevent overwriting real dir
		exenam			= "ptest";
		downradiuspix	= 0;
		zmin			= 0;
		zmax			= 32768;
		xml_type		= -999;
		xml_min			= -999;
		xml_max			= -999;
		NoFolds			= false;
		NoDirs			= false;
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
		"Usage: makealn <idbpath> -d=temp -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			idbpath = argv[i];
		else if( GetArgStr( outdir, "-d=", argv[i] ) )
			;
		else if( GetArgStr( exenam, "-exe=", argv[i] ) )
			;
		else if( GetArg( &minareafrac, "-minareafrac=%lf", argv[i] ) )
			;
		else if( GetArg( &downradiuspix, "-downradiuspix=%d", argv[i] ) )
			;
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
/* WriteRunlsqFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void _WriteRunlsqFile( const char *path )
{
	FILE	*f = FileOpenOrDie( path, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	char	xmlprms[256] = "";
	int		L = 0;

	if( gArgs.xml_type != -999 )
		L += sprintf( xmlprms + L, "-xmltype=%d ", gArgs.xml_type );

	if( gArgs.xml_min != -999 )
		L += sprintf( xmlprms + L, "-xmlmin=%d ", gArgs.xml_min );

	if( gArgs.xml_max != -999 )
		L += sprintf( xmlprms + L, "-xmlmax=%d ", gArgs.xml_max );

	fprintf( f,
	"lsq pts.all -scale=.1 -square=.1 %s> lsq.txt\n\n",
	xmlprms );

	fclose( f );
	FileScriptPerms( path );
}


static void WriteRunlsqFile()
{
	char	buf[2048];
	sprintf( buf, "%s/stack/runlsq.sht", gArgs.outdir );
	_WriteRunlsqFile( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmosFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmosFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/mosaic/submos.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "export MRC_TRIM=12\n\n" );

	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n\n" );

	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );

	fprintf( f,
	"\t\tqsub -N mos-$lyr -cwd -V -b y -pe batch 8"
	" \"mos ../stack/simple 0,0,-1,-1 $lyr,$lyr -warp%s"
	" > mos_$lyr.txt\"\n",
	(gArgs.NoFolds ? " -nf" : "") );

	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubNFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubNFile( int njobs )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sub%d.sht", gArgs.outdir, njobs );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "export MRC_TRIM=12\n\n" );

	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n\n" );

	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr/S0_0\n" );
	fprintf( f, "\t\tqsub -N qS0_0-$lyr -cwd -V -b y -pe batch 8 make -f make.same -j %d EXTRA='\"\"'\n", njobs );
	fprintf( f, "\t\tif (($lyr > $1))\n" );
	fprintf( f, "\t\tthen\n" );
	fprintf( f, "\t\t\tcd ../D0_0\n" );
	fprintf( f, "\t\t\tqsub -N qD0_0-$lyr -cwd -V -b y -pe batch 8 make -f make.down -j %d EXTRA='\"\"'\n", njobs );
	fprintf( f, "\t\tfi\n" );
	fprintf( f, "\t\tcd ../..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/report.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "ls -l */S0_0/qS0_0-*.e* > SameErrs.txt\n" );
	fprintf( f, "ls -l */D0_0/qD0_0-*.e* > DownErrs.txt\n\n" );

	fprintf( f, "ls -l */S0_0/pts.same > SamePts.txt\n" );
	fprintf( f, "ls -l */D0_0/pts.down > DownPts.txt\n\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteCombineFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCombineFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/combine.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "# Gather point pair files into stack/pts.all\n" );
	fprintf( f, "#\n\n" );

	fprintf( f, "rm -f pts.all\n\n" );

	fprintf( f, "#get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
	fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < imageparams.txt)\n\n" );

	fprintf( f, "cp imageparams.txt pts.all\n\n" );

	fprintf( f, "for lyr in $(seq $1 $2)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tcat $idb/$lyr/fm.same >> pts.all\n" );
	fprintf( f, "done\n\n" );

	fprintf( f, "for lyr in $(seq $1 $2)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif (($lyr == $1))\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcat $lyr/S0_0/pts.same >> pts.all\n" );
	fprintf( f, "\telse\n" );
	fprintf( f, "\t\tcat $lyr/{S0_0/pts.same,D0_0/pts.down} >> pts.all\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

	fprintf( f, "mv pts.all stack\n\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteFinishFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFinishFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/finish.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "./combine.sht %d %d\n", gArgs.zmin, gArgs.zmax );
	fprintf( f, "cd stack\n" );
	fprintf( f, "./runlsq.sht\n\n" );

	fclose( f );
	FileScriptPerms( buf );
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
/* CreateJobSubdirs ---------------------------------------------- */
/* --------------------------------------------------------------- */

// S0_0 subdir is working dir for same layer alignments.
// D0_0 subdir is working dir for all down jobs.
//
static void CreateJobSubdirs( const char *lyrdir, int is0, int id0 )
{
	int	zs = TS.vtil[is0].z,
		zd = (id0 >= 0 ? TS.vtil[id0].z : -1);

	fprintf( flog, "--CreateJobSubdirs: layer %d\n", zs );

	CreateJobsDir( lyrdir, 0, 0, zs, zs, flog );
	CreateJobsDir( lyrdir, 0, 0, zs, zd, flog );
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

	fprintf( flog, "----ABOlap: Tile %3d - %3d; area frac %f\n",
	TS.vtil[a].id, TS.vtil[b].id, A );

	return A > gArgs.minareafrac;
}

/* --------------------------------------------------------------- */
/* WriteThumbMakeFile -------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell thumbs to process the pairs
// of images described by (P).
//
static void WriteThumbMakeFile(
	const char			*lyrdir,
	int					SD,
	int					ix,
	int					iy,
	const vector<Pair>	&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/%c%d_%d/make.%s",
	lyrdir, SD, ix, iy, (SD == 'S' ? "same" : "down") );

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

	WriteThumbMakeFile( lyrdir, 'S', 0, 0, P );
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
	WriteThumbMakeFile( lyrdir, 'D', 0, 0, P );
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P).
//
static void WriteMakeFile(
	const char			*lyrdir,
	int					SD,
	int					ix,
	int					iy,
	const vector<Pair>	&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/%c%d_%d/make.%s",
	lyrdir, SD, ix, iy, (SD == 'S' ? "same" : "down") );

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

	WriteMakeFile( lyrdir, 'S', 0, 0, P );
}

/* --------------------------------------------------------------- */
/* Make_MakeDown ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[below layer]).
//
static void Make_MakeDownTight(
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
	WriteMakeFile( lyrdir, 'D', 0, 0, P );
}


static void Make_MakeDownLoose(
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

			double	D;
			int		w, h;

			TS.GetTileDims( w, h );

			Point	pa( w/2, h/2 ), pb = pa;

			TS.vtil[a].T.Transform( pa );
			TS.vtil[b].T.Transform( pb );

			D = pb.Dist( pa );

			fprintf( flog,
			"----Loose: Tile %3d - %3d; dist %d\n",
			TS.vtil[a].id, TS.vtil[b].id, (int)D );

			if( D < gArgs.downradiuspix )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteMakeFile( lyrdir, 'D', 0, 0, P );
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer()
{
	int		id0, idN, is0, isN;

	id0 = -1;
	idN = -1;
	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 ) {

		char	lyrdir[2048];

		CreateLayerDir( lyrdir, TS.vtil[is0].z );

		CreateJobSubdirs( lyrdir, is0, id0 );

		if( !gArgs.NoDirs )
			CreateTileSubdirs( lyrdir, is0, isN );

		//Make_ThumbsSame( lyrdir, is0, isN );
		//Make_ThumbsDown( lyrdir, is0, isN, id0, idN );

		Make_MakeSame( lyrdir, is0, isN );

		if( gArgs.downradiuspix )
			Make_MakeDownLoose( lyrdir, is0, isN, id0, idN );
		else
			Make_MakeDownTight( lyrdir, is0, isN, id0, idN );

		id0 = is0;
		idN = isN;
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

	WriteSubNFile( 4 );
	WriteSubNFile( 8 );
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


