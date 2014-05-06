//
// makeidb reads either a 'Rick' txt file or TrakEM2 xml
// file and creates an 'image database' with this structure:
//
//	folder 'idbname'		// top folder
//		imageparams.txt		// IDBPATH, IMAGESIZE tags
//		folder '0'			// folder per layer, here, '0'
//			TileToImage.txt	// TForm, image path from id
//			folder 'nmrc'	// mrc_to_png folder if needed
//	<if -nf option set...>
//			fm.same			// FOLDMAP2 entries {z,id,nrgn=1}
//	<if -nf option not set...>
//			make.fm			// make file for 'tiny'
//			TileToFM.txt	// fm path from id
//			TileToFMD.txt	// fmd path from id
//			folder 'fm'		// the foldmasks
//			folder 'fmd'	// the drawing foldmasks
//
// Notes:
//
// - In TileToImage, we adopt the paths from the input file,
// unless the input paths are mrc files. In that case, we
// expect that 'tiny' will be run to convert the mrc files
// into nmrc png files and those will reside in the database
// hierarchy.
//
// - TileToFM and TileToFMD are written by makeIDB but the
// foldmasks themselves will not exist until 'tiny' runs.
//
// - All entries in TileToXXX files are in tile id order.
//
// If output directory (idbname) is unspecified, either by
// omitting '-idb' option entirely, or by using '-idb' with no
// name, then no idb is generated. Rather, we write TrakEM2
// generator files.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"CTileSet.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* cArgs_idb ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class cArgs_idb {

public:
	// xml_type values: these are ImagePlus codes:
	// AUTO			= -1
	// GRAY8		= 0
	// GRAY16		= 1
	// GRAY32		= 2
	// COLOR_256	= 3
	// COLOR_RGB	= 4
	//
	const char	*infile,
				*outdir,
				*crop,
				*lens,
				*clk;
	int			zmin,
				zmax,
				xml_type,
				xml_min,
				xml_max;
	bool		NoFolds;

public:
	cArgs_idb()
	{
		infile		=
		outdir		= "NoSuch";	// prevent overwriting real dir
		crop		= NULL;
		lens		= NULL;
		clk			= NULL;
		zmin		= 0;
		zmax		= 32768;
		xml_type	= 0;
		xml_min		= 0;
		xml_max		= 0;
		NoFolds		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char				gtopdir[2048];
static cArgs_idb		gArgs;
static CTileSet			TS;
static vector<string>	vname0;
static FILE*			flog	= NULL;
static int				ismrc	= false;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void cArgs_idb::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "makeidb.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make img database: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: makeidb <source-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( outdir, "-idb=", argv[i] ) )
			;
		else if( GetArgStr( crop, "-crop=", argv[i] ) )
			;
		else if( GetArgStr( lens, "-lens=", argv[i] ) )
			;
		else if( GetArgStr( clk, "-k=", argv[i] ) )
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
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Make_nmrc_paths ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void Make_nmrc_paths()
{
// topdir gets the full path to the top directory

	char	topdir[2048];

	DskAbsPath( topdir, sizeof(topdir), gArgs.outdir, flog );

// save mrc names in global vector and replace them in TS

	int	nt = TS.vtil.size();

	for( int i = 0; i < nt; ++i ) {

		CUTile&	U = TS.vtil[i];
		char	path[2048];

		sprintf( path, "%s/%d/nmrc/nmrc_%d_%d.png",
			topdir, U.z, U.z, U.id );

		vname0.push_back( U.name );
		U.name = path;
	}
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
}

/* --------------------------------------------------------------- */
/* WriteImageparamsFile ------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteImageparamsFile()
{
	char	name[2048];
	FILE	*f;
	int		w, h;

	TS.GetTileDims( w, h );

	sprintf( name, "%s/imageparams.txt", gArgs.outdir );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "IDBPATH %s\n", gtopdir );
	fprintf( f, "IMAGESIZE %d %d\n", w, h );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* CopyCropFile -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CopyCropFile()
{
	if( gArgs.crop ) {
		char	buf[2048];
		sprintf( buf, "cp %s %s/crop.txt", gArgs.crop, gArgs.outdir );
		system( buf );
	}
}

/* --------------------------------------------------------------- */
/* CopyLensFile -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CopyLensFile()
{
	if( gArgs.lens ) {
		char	buf[2048];
		sprintf( buf, "cp %s %s/lens.txt", gArgs.lens, gArgs.outdir );
		system( buf );
	}
}

/* --------------------------------------------------------------- */
/* WriteSubfmFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubfmFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/subfm.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, generate foldmasks for images.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./subfm.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "nthr=4\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\t\tqsub -N makefm-$lyr -cwd -V -b y -pe batch $nthr make -f make.fm -j $nthr EXTRA='\"\"'\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

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

	sprintf( buf, "%s/report_fm.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# For layer range, tabulate sizes of all cluster stderr logs\n" );
	fprintf( f, "# from foldmask generation jobs.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./report_fm.sht <zmin> [zmax]\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l */makefm*.e* > FmErrs.txt\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* CreateLayerDir ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Each layer gets a directory named by its z-index.
// All content pertains to this layer.
//
static void CreateLayerDir( char *lyrdir, int L )
{
	fprintf( flog, "\n\nCreateLayerDir: layer %d\n", L );

	sprintf( lyrdir, "%s/%d", gArgs.outdir, L );
	DskCreateDir( lyrdir, flog );
}

/* --------------------------------------------------------------- */
/* Make_TileToFM ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write file=TileToFM (or file=TileToFMD) for this layer.
// Each row gives {tile, image path}.
// Create corresponding folder=fm (or folder=fmd).
//
static void Make_TileToFM(
	const char	*lyrdir,
	const char	*file,
	const char	*folder,
	int			is0,
	int			isN )
{
// Open file

	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/%s.txt", lyrdir, file );

	f = FileOpenOrDie( name, "w", flog );

// Header

	fprintf( f, "ID\tPath\n" );

// Create fm dir

	char	fmpath[2048];

	sprintf( fmpath, "%s/%d/%s", gtopdir, TS.vtil[is0].z, folder );
	DskCreateDir( fmpath, flog );

// Write sorted entries

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		fprintf( f,
			"%d\t%s/%s_%d_%d.png\n",
			U.id, fmpath, folder, U.z, U.id );
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

	sprintf( name, "%s/fm.same", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Create fmpath

	char	fmpath[2048];

	sprintf( fmpath, "%s/%d/fm", gtopdir, TS.vtil[is0].z );

// FOLDMAP2 entries

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		fprintf( f, "FOLDMAP2 %d.%d 1\n", U.z, U.id );
	}

	fclose( f );
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
/* Make_MakeFM --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write script to tell program tiny to calculate foldmasks for
// each tile (and create nmrc image if needed).
//
static void Make_MakeFM( const char *lyrdir, int is0, int isN )
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/make.fm", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Master target depends on all others

	fprintf( f, "all: " );

	for( int i = is0; i < isN; ++i )
		fprintf( f, "fm/fm_%d_%d.png ", TS.vtil[is0].z, TS.vtil[i].id );

	fprintf( f, "\n\n" );

// Subtargets and rules

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

#if 0
// target with dependency
		char dep[2048];
		ConvertSpaces( dep, vname0[i].c_str() );
		fprintf( f, "fm/fm_%d_%d.png: %s\n", U.z, U.id, dep );
#else
// target only
		fprintf( f, "fm/fm_%d_%d.png:\n", U.z, U.id );
#endif

		if( gArgs.NoFolds ) {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=nmrc/nmrc_%d_%d.png'"
				" -nf ${EXTRA}\n",
				U.z, U.id, vname0[i].c_str(),
				U.z, U.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" -nf ${EXTRA}\n",
				U.z, U.id, vname0[i].c_str() );
			}
		}
		else {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=nmrc/nmrc_%d_%d.png'"
				" '-fm=fm/fm_%d_%d.png'"
				" '-fmd=fmd/fmd_%d_%d.png'"
				" ${EXTRA}\n",
				U.z, U.id, vname0[i].c_str(),
				U.z, U.id,
				U.z, U.id,
				U.z, U.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-fm=fm/fm_%d_%d.png'"
				" '-fmd=fmd/fmd_%d_%d.png'"
				" ${EXTRA}\n",
				U.z, U.id, vname0[i].c_str(),
				U.z, U.id,
				U.z, U.id );
			}
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer()
{
	int	is0, isN;

	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 ) {

		char	lyrdir[2048];

		CreateLayerDir( lyrdir, TS.vtil[is0].z );

		TS.WriteTileToImage( gtopdir, false, ismrc, is0, isN );

		if( gArgs.NoFolds ) {

			if( ismrc )
				Make_MakeFM( lyrdir, is0, isN );
			else
				Make_fmsame( lyrdir, is0, isN );
		}
		else {
			Make_TileToFM( lyrdir, "TileToFM",  "fm",  is0, isN );
			Make_TileToFM( lyrdir, "TileToFMD", "fmd", is0, isN );
			Make_MakeFM( lyrdir, is0, isN );
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

/* ---------------- */
/* Read source file */
/* ---------------- */

	int		isrickfile = !FileIsExt( gArgs.infile, ".xml" );

	if( isrickfile )
		TS.FillFromRickFile( gArgs.infile, gArgs.zmin, gArgs.zmax );
	else
		TS.FillFromTrakEM2( gArgs.infile, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", (int)TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	if( isrickfile )
		TS.SetTileDimsFromImageFile();

	TS.SortAll_z_id();

	ismrc = FileIsExt( TS.vtil[0].name.c_str(), ".mrc" );

	if( ismrc )
		Make_nmrc_paths();

/* ----------- */
/* Diagnostics */
/* ----------- */

	if( isrickfile || ismrc ) {

		TS.WriteTrakEM2_EZ( "PreClicks.xml",
			gArgs.xml_type, gArgs.xml_min, gArgs.xml_max );
	}

	if( gArgs.clk ) {

		TS.ApplyClix( tsClixAffine, gArgs.clk );

		TS.WriteTrakEM2_EZ( "PostClicks.xml",
			gArgs.xml_type, gArgs.xml_min, gArgs.xml_max );
	}

/* ----------------------- */
/* Just make generator xml */
/* ----------------------- */

	if( !gArgs.outdir[0] || !strcmp( gArgs.outdir, "NoSuch" ) )
		goto exit;

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteImageparamsFile();
	CopyCropFile();
	CopyLensFile();

	if( !gArgs.NoFolds || ismrc ) {
		WriteSubfmFile();
		WriteReportFile();
	}

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


