//
// MakeIDB reads either a Leginon (simple) file or TrakEM2 xml
// file and creates an 'image database' with this structure:
//
//	folder 'idbname'		// top folder
//		imagesize.txt		// has IMAGESIZE tag
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


#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"CTForm.h"

#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;	// file name
	double	r;		// inlayer radius from center
	int		z;		// Z layer
	int		id;		// inlayer id
	TForm	tr;		// local to global
	TForm	inv;	// global to local

public:
	Picture()	{id = -1;};
};

/* --------------------------------------------------------------- */
/* cArgs_idb ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class cArgs_idb {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images.
	CRegexID	re_id;

public:
	char	*infile,
			*outdir;
	int		zmin,
			zmax;
	bool	Simple,
			NoFolds;

public:
	cArgs_idb()
	{
		infile		=
		outdir		= "NoSuch";	// prevent overwriting real dir
		zmin		= 0;
		zmax		= 32768;
		Simple		= false;
		NoFolds		= false;
	};

	void SetCmdLine( int argc, char* argv[] );

	int DecodeID( const char *name );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char			gtopdir[2048];
static cArgs_idb	gArgs;
static FILE*		flog	= NULL;
static uint32		gW		= 0,	// universal pic dims
					gH		= 0;
static int			ismrc	= false;






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

	char	*pat;

	re_id.Set( "/N" );

	if( argc < 3 ) {
		printf( "Usage: makeidb <source-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( outdir, "-d", argv[i] ) )
			;
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( IsArg( "-simple", argv[i] ) )
			Simple = true;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
}

/* --------------------------------------------------------------- */
/* DecodeID ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int cArgs_idb::DecodeID( const char *name )
{
	const char	*s = strrchr( name, '/' );
	int			id;

	if( !s ) {
		fprintf( flog, "No '/' in [%s].\n", name );
		exit( 42 );
	}

	if( !re_id.Decode( id, ++s ) ) {
		printf( "No tile-id found in [%s].\n", s );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* ParseSimple --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseSimple( vector<Picture> &vp )
{
	FILE	*fp = FileOpenOrDie( gArgs.infile, "r", flog );

/* ---------- */
/* Scan lines */
/* ---------- */

	for( ;; ) {

		Picture	p;
		char	name[2048];
		int		x, y, z;

		/* ---------- */
		/* Get a line */
		/* ---------- */

		if( fscanf( fp, "%s %d %d %d", name, &x, &y, &z ) != 4 )
			break;

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		/* ---------------------- */
		/* Read actual dimensions */
		/* ---------------------- */

		if( !gW ) {

			uint8 *ras = Raster8FromAny( name, gW, gH, flog );

			if( !ras || !gW ) {
				fprintf( flog, "Error loading [%s].\n", name );
				exit( 42 );
			}

			RasterFree( ras );
		}

		/* ----------------- */
		/* Set picture entry */
		/* ----------------- */

		p.fname	= name;
		p.z		= z;
		p.id	= gArgs.DecodeID( name );
		p.tr.SetXY( x, y );

		vp.push_back( p );
	}

/* ----- */
/* Close */
/* ----- */

	fclose( fp );
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<Picture> &vp )
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	layer = hDoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

	//fprintf( flog, "Child element value %s.\n", layer->Value() );

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

		//fprintf( flog, "Got a <t2_layer>.\n" );

		const char	*sz = layer->Attribute( "z" );
		int			z	= int(atof(sz) + 0.5);

		//fprintf( flog, "z = %s.\n", sz );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		/* ------------------------------ */
		/* For each patch (tile) in layer */
		/* ------------------------------ */

		TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

		for( ; ptch; ptch = ptch->NextSiblingElement() ) {

			//fprintf( flog, "Got a <t2_patch>.\n" );

			Picture		p;
			const char	*name = ptch->Attribute( "file_path" );

			/* ---- */
			/* Dims */
			/* ---- */

			if( !gW ) {

				uint8 *ras = Raster8FromAny( name, gW, gH, flog );

				if( !ras || !gW ) {
					fprintf( flog, "Error loading [%s].\n", name );
					exit( 42 );
				}

				RasterFree( ras );
			}

			/* ----------------- */
			/* Set picture entry */
			/* ----------------- */

			p.fname	= name;
			p.z		= z;
			p.id	= gArgs.DecodeID( name );

			p.tr.ScanTrackEM2( ptch->Attribute( "transform" ) );

			vp.push_back( p );
		}
	}
}

/* --------------------------------------------------------------- */
/* Sort_z_inc ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_z_inc( const Picture &A, const Picture &B )
{
	return A.z < B.z || (A.z == B.z && A.id < B.id);
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
	char	name[2048];

// gtopdir gets the full path to the top directory
	getcwd( gtopdir, sizeof(gtopdir) );
	sprintf( name, "/%s", gArgs.outdir );
	strcat( gtopdir, name );

// create the top dir
	DskCreateDir( gArgs.outdir, flog );
}

/* --------------------------------------------------------------- */
/* WriteImageSizeFile -------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteImageSizeFile()
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/imagesize.txt", gArgs.outdir );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "IMAGESIZE %d %d\n", gW, gH );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteSubfmFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubfmFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/subfm", gArgs.outdir );

	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "setenv MRC_TRIM 12\n\n" );

	fprintf( f, "foreach i (`seq $1 $2`)\n" );
	fprintf( f, "cd $i\n" );

	fprintf( f,
	"qsub -N lou-f-$i"
	" -cwd -V -b y -pe batch 8"
	" make -f make.fm -j 8"
	" EXTRA='\"\"'\n" );

	fprintf( f, "cd ..\n" );
	fprintf( f, "end\n\n" );

	fclose( f );

	sprintf( buf, "chmod ug=rwx,o=rx %s/subfm", gArgs.outdir );
	system( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/report_fm", gArgs.outdir );

	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/csh\n\n" );

	fprintf( f, "ls -l */lou-f*.e* > FmErrs.txt\n\n" );

	fclose( f );

	sprintf( buf, "chmod ug=rwx,o=rx %s/report_fm", gArgs.outdir );
	system( buf );
}

/* --------------------------------------------------------------- */
/* GetLayerLimits ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Given starting index i0, which selects a layer z, set iN to be
// one beyond the highest index of a picture having same z. This
// makes loop limits [i0,iN) exclusive.
//
// If i0 or iN are out of bounds, both are set to -1.
//
static void GetLayerLimits(
	const vector<Picture>	&vp,
	int						&i0,
	int						&iN )
{
	int	np = vp.size();

	if( i0 < 0 || i0 >= np ) {

		i0 = -1;
		iN = -1;
	}
	else {

		int	Z = vp[i0].z;

		for( iN = i0 + 1; iN < np && vp[iN].z == Z; ++iN )
			;
	}
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
/* Make_TileToImage ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Write TileToImage.txt for this layer.
// Each row gives {tile, global-Tr, image path}.
//
static void Make_TileToImage(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
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

		// write the entries

		for( int i = is0; i < isN; ++i ) {

			const Picture&	P = vp[i];
			const double*	T = P.tr.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
				P.id, T[0], T[1], T[2], T[3], T[4], T[5],
				P.fname.c_str() );
		}
	}
	else {

		// create nmrc dir

		char	nmrcpath[2048];

		sprintf( nmrcpath, "%s/%d/nmrc", gtopdir, vp[is0].z );
		DskCreateDir( nmrcpath, flog );

		// write the entries

		for( int i = is0; i < isN; ++i ) {

			const Picture&	P = vp[i];
			const double*	T = P.tr.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s/nmrc_%d_%d.png\n",
				P.id, T[0], T[1], T[2], T[3], T[4], T[5],
				nmrcpath, P.z, P.id );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_TileToFM ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write file=TileToFM (or file=TileToFMD) for this layer.
// Each row gives {tile, image path}.
// Create corresponding folder=fm (or folder=fmd).
//
static void Make_TileToFM(
	const char				*lyrdir,
	const char				*file,
	const char				*folder,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
// Open file

	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/%s.txt", lyrdir, file );

	f = FileOpenOrDie( name, "w", flog );

// Header

	fprintf( f, "Tile\tPath\n" );

// Create fm dir

	char	fmpath[2048];

	sprintf( fmpath, "%s/%d/%s", gtopdir, vp[is0].z, folder );
	DskCreateDir( fmpath, flog );

// Write sorted entries

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];

		fprintf( f,
			"%d\t%s/%s_%d_%d.png\n",
			P.id, fmpath, folder, P.z, P.id );
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
static void Make_fmsame(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/fm.same", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Create fmpath

	char	fmpath[2048];

	sprintf( fmpath, "%s/%d/fm", gtopdir, vp[is0].z );

// FOLDMAP2 entries

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];

		fprintf( f, "FOLDMAP2 %d %d 1\n", P.z, P.id );
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
static void Make_MakeFM(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/make.fm", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Master target depends on all others

	fprintf( f, "all: " );

	for( int i = is0; i < isN; ++i )
		fprintf( f, "fm/fm_%d_%d.png ", vp[is0].z, vp[i].id );

	fprintf( f, "\n\n" );

// Subtargets and rules

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];
		char dep[2048];

		ConvertSpaces( dep, P.fname.c_str() );

		fprintf( f, "fm/fm_%d_%d.png: %s\n", P.z, P.id, dep );

		if( gArgs.NoFolds ) {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=nmrc/nmrc_%d_%d.png'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str(),
				P.z, P.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str() );
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
				P.z, P.id, P.fname.c_str(),
				P.z, P.id,
				P.z, P.id,
				P.z, P.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-fm=fm/fm_%d_%d.png'"
				" '-fmd=fmd/fmd_%d_%d.png'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str(),
				P.z, P.id,
				P.z, P.id );
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
static void ForEachLayer( const vector<Picture> &vp )
{
	int		is0, isN;

	GetLayerLimits( vp, is0 = 0, isN );

	while( isN != -1 ) {

		char	lyrdir[2048];

		CreateLayerDir( lyrdir, vp[is0].z );

		Make_TileToImage( lyrdir, vp, is0, isN );

		if( gArgs.NoFolds ) {

			if( ismrc )
				Make_MakeFM( lyrdir, vp, is0, isN );
			else
				Make_fmsame( lyrdir, vp, is0, isN );
		}
		else {
			Make_TileToFM( lyrdir, "TileToFM",  "fm",  vp, is0, isN );
			Make_TileToFM( lyrdir, "TileToFMD", "fmd", vp, is0, isN );
			Make_MakeFM( lyrdir, vp, is0, isN );
		}

		GetLayerLimits( vp, is0 = isN, isN );
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Picture>	vp;

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ---------------- */
/* Read source file */
/* ---------------- */

	if( gArgs.Simple )
		ParseSimple( vp );
	else
		ParseTrakEM2( vp );

	fprintf( flog, "Got %d images.\n", vp.size() );

	if( !vp.size() )
		goto exit;

	ismrc = strstr( vp[0].fname.c_str(), ".mrc" ) != NULL;

/* ------------------------- */
/* Sort tiles by layer, tile */
/* ------------------------- */

	sort( vp.begin(), vp.end(), Sort_z_inc );

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteImageSizeFile();

	if( !gArgs.NoFolds || ismrc ) {
		WriteSubfmFile();
		WriteReportFile();
	}

	ForEachLayer( vp );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


