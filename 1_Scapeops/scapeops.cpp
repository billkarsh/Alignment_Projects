//
// scapeops
//
// Given file, working dir, id pattern:
//
// xmlfile -d%s -p%s	# e.g. -dtemp -p_N_
//
// Perform montage drawing and/or strip aligning as follows:
//
//	If drawing a montage...
//
//	-mb -zb=%d -mbscl=%d
//
//		[-mbsdev=%d]
//
// If aligning strips...
//
//	-ab -za=%d -zb=%d -abwide=%d -abscl=%d -abcorr=%lf
//
//		[-absdev=%d] [-abdbg] [-abctr=%lf]
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"Scape.h"
#include	"CThmScan.h"
#include	"Geometry.h"
#include	"Maths.h"
#include	"ImageIO.h"
#include	"Timer.h"
#include	"Debug.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	minolap	0.020

// Typical for most data
//#define	minolap	0.025

// Special override for cur Davi
//#define	minolap	0.020

// Special override for old Davi: allow tiny overlap
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


class Block {

public:
	vector<Pair>	P;
};


class BlockSet {

private:
	enum {
		klowcount = 12
	};

public:
	vector<Block>	K;
	int				w, h,
					kx, ky,
					dx, dy,
					nb;

private:
	void OrientLayer( int is0, int isN );
	void SetDims();
	void PartitionJobs( int is0, int isN );
	void Consolidate();
	void ReportBlocks( int z );

public:
	void CarveIntoBlocks( int is0, int isN );
	void MakeJobs( const char *lyrdir, int z );
};

/* --------------------------------------------------------------- */
/* Superscape ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSuperscape {

public:
	DBox	B;			// oriented bounding box
	double	x0, y0;		// scape corner in oriented system
	uint8	*ras;		// scape pixels
	uint32	ws, hs;		// scape dims
	int		is0, isN,	// layer index range
			Bxc, Byc,	// oriented layer center
			Bxw, Byh,	// oriented layer span
			deg;		// rotate this much to orient

public:
	CSuperscape()
		{ras = NULL;};

	virtual ~CSuperscape()
		{KillRas();};

	void KillRas()
		{
			if( ras ) {
				RasterFree( ras );
				ras = NULL;
			}
		};

	void DrawRas( const char *name )
		{
			if( ras )
				Raster8ToPng8( name, ras, ws, hs );
		};

	bool Load( const char *name, FILE* flog )
		{
			x0	= y0 = 0.0;
			ras	= Raster8FromAny( name, ws, hs, flog );
			return (ras != NULL);
		};

	void FindLayerIndices( int z );
	void OrientLayer();

	bool MakeWholeRaster();
	bool MakeRasV();
	bool MakeRasH();

	void WriteMeta( char clbl, int z, int scl );

	void MakePoints( vector<double> &v, vector<Point> &p );
};

/* --------------------------------------------------------------- */
/* CArgs_scp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scp {

public:
	double	inv_abscl,
			abcorr,
			abctr;
	char	*infile,
			*outdir,
			*pat,
			*exenam;
	int		za,
			zb,
			abwide,
			mbscl,
			abscl,
			mbsdev,
			absdev,
			blksize;
	bool	ismb,
			isab,
			abdbg,
			NoFolds,
			NoDirs;

public:
	CArgs_scp()
	{
		abcorr		= 0.20;
		abctr		= 0.0;
		infile		= NULL;
		outdir		= "NoSuch";	// prevent overwriting real dir
		pat			= "/N";
		exenam		= "ptest";
		za			= -1;
		zb			= -1;
		abwide		= 5;
		mbscl		= 200;
		abscl		= 200;
		mbsdev		= 0;	// 12 useful for Davi EM
		absdev		= 0;	// 12 useful for Davi EM
		blksize		= 60;
		ismb		= false;
		isab		= false;
		abdbg		= false;
		NoFolds		= false;
		NoDirs		= false;

		inv_abscl	= 1.0/abscl;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scp	gArgs;
static CTileSet		TS;
static FILE*		flog	= NULL;
static int			gW		= 0,	// universal pic dims
					gH		= 0,
					gW2		= 0,
					gH2		= 0;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scp::SetCmdLine( int argc, char* argv[] )
{
// label output by layer b

	for( int i = 1; i < argc; ++i ) {

		if( GetArg( &zb, "-zb=%d", argv[i] ) )
			break;
	}

	if( zb < 0 ) {
		printf( "scapeops: Missing -zb option!!\n" );
		exit( 42 );
	}

// start log

	char	buf[256];

	sprintf( buf, "scplogs/scp_%d.log", zb );
	flog = FileOpenOrDie( buf, "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Scapeops start: %s ", atime );

// parse command line args

	if( argc < 6 ) {
		printf(
		"Usage: See scapeops.cpp comments.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &za, "-za=%d", argv[i] ) )
			;
		else if( GetArg( &zb, "-zb=%d", argv[i] ) )
			;
		else if( GetArgStr( outdir, "-d", argv[i] ) )
			;
		else if( GetArgStr( pat, "-p", argv[i] ) )
			;
		else if( GetArgStr( exenam, "-exe=", argv[i] ) )
			;
		else if( GetArg( &abwide, "-abwide=%d", argv[i] ) )
			;
		else if( GetArg( &mbscl, "-mbscl=%d", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			inv_abscl = 1.0/abscl;
		else if( GetArg( &mbsdev, "-mbsdev=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &abctr, "-abctr=%lf", argv[i] ) )
			;
		else if( GetArg( &blksize, "-b=%d", argv[i] ) )
			;
		else if( IsArg( "-mb", argv[i] ) )
			ismb = true;
		else if( IsArg( "-ab", argv[i] ) )
			isab = true;
		else if( IsArg( "-abdbg", argv[i] ) )
			abdbg = true;
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

	if( !ismb && !isab ) {
		fprintf( flog, "No operations specified.\n" );
		exit( 0 );
	}
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
	char	name[2048];

// create the top dir
	DskCreateDir( gArgs.outdir, flog );

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

	sprintf( buf, "%s/stack/runlsq.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

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
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubsNFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubsNFile( int njobs )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/subs%d.sht", gArgs.outdir, njobs );
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
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcd $jb\n" );
	fprintf( f, "\t\t\tqsub -N q$jb-$lyr -cwd -V -b y -pe batch 8 make -f make.same -j %d EXTRA='\"\"'\n", njobs );
	fprintf( f, "\t\t\tcd ..\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

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

	sprintf( buf, "%s/report.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "ls -l */S*/qS*.e* > SameErrs.txt\n\n" );

	fprintf( f, "ls -l */S*/pts.same > SamePts.txt\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteMontage1File --------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMontage1File()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/montage1.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "cd $1\n\n" );

	fprintf( f, "rm -f pts.all\n\n" );

	fprintf( f, "#get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
	fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < ../imageparams.txt)\n\n" );

	fprintf( f, "cp ../imageparams.txt pts.all\n" );
	fprintf( f, "cat $idb/$1/fm.same >> pts.all\n\n" );

	fprintf( f, "for jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tcat $jb/pts.same >> pts.all\n" );
	fprintf( f, "done\n\n" );

	fprintf( f, "mv pts.all montage\n\n" );

	fprintf( f, "cd montage\n" );
	fprintf( f, "./runlsq.sht\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmonFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmonFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/submon.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

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
	"\t\tqsub -N mon-$lyr -cwd -V -b y -pe batch 8"
	" \"./montage1.sht $lyr\"\n" );

	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

	fclose( f );
	ScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportMonsFile ------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportMonsFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/sumymons.sht", gArgs.outdir );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n\n" );

	fprintf( f, "rm -rf MonSumy.txt\n\n" );

	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tlog=$lyr/montage/lsq.txt\n" );
	fprintf( f, "\tif [ -f \"$log\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\techo Z $lyr `grep -e \"FINAL*\" $log` >> MonSumy.txt\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

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

// Create layer dir
	sprintf( lyrdir, "%s/%d", gArgs.outdir, L );
	DskCreateDir( lyrdir, flog );

// Create montage subdir
	char	buf[2048];
	int		len;

	len = sprintf( buf, "%s/montage", lyrdir );
	DskCreateDir( buf, flog );

// Create montage script
	sprintf( buf + len, "/runlsq.sht" );
	FILE	*f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n\n" );

	fprintf( f, "lsq pts.all -scale=.1 -square=.1 > lsq.txt\n\n" );

	fclose( f );
	ScriptPerms( buf );
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
/* OrientLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rotate layer to have smallest footprint.
//
void BlockSet::OrientLayer( int is0, int isN )
{
// Collect all tile corners

	vector<Point>	C;

	for( int i = is0; i < isN; ++i ) {

		vector<Point>	p( 4 );

		p[0] = Point( 0.0 , 0.0 );
		p[1] = Point( gW-1, 0.0 );
		p[2] = Point( gW-1, gH-1 );
		p[3] = Point( 0.0 , gH-1 );

		TS.vtil[i].T.Transform( p );

		for( int i = 0; i < 4; ++i )
			C.push_back( p[i] );
	}

// Rotate layer upright and translate to (0,0)

	TForm	R;
	DBox	B;
	int		deg = TightestBBox( B, C );

	R.NUSetRot( deg*PI/180 );

	for( int i = is0; i < isN; ++i ) {

		TForm&	T = TS.vtil[i].T;

		MultiplyTrans( T, R, TForm( T ) );
		T.AddXY( -B.L, -B.B );
	}

	w = int(B.R - B.L) + 1;
	h = int(B.T - B.B) + 1;
}

/* --------------------------------------------------------------- */
/* SetDims ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Divide layer bbox into (kx,ky) cells of size (dx,dy).
//
void BlockSet::SetDims()
{
	dx = (int)sqrt( gArgs.blksize );
	dy = dx;
	dx *= gW;
	dy *= gH;
	kx = (int)ceil( (double)w / dx );
	ky = (int)ceil( (double)h / dy );
	dx = w / kx;
	dy = h / ky;
	nb = kx * ky;
}

/* --------------------------------------------------------------- */
/* PartitionJobs ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fill pair jobs into blocks according to a-tile center.
//
void BlockSet::PartitionJobs( int is0, int isN )
{
	K.clear();
	K.resize( nb );

	for( int a = is0; a < isN; ++a ) {

		Point	pa( gW2, gH2 );
		int		ix, iy;

		TS.vtil[a].T.Transform( pa );

		ix = int(pa.x / dx);

		if( ix < 0 )
			ix = 0;
		else if( ix >= kx )
			ix = kx - 1;

		iy = int(pa.y / dy);

		if( iy < 0 )
			iy = 0;
		else if( iy >= ky )
			iy = ky - 1;

		for( int b = a + 1; b < isN; ++b ) {

			if( TS.ABOlap( a, b ) > minolap )
				K[ix + kx*iy].P.push_back( Pair( a, b ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* Consolidate --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Append small job sets into more central neighbors.
//
void BlockSet::Consolidate()
{
	if( nb <= 1 )
		return;

	bool	changed;

	do {

		changed = false;

		for( int i = 0; i < nb; ++i ) {

			int	ic = K[i].P.size();

			if( !ic || ic >= klowcount )
				continue;

			int	iy = i / kx,
				ix = i - kx * iy,
				lowc = 0,
				lowi, c;

			// find lowest count neib

			if( iy > 0 && (c = K[i-kx].P.size()) ) {
				lowc = c;
				lowi = i-kx;
			}

			if( iy < ky-1 && (c = K[i+kx].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i+kx;
			}

			if( ix > 0 && (c = K[i-1].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i-1;
			}

			if( ix < kx-1 && (c = K[i+1].P.size()) &&
				(!lowc || c < lowc) ) {

				lowc = c;
				lowi = i+1;
			}

			// merge

			if( !lowc )
				continue;

			changed = true;

			for( int j = 0; j < ic; ++j )
				K[lowi].P.push_back( K[i].P[j] );

			K[i].P.clear();
		}

	} while( changed );
}

/* --------------------------------------------------------------- */
/* ReportBlocks -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print job count array.
//
void BlockSet::ReportBlocks( int z )
{
	int	njobs = 0;

	fprintf( flog, "\nZ %d, Array %dx%d, Jobs(i,j):\n", z, kx, ky );

	for( int i = 0; i < nb; ++i ) {

		int	iy = i / kx,
			ix = i - kx * iy,
			ij = K[i].P.size();

		fprintf( flog, "%d%c", ij, (ix == kx - 1 ? '\n' : '\t') );
		njobs += ij;
	}

	fprintf( flog, "Total = %d\n", njobs );
}

/* --------------------------------------------------------------- */
/* CarveIntoBlocks ----------------------------------------------- */
/* --------------------------------------------------------------- */

void BlockSet::CarveIntoBlocks( int is0, int isN )
{
	OrientLayer( is0, isN );
	SetDims();
	PartitionJobs( is0, isN );
	Consolidate();
	ReportBlocks( TS.vtil[is0].z );
}

/* --------------------------------------------------------------- */
/* MakeJobs ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void BlockSet::MakeJobs( const char *lyrdir, int z )
{
	for( int i = 0; i < nb; ++i ) {

		if( K[i].P.size() ) {

			int	iy = i / kx,
				ix = i - kx * iy;

			CreateJobsDir( lyrdir, ix, iy, z, z, flog );
			WriteMakeFile( lyrdir, 'S', ix, iy, K[i].P );
		}
	}
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer()
{
	int		is0, isN;

	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 ) {

		char		lyrdir[2048];
		BlockSet	BS;
		int			z = TS.vtil[is0].z;

		CreateLayerDir( lyrdir, z );

		if( !gArgs.NoDirs )
			CreateTileSubdirs( lyrdir, is0, isN );

		BS.CarveIntoBlocks( is0, isN );
		BS.MakeJobs( lyrdir, z );

		TS.GetLayerLimits( is0 = isN, isN );
	}
}

/* --------------------------------------------------------------- */
/* FindLayerIndices ---------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::FindLayerIndices( int z )
{
	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 && TS.vtil[is0].z != z )
		TS.GetLayerLimits( is0 = isN, isN );
}

/* --------------------------------------------------------------- */
/* OrientLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rotate layer to have smallest footprint.
//
void CSuperscape::OrientLayer()
{
// Collect all tile corners

	vector<Point>	C;

	for( int i = is0; i < isN; ++i ) {

		vector<Point>	p( 4 );

		p[0] = Point( 0.0 , 0.0 );
		p[1] = Point( gW-1, 0.0 );
		p[2] = Point( gW-1, gH-1 );
		p[3] = Point( 0.0 , gH-1 );

		TS.vtil[i].T.Transform( p );

		for( int i = 0; i < 4; ++i )
			C.push_back( p[i] );
	}

// Rotate layer upright and translate to (0,0)

	TForm	R;

	deg = TightestBBox( B, C );

	R.NUSetRot( deg*PI/180 );

	for( int i = is0; i < isN; ++i ) {

		TForm&	T = TS.vtil[i].T;

		MultiplyTrans( T, R, TForm( T ) );
	}

	Bxc = int((B.R + B.L)/2.0);
	Byc = int((B.B + B.T)/2.0);
	Bxw = int(B.R - B.L);
	Byh = int(B.T - B.B);
}

/* --------------------------------------------------------------- */
/* MakeWholeRaster ----------------------------------------------- */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeWholeRaster()
{
	vector<ScpTile>	S;

	for( int i = is0; i < isN; ++i ) {

		const CUTile& U = TS.vtil[i];

		S.push_back( ScpTile( U.name, U.T ) );
	}

	ras = Scape( ws, hs, x0, y0, S, gW, gH,
			1.0/gArgs.mbscl, 1, 0, gArgs.mbsdev, flog );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasV ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasV()
{
// Collect strip tiles

	vector<ScpTile>	S;
	int				w1, w2, h1, h2;

	w1 = (gArgs.abwide * gW)/2;
	w2 = Bxc + w1;
	w1 = Bxc - w1;

	h1 = int(Byh * 0.45);
	h2 = Byc + h1;
	h1 = Byc - h1;

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		if( U.T.t[2] >= w1 && U.T.t[2] <= w2 &&
			U.T.t[5] >= h1 && U.T.t[5] <= h2 ) {

			S.push_back( ScpTile( U.name, U.T ) );
		}
	}

	ras = Scape( ws, hs, x0, y0, S, gW, gH,
			gArgs.inv_abscl, 1, 0, gArgs.absdev, flog );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasH ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasH()
{
// Collect strip tiles

	vector<ScpTile>	S;
	int				w1, w2, h1, h2;

	w1 = int(Bxw * 0.45);
	w2 = Bxc + w1;
	w1 = Bxc - w1;

	h1 = (gArgs.abwide * gH)/2;
	h2 = Byc + h1;
	h1 = Byc - h1;

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		if( U.T.t[2] >= w1 && U.T.t[2] <= w2 &&
			U.T.t[5] >= h1 && U.T.t[5] <= h2 ) {

			S.push_back( ScpTile( U.name, U.T ) );
		}
	}

	ras = Scape( ws, hs, x0, y0, S, gW, gH,
			gArgs.inv_abscl, 1, 0, gArgs.absdev, flog );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* WriteMeta ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::WriteMeta( char clbl, int z, int scl )
{
	fprintf( flog,
	"*%c: z deg [l,r,b,t] scl [ws,hs] [x0,y0]\n", clbl );

	fprintf( flog,
	"%d %d [%g,%g,%g,%g] %d [%d,%d] [%g,%g]\n",
	z, deg, B.L, B.R, B.B, B.T, scl, ws, hs, x0, y0 );
}

/* --------------------------------------------------------------- */
/* MakePoints ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::MakePoints( vector<double> &v, vector<Point> &p )
{
	int		np = ws * hs;

	for( int i = 0; i < np; ++i ) {

		if( ras[i] ) {

			int	iy = i / ws,
				ix = i - ws * iy;

			v.push_back( ras[i] );
			p.push_back( Point( ix, iy ) );
		}
	}

	Normalize( v );

	KillRas();
}

/* --------------------------------------------------------------- */
/* MakeStripRasters ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeStripRasters( CSuperscape &A, CSuperscape &B )
{
	char	buf[256];

#if 1
	A.MakeRasH();
	sprintf( buf, "strips/A_%d.png", gArgs.za );
	A.DrawRas( buf );

	B.MakeRasV();
	sprintf( buf, "strips/B_%d.png", gArgs.zb );
	B.DrawRas( buf );
#else
// simple debugging - load existing image, but does
// NOT acquire {x0,y0,B,...} meta data!!

	sprintf( buf, "Astrip_%d.png", gArgs.za );
	A.Load( buf, flog );

	sprintf( buf, "Bstrip_%d.png", gArgs.zb );
	B.Load( buf, flog );
#endif
}

/* --------------------------------------------------------------- */
/* ScapeStuff ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// The strip alignment from S.DenovoBestAngle (followed by SetXY)
// produces a transform A->B = best.T for the scaled down scapes.
// Here's how to convert that to a transform between the original
// montages:
//
// Recapitulate total process...
//	TForm	Rbi, Ra, t;
//
// Scale back up
//	best.T.MulXY( gArgs.abscl );
//	A.x0 *= gArgs.abscl;
//	A.y0 *= gArgs.abscl;
//	B.x0 *= gArgs.abscl;
//	B.y0 *= gArgs.abscl;
//
// A-montage -> A-oriented
//	Ra.NUSetRot( A.deg*PI/180 );
//
// A-oriented -> A-scape
//	Ra.AddXY( -A.x0, -A.y0 );
//
// A-scape -> B-scape
//	MultiplyTrans( t, best.T, Ra );
//
// B-scape -> B->oriented
//	t.AddXY( B.x0, B.y0 );
//
// B->oriented -> B-montage
//	Rbi.NUSetRot( -B.deg*PI/180 );
//	MultiplyTrans( best.T, Rbi, t );
//
static void ScapeStuff()
{
	clock_t		t0 = StartTiming();
	CSuperscape	A, B;
	ThmRec		thm;
	CThmScan	S;
	CorRec		best;

	B.FindLayerIndices( gArgs.zb );
	B.OrientLayer();

	if( gArgs.ismb ) {

		char	buf[256];

		B.MakeWholeRaster();
		sprintf( buf, "montages/M_%d_0.png", gArgs.zb );
		B.DrawRas( buf );
		B.KillRas();
		B.WriteMeta( 'M', gArgs.zb, gArgs.mbscl );
		StopTiming( flog, "MakeMontage", t0 );
		t0 = StartTiming();
	}

	if( !gArgs.isab )
		return;

	A.FindLayerIndices( gArgs.za );
	A.OrientLayer();

	MakeStripRasters( A, B );
	StopTiming( flog, "MakeStrips", t0 );
	t0 = StartTiming();

	A.MakePoints( thm.av, thm.ap );
	A.KillRas();
	A.WriteMeta( 'A', gArgs.za, gArgs.abscl );

	B.MakePoints( thm.bv, thm.bp );
	B.KillRas();
	B.WriteMeta( 'B', gArgs.zb, gArgs.abscl );

	thm.ftc.clear();
	thm.reqArea	= int(gW * gH * gArgs.inv_abscl * gArgs.inv_abscl);
	thm.olap1D	= int(gW * gArgs.inv_abscl * 0.20);
	thm.scl		= 1;

	int	Ox	= -int(A.ws/2),
		Oy	=  int(B.hs/2),
		Rx	=  int(-Ox * 0.80),
		Ry	=  int( Oy * 0.80);

	S.Initialize( flog, best );
	S.SetTuser( 1, 1, 1, 0, 0 );
	S.SetRThresh( gArgs.abcorr );
	S.SetNbMaxHt( 0.99 );
	S.SetSweepType( true, false );
	S.SetUseCorrR( false );
	S.SetDisc( Ox, Oy, Rx, Ry );

	if( gArgs.abdbg ) {

		dbgCor = true;
		S.RFromAngle( best, gArgs.abctr, thm );
	}
	else {

		S.DenovoBestAngle( best, 0, 175, 5, thm );
		best.T.SetXY( best.X, best.Y );

		fprintf( flog, "*T: [0,1,2,3,4,5] (strip-strip)\n" );
		fprintf( flog, "[%g,%g,%g,%g,%g,%g]\n",
		best.T.t[0], best.T.t[1], best.T.t[2],
		best.T.t[3], best.T.t[4], best.T.t[5] );
	}

	StopTiming( flog, "Corr", t0 );
	t0 = StartTiming();
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	clock_t		t0 = StartTiming();

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

	TS.SetLogFile( flog );
	TS.SetDecoderPat( gArgs.pat );

/* ---------------- */
/* Read source data */
/* ---------------- */

	if( gArgs.zb >= 0 && gArgs.za < 0 )
		gArgs.za = gArgs.zb;

	TS.FillFromTrakEM2( gArgs.infile, gArgs.zb, gArgs.za );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.GetTileDims( gW, gH );
	gW2 = gW/2;
	gH2 = gH/2;

	StopTiming( flog, "ReadFile", t0 );
	t0 = StartTiming();

/* ------------- */
/* Sort by layer */
/* ------------- */

	TS.SortAll_z();

/* ----- */
/* Stuff */
/* ----- */

	ScapeStuff();

/* --------------- */
/* Create dir tree */
/* --------------- */

#if 0
	CreateTopDir();

	WriteRunlsqFile();
	WriteSubmosFile();

	WriteSubsNFile( 4 );
	WriteSubsNFile( 8 );
	WriteReportFile();
	WriteMontage1File();
	WriteSubmonFile();
	WriteReportMonsFile();

	ForEachLayer();
#endif

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


