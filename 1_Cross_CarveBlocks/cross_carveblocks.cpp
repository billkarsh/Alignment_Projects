//
// Carve cross layer work into blocks; write bsub.sht.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"Geometry.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Block {
public:
	vector<int>	vID;
};


class BlockSet {
public:
	vector<Block>	K;
	int				w, h,
					kx, ky,
					dx, dy,
					nb;
private:
	void OrientLayer( int is0, int isN );
	void SetDims();
	void PartitionTiles( int is0, int isN );
	void Consolidate();
	void ReportBlocks( int z );
public:
	void CarveIntoBlocks( int is0, int isN );
	void WriteParams( int za, int zb );
};

/* --------------------------------------------------------------- */
/* CArgs_alnmon -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_alnmon {

public:
	const char	*srcscaf,
				*script;
	int			zmin,
				zmax;
public:
	CArgs_alnmon()
	: srcscaf(NULL), script(NULL), zmin(0), zmax(32768) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_alnmon	gArgs;
static ScriptParams	scr;
static CTileSet		TS;
static FILE*		flog	= NULL;
static int			gW		= 0,	// universal pic dims
					gH		= 0;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_alnmon::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "cross_carveblocks.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Carve blocks: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: cross_carveblocks srcscaf"
		" -script=scriptpath -z=i,j.\n" );
		exit( 42 );
	}

	vector<int>	vi;

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			srcscaf = argv[i];
		else if( GetArgStr( script, "-script=", argv[i] ) )
			;
		else if( GetArgList( vi, "-z=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zmin = vi[0];
				zmax = vi[1];
			}
			else {
				fprintf( flog,
				"Bad format in -z [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* WriteTestblockFile -------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTestblockFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "testblock.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Put this script into a Dx_y folder to try or debug block.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -evalalldz\t\t;force evaluation of all maxdz layers\n" );
	fprintf( f, "# -abdbg\t\t\t;make diagnostic images and exit (Z^Z-1)\n" );
	fprintf( f, "# -abdbg=k\t\t\t;make diagnostic images and exit (Z^k)\n" );
	fprintf( f, "# -abctr=0\t\t\t;debug at this a-to-b angle\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "qsub -N x -cwd -V -b y -pe batch %d cross_thisblock -script=%s\n",
	scr.blockslots, gArgs.script );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteBSubFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteBSubFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "bsub.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Fifth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_thisblock -script=scriptpath\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Do one block alignment job, data read from 'blockdat.txt'.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -evalalldz\t\t;force evaluation of all maxdz layers\n" );
	fprintf( f, "# -abdbg\t\t\t;make diagnostic images and exit (Z^Z-1)\n" );
	fprintf( f, "# -abdbg=k\t\t\t;make diagnostic images and exit (Z^k)\n" );
	fprintf( f, "# -abctr=0\t\t\t;debug at this a-to-b angle\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "cd ..\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcd $jb\n" );
	fprintf( f, "\t\t\tqsub -N x$jb-$lyr -cwd -V -b y -pe batch %d cross_thisblock -script=%s\n",
	scr.blockslots, gArgs.script );
	fprintf( f, "\t\t\tcd ..\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteCountdowndirsFile ---------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCountdowndirsFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "countdowndirs.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Count all 'Dx_y' dirs in layer range\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./countdowndirs.sht i j\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "export MRC_TRIM=12\n" );
	fprintf( f, "\n" );
	fprintf( f, "if (($# == 1))\n" );
	fprintf( f, "then\n" );
	fprintf( f, "\tlast=$1\n" );
	fprintf( f, "else\n" );
	fprintf( f, "\tlast=$2\n" );
	fprintf( f, "fi\n" );
	fprintf( f, "\n" );
	fprintf( f, "cd ..\n" );
	fprintf( f, "\n" );
	fprintf( f, "cnt=0\n" );
	fprintf( f, "\n" );
	fprintf( f, "for lyr in $(seq $1 $last)\n" );
	fprintf( f, "do\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcnt=$(($cnt+1))\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\techo z= $lyr  cum= $cnt\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteReportFiles ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteReportFiles()
{
	char	buf[2048];
	FILE	*f;

// errors

	sprintf( buf, "breport.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Sixth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Run this after bsub completes to compile tables\n" );
	fprintf( f, "# of block alignment errors, FAILs and make sizes.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./breport.sht\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l ../*/D*/xD*.e* > BlockErrs.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l ../*/D*/xD*.o* > BlockOuts.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l ../*/D*/make.down > BlockMakes.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "grep FAIL ../*/D*/cross_thisblock.log > BlockFAIL.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "grep -e \"Final coverage\" ../*/D*/cross_thisblock.log > BlockCoverage.txt\n" );
	fprintf( f, "\n" );
	fprintf( f, "grep -e \"PeakHunt: Best\" ../*/D*/cross_thisblock.log > BlockTForms.txt\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* OrientLayer --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rotate layer to have smallest footprint.
//
void BlockSet::OrientLayer( int is0, int isN )
{
// Collect all tile corners in (C)

	vector<Point>	C, cnr;
	Set4Corners( cnr, gW, gH );

	for( int i = is0; i < isN; ++i ) {

		vector<Point>	c( 4 );
		memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
		TS.vtil[i].T.Transform( c );

		for( int i = 0; i < 4; ++i )
			C.push_back( c[i] );
	}

// Rotate layer upright and translate to (0,0)

	TAffine	R;
	DBox	B;
	int		deg = TightestBBox( B, C );

	R.NUSetRot( deg*PI/180 );

	for( int i = is0; i < isN; ++i ) {

		TAffine&	T = TS.vtil[i].T;

		T = R * T;
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
	dx = scr.crossblocksize;
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
/* PartitionTiles ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Fill tiles into blocks according to tile center.
//
void BlockSet::PartitionTiles( int is0, int isN )
{
	K.clear();
	K.resize( nb );

	int	W2 = gW/2,
		H2 = gH/2;

	for( int i = is0; i < isN; ++i ) {

		Point	p( W2, H2 );
		int		ix, iy;

		TS.vtil[i].T.Transform( p );

		ix = int(p.x / dx);

		if( ix < 0 )
			ix = 0;
		else if( ix >= kx )
			ix = kx - 1;

		iy = int(p.y / dy);

		if( iy < 0 )
			iy = 0;
		else if( iy >= ky )
			iy = ky - 1;

		K[ix + kx*iy].vID.push_back( i );
	}
}

/* --------------------------------------------------------------- */
/* Consolidate --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Distribute small tile sets into near neighbor cells.
//
void BlockSet::Consolidate()
{
	if( nb <= 1 )
		return;

	vector<Point>	pc;

	int	lowcount	=
			int(scr.crossblocksize * scr.crossblocksize * 0.25),
		W2			= gW/2,
		H2			= gH/2;

	bool	changed;

	do {	// until no changes

		changed = false;

		// find lowest occupancy

		int	cmin = 0, imin;

		for( int i = 0; i < nb; ++i ) {

			int	ic = K[i].vID.size();

			if( !ic || ic >= lowcount )
				continue;

			if( !cmin || ic < cmin ) {
				cmin	= ic;
				imin	= i;
			}
		}

		if( !cmin )
			break;

		// calc cell centers if not yet done

		if( !pc.size() ) {

			pc.resize( nb );

			for( int i = 0; i < nb; ++i ) {

				int	iy = i / kx,
					ix = i - kx * iy;

				pc[i].x = ix * dx + dx/2;
				pc[i].y = iy * dy + dy/2;
			}
		}

		// regroup each of the tiles in imin to nearest 4-neib

		for( int k = 0; k < cmin; ++k ) {

			Point	P( W2, H2 );

			TS.vtil[K[imin].vID[k]].T.Transform( P );

			double	dnear = 0, d;
			int		iy = imin / kx,
					ix = imin - kx * iy,
					inear, ii;

			if( iy > 0 && K[ii = imin-kx].vID.size() ) {

				dnear = P.y - pc[ii].y;
				inear = ii;
			}

			if( iy < ky-1 && K[ii = imin+kx].vID.size() &&
				(!dnear || (d = pc[ii].y - P.y) < dnear) ) {

				dnear = d;
				inear = ii;
			}

			if( ix > 0 && K[ii = imin-1].vID.size() &&
				(!dnear || (d = P.x - pc[ii].x) < dnear) ) {

				dnear = d;
				inear = ii;
			}

			if( ix < kx-1 && K[ii = imin+1].vID.size() &&
				(!dnear || (d = pc[ii].x - P.x) < dnear) ) {

				dnear = d;
				inear = ii;
			}

			if( !dnear )
				break;

			// merge

			changed = true;

			K[inear].vID.push_back( K[imin].vID[k] );
		}

		if( changed )
			K[imin].vID.clear();

	} while( changed );
}

/* --------------------------------------------------------------- */
/* ReportBlocks -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print tile count array.
//
void BlockSet::ReportBlocks( int z )
{
	int	ntiles = 0;

	fprintf( flog, "\nZ %d, Array %dx%d, Tiles(i,j):\n", z, kx, ky );

// Report occupancy

	for( int i = 0; i < nb; ++i ) {

		int	iy = i / kx,
			ix = i - kx * iy,
			ij = K[i].vID.size();

		fprintf( flog, "%d%c", ij, (ix == kx - 1 ? '\n' : '\t') );
		ntiles += ij;
	}

	fprintf( flog, "Total = %d\n", ntiles );

// Report block centers

	fprintf( flog, "Block centers:\n" );

	for( int i = 0; i < nb; ++i ) {

		int	iy = i / kx,
			ix = i - kx * iy;

		if( K[i].vID.size() ) {

			fprintf( flog, "D%d_%d %d %d\n",
			ix, iy,
			ix * dx + dx/2,
			iy * dy + dy/2 );
		}
	}
}

/* --------------------------------------------------------------- */
/* CarveIntoBlocks ----------------------------------------------- */
/* --------------------------------------------------------------- */

void BlockSet::CarveIntoBlocks( int is0, int isN )
{
	OrientLayer( is0, isN );
	SetDims();
	PartitionTiles( is0, isN );
	Consolidate();
	ReportBlocks( TS.vtil[is0].z );
}

/* --------------------------------------------------------------- */
/* WriteParams --------------------------------------------------- */
/* --------------------------------------------------------------- */

void BlockSet::WriteParams( int za, int zb )
{
// In simplest terms, the lowest B we should match to is
// zmin = min( za - blockmaxdz, gArgs.zmin ). Since zb is
// usually just za-1, we could write that equivalently as
// zmin = min( zb + 1 - blockmaxdz, gArgs.zmin ). This
// second form is preferred in cases where layers are missing
// so that zb < za-1, but we still want to give cross_thisblock
// several tries to get a good match (yes, even if we exceed
// blockmaxdz to get those tries).

	int zmin = max( zb + 1 - scr.blockmaxdz, gArgs.zmin );

	for( int i = 0; i < nb; ++i ) {

		int	nID;

		if( nID = K[i].vID.size() ) {

			char	path[256];
			int		iy = i / kx,
					ix = i - kx * iy;

			// make Dx_y folder;
			// cross_thisblock will make ThmPair files as needed
			sprintf( path, "../%d", za );
			CreateJobsDir( path, ix, iy, za, -1, flog );

			// write params file
			sprintf( path, "../%d/D%d_%d/blockdat.txt", za, ix, iy );
			FILE	*f = FileOpenOrDie( path, "w", flog );

			fprintf( f, "scaf=%s\n", gArgs.srcscaf );

			fprintf( f, "ZaZmin=%d,%d\n", za, zmin );

			// list actual tile-IDs (TS.vtil[].id)
			fprintf( f, "nIDs=%d\n", nID );
			for( int j = 0; j < nID; ++j )
				fprintf( f, "%d\n", TS.vtil[K[i].vID[j]].id );

			fclose( f );
		}
	}
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ForEachLayer()
{
	int		is0, isN, zb;

// Lowest layer - just get zb

	TS.GetLayerLimits( is0 = 0, isN );

	if( isN == -1 )
		return;

	zb = TS.vtil[is0].z;

// Loop over layers above first

	TS.GetLayerLimits( is0 = isN, isN );

	while( isN != -1 ) {

		BlockSet	BS;
		int			za = TS.vtil[is0].z;

		BS.CarveIntoBlocks( is0, isN );
		BS.WriteParams( za, zb );

		zb = za;
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

	if( !ReadScriptParams( scr, gArgs.script, flog ) )
		exit( 42 );

/* ------------- */
/* Read src file */
/* ------------- */

	string	idb;

	IDBFromTemp( idb, "..", flog );

	if( idb.empty() )
		exit( 42 );

	TS.FillFromRgns( gArgs.srcscaf, idb, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", (int)TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SetTileDimsFromImageFile();
	TS.GetTileDims( gW, gH );

	TS.SortAll_z_id();

/* ------------- */
/* Driver script */
/* ------------- */

	WriteTestblockFile();
	WriteBSubFile();
	WriteCountdowndirsFile();
	WriteReportFiles();

/* ----- */
/* Carve */
/* ----- */

	ForEachLayer();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


