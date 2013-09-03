//
// Carve cross layer work into blocks; write subblocks.sht.
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
	double	blkmincorr,
			blknomcorr,
			xyconf;		// search radius = (1-conf)(blockwide)
	char	xml_hires[2048];
	int		zmin,
			zmax,
			blksize,
			abscl,
			ablgord,
			absdev,
			maxDZ;
	bool	NoFolds;

public:
	CArgs_alnmon()
	{
		blkmincorr		= 0.15;
		blknomcorr		= 0.40;
		xyconf			= 0.75;
		xml_hires[0]	= 0;
		zmin			= 0;
		zmax			= 32768;
		blksize			= 10;
		abscl			= 50;
		ablgord			= 1;	// 1  probably good for Davi EM
		absdev			= 42;	// 42 useful for Davi EM
		maxDZ			= 10;
		NoFolds			= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_alnmon	gArgs;
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
		"Usage: cross_carveblocks <xml-file> -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			DskAbsPath( xml_hires, sizeof(xml_hires), argv[i], flog );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &blksize, "-b=%d", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			;
		else if( GetArg( &ablgord, "-ablgord=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &blkmincorr, "-blkmincorr=%lf", argv[i] ) )
			;
		else if( GetArg( &blknomcorr, "-blknomcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &xyconf, "-xyconf=%lf", argv[i] ) ) {

			if( xyconf < 0.0 || xyconf > 1.0 )
				xyconf = 0.5;
		}
		else if( GetArg( &maxDZ, "-maxDZ=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* WriteSubblocksFile -------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubblocksFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "subblocks.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Purpose:\n" );
	fprintf( f, "# Sixth step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > cross_thisblock [options]\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Do one block alignment job, data read from 'blockdat.txt'.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -nf\t\t\t\t;no foldmasks\n" );
	fprintf( f, "# -abscl=50\t\t\t;integer scale reduction\n" );
	fprintf( f, "# -ablgord=1\t\t;Legendre poly field-flat max int order\n" );
	fprintf( f, "# -absdev=42\t\t;int: if > 0, img normed to mean=127, sd=sdev (recmd 42)\n" );
	fprintf( f, "# -blkmincorr=0.15\t;required min corr for alignment\n" );
	fprintf( f, "# -blknomcorr=0.40\t;nominal corr for alignment\n" );
	fprintf( f, "# -xyconf=0.75\t\t;search radius = (1-xyconf)*blkwide\n" );
	fprintf( f, "# -abdbg\t\t\t;make diagnostic images and exit (Z@Z-1)\n" );
	fprintf( f, "# -abdbg=k\t\t\t;make diagnostic images and exit (Z@k)\n" );
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
	fprintf( f, "\t\t\tqsub -N x$jb-$lyr -cwd -V -b y -pe batch 8 cross_thisblock%s -abscl=%d -ablgord=%d -absdev=%d -blkmincorr=%g -blknomcorr=%g -xyconf=%g\n",
	(gArgs.NoFolds ? " -nf" : ""),
	gArgs.abscl, gArgs.ablgord, gArgs.absdev,
	gArgs.blkmincorr, gArgs.blknomcorr, gArgs.xyconf );
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
	fprintf( f, "\techo $lyr\n" );
	fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
	fprintf( f, "\tthen\n" );
	fprintf( f, "\t\tcd $lyr\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcnt=$(($cnt+1))\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n" );
	fprintf( f, "\n" );
	fprintf( f, "echo down dirs = $cnt\n" );
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
	fprintf( f, "# Seventh step in cross-layer alignment.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Run this after subblocks completes to compile tables\n" );
	fprintf( f, "# of block alignment errors, FAILs and make sizes.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > ./breports.sht\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "ls -l ../*/D*/xD*.e* > BlockErrs.txt\n" );
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
	dx = gArgs.blksize;
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

	int	lowcount	= int(gArgs.blksize * gArgs.blksize * 0.25),
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

	for( int i = 0; i < nb; ++i ) {

		int	iy = i / kx,
			ix = i - kx * iy,
			ij = K[i].vID.size();

		fprintf( flog, "%d%c", ij, (ix == kx - 1 ? '\n' : '\t') );
		ntiles += ij;
	}

	fprintf( flog, "Total = %d\n", ntiles );
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
// zmin = min( za - gArgs.maxDZ, gArgs.zmin ). Since zb is
// usually just za-1, we could write that equivalently as
// zmin = min( zb + 1 - gArgs.maxDZ, gArgs.zmin ). This
// second form is preferred in cases where layers are missing
// so that zb < za-1, but we still want to give cross_thisblock
// several tries to get a good match (yes, even if we exceed
// maxDZ to get those tries).

	int zmin = max( zb + 1 - gArgs.maxDZ, gArgs.zmin );

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

			fprintf( f, "file=%s\n", gArgs.xml_hires );

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

/* ------------- */
/* Read src file */
/* ------------- */

	TS.FillFromTrakEM2( gArgs.xml_hires, gArgs.zmin, gArgs.zmax );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.GetTileDims( gW, gH );

	TS.SortAll_z_id();

/* ------------- */
/* Driver script */
/* ------------- */

	WriteSubblocksFile();
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


