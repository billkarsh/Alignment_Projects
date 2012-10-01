//
// Carve cross layer work into blocks; write subblocks.sht.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"Geometry.h"


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
	double	abcorr;
	char	xml_hires[2048];
	char	*pat;
	int		zmin,
			zmax,
			blksize,
			abscl,
			absdev;

public:
	CArgs_alnmon()
	{
		abcorr			= 0.20;
		xml_hires[0]	= 0;
		pat				= "/N";
		zmin			= 0;
		zmax			= 32768;
		blksize			= 10;
		abscl			= 200;
		absdev			= 0;	// 12 useful for Davi EM
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
		else if( GetArgStr( pat, "-p", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArg( &blksize, "-b=%d", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
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
/* WriteSubblocksFile -------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubblocksFile()
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "../subblocks.sht" );
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
	fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
	fprintf( f, "\t\tdo\n" );
	fprintf( f, "\t\t\tcd $jb\n" );
	fprintf( f, "\t\t\tqsub -N q$jb-$lyr -cwd -V -b y -pe batch 8"
		" cross_thisblock -p%s -abscl=%d -absdev=%d -abcorr=%g\n",
		gArgs.pat, gArgs.abscl, gArgs.absdev, gArgs.abcorr );
	fprintf( f, "\t\t\tcd ..\n" );
	fprintf( f, "\t\tdone\n" );
	fprintf( f, "\t\tcd ..\n" );
	fprintf( f, "\tfi\n" );
	fprintf( f, "done\n\n" );

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
	for( int i = 0; i < nb; ++i ) {

		int	nID;

		if( nID = K[i].vID.size() ) {

			char	path[256];
			int		iy = i / kx,
					ix = i - kx * iy;

			// make Dx_y folder
			sprintf( path, "../%d", za );
			CreateJobsDir( path, ix, iy, za, zb, flog );

			// write params file
			sprintf( path, "../%d/D%d_%d/blockdat.txt", za, ix, iy );
			FILE	*f = FileOpenOrDie( path, "w", flog );

			fprintf( f, "file=%s\n", gArgs.xml_hires );
			fprintf( f, "ZaZb=%d,%d\n", za, zb );

			fprintf( f, "nIDs=%d\n", nID );
			for( int j = 0; j < nID; ++j )
				fprintf( f, "%d\n", K[i].vID[j] );

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
	TS.SetDecoderPat( gArgs.pat );

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


