//
// Modeled after scapeops, align down subblocks.
//


#include	"Cmdline.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"CThmScan.h"
#include	"Geometry.h"
#include	"Maths.h"
#include	"ImageIO.h"
#include	"Timer.h"
#include	"Debug.h"

#include	<string.h>

#include	<set>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	kPairOlapLo		0.02
#define	kTileAnchorLo	0.30
#define	kTileAnchorHi	0.80
#define	kBlockAnchorHi	0.90

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Pair {

public:
	int	a, b, it;

public:
	Pair( int _a, int _b, int _it )	{a=_a; b=_b; it=_it;};
};

/* --------------------------------------------------------------- */
/* Superscape ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSuperscape {

public:
	vector<int>	vID;
	DBox	B;			// oriented bounding box
	Point	Opts;		// origin of aligned point list
	double	x0, y0;		// scape corner in oriented system
	uint8	*ras;		// scape pixels
	uint32	ws, hs;		// scape dims
	int		is0, isN;	// layer index range
	int		clbl, rsvd;	// 'A' or 'B'

public:
	CSuperscape()
		{ras = NULL;};

	virtual ~CSuperscape()
		{KillRas();};

	void SetLabel( int clbl )
		{this->clbl = clbl;};

	void KillRas()
		{
			if( ras ) {
				RasterFree( ras );
				ras = NULL;
			}
		};

	int  FindLayerIndices( int next_isN );
	void vID_From_sID();
	void CalcBBox();

	bool MakeRasA();
	bool MakeRasB( const DBox &A );

	void DrawRas();
	bool Load( FILE* flog );

	bool MakePoints( vector<double> &v, vector<Point> &p );
	void WriteMeta();
};

/* --------------------------------------------------------------- */
/* CBlockDat ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CBlockDat {

public:
	char		xmlfile[2048];
	int			za, zmin;
	int			ntil;
	set<int>	sID;

public:
	void ReadFile();
};

/* --------------------------------------------------------------- */
/* CArgs_scp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scp {

public:
	double	inv_abscl,
			abcorr,
			abctr,
			xyconf;		// search radius = (1-conf)(blockwide)
	int		abscl,
			ablgord,
			absdev;
	bool	abdbg,
			NoFolds;

public:
	CArgs_scp()
	{
		abcorr		= 0.20;
		abctr		= 0.0;
		xyconf		= 0.50;
		abscl		= 200;
		ablgord		= 1;	// 3  probably good for Davi EM
		absdev		= 0;	// 42 useful for Davi EM
		abdbg		= false;
		NoFolds		= false;

		inv_abscl	= 1.0/abscl;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scp	gArgs;
static CBlockDat	gDat;
static CTileSet		TS;
static FILE*		flog	= NULL;
static int			gW		= 0,	// universal pic dims
					gH		= 0;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scp::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "cross_thisblock.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Align this block: %s ", atime );

// parse command line args

	if( argc < 1 ) {
		printf(
		"Usage: cross_thisblock [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			inv_abscl = 1.0/abscl;
		else if( GetArg( &ablgord, "-ablgord=%d", argv[i] ) )
			;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &abctr, "-abctr=%lf", argv[i] ) )
			;
		else if( GetArg( &xyconf, "-xyconf=%lf", argv[i] ) ) {

			if( xyconf < 0.0 || xyconf > 1.0 )
				xyconf = 0.5;
		}
		else if( IsArg( "-abdbg", argv[i] ) )
			abdbg = true;
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
/* ReadFile ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CBlockDat::ReadFile()
{
	FILE		*f = FileOpenOrDie( "blockdat.txt", "r", flog );
	CLineScan	LS;
	int			item = 0;

	while( LS.Get( f ) > 0 ) {

		if( item == 0 ) {
			sscanf( LS.line, "file=%s", xmlfile );
			item = 1;
		}
		else if( item == 1 ) {
			sscanf( LS.line, "ZaZmin=%d,%d", &za, &zmin );
			item = 2;
		}
		else if( item == 2 ) {
			sscanf( LS.line, "nIDs=%d", &ntil );
			item = 3;
		}
		else {

			for( int i = 0; i < ntil; ++i ) {

				int	id;

				sscanf( LS.line, "%d", &id );
				sID.insert( id );
				LS.Get( f );
			}
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* FindLayerIndices ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Find is0 such that [is0,isN) is whole layer.
//
// Return is0.
//
int CSuperscape::FindLayerIndices( int next_isN )
{
	TS.GetLayerLimitsR( is0, isN = next_isN );
	return is0;
}

/* --------------------------------------------------------------- */
/* vID_From_sID -------------------------------------------------- */
/* --------------------------------------------------------------- */

// BlockDat.sID are actual tile-IDs, which are a robust way
// to specify <which tiles> across programs.
//
// vID are TS.vtil[] index numbers (like is0) which are an
// efficient way to walk TS data.
//
void CSuperscape::vID_From_sID()
{
	int	n = 0;

	vID.resize( gDat.ntil );

	for( int i = is0; i < isN; ++i ) {

		if( gDat.sID.find( TS.vtil[i].id ) != gDat.sID.end() ) {

			vID[n++] = i;

			if( n == gDat.ntil )
				break;
		}
	}
}

/* --------------------------------------------------------------- */
/* CalcBBox ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CSuperscape::CalcBBox()
{
	B.L =  BIGD;
	B.R = -BIGD;
	B.B =  BIGD;
	B.T = -BIGD;

	for( int i = 0; i < gDat.ntil; ++i ) {

		vector<Point>	cnr( 4 );

		cnr[0] = Point( 0.0 , 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point( 0.0 , gH-1 );

		TS.vtil[vID[i]].T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			B.L = fmin( B.L, cnr[k].x );
			B.R = fmax( B.R, cnr[k].x );
			B.B = fmin( B.B, cnr[k].y );
			B.T = fmax( B.T, cnr[k].y );
		}
	}
}

/* --------------------------------------------------------------- */
/* MakeRasA ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasA()
{
	ras = TS.Scape( ws, hs, x0, y0,
			vID, gArgs.inv_abscl, 1, 0,
			gArgs.ablgord, gArgs.absdev );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasB( const DBox &A )
{
	vector<int>	vid;
	int			W2 = gW/2, H2 = gH/2;

	for( int i = is0; i < isN; ++i ) {

		const CUTile& U = TS.vtil[i];

		Point	p( W2, H2 );
		U.T.Transform( p );

		if( p.x >= A.L && p.x <= A.R && p.y >= A.B && p.y <= A.T )
			vid.push_back( i );
	}

	ras = TS.Scape( ws, hs, x0, y0,
			vid, gArgs.inv_abscl, 1, 0,
			gArgs.ablgord, gArgs.absdev );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* DrawRas ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::DrawRas()
{
	if( ras ) {
		char	name[128];
		sprintf( name, "Ras_%c_%d.png", clbl, TS.vtil[is0].z );
		Raster8ToPng8( name, ras, ws, hs );
	}
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CSuperscape::Load( FILE* flog )
{
	char	name[128];

	sprintf( name, "Ras_%c_%d.png", clbl, TS.vtil[is0].z );
	x0 = y0 = 0.0;

	ras	= Raster8FromAny( name, ws, hs, flog );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakePoints ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CSuperscape::MakePoints( vector<double> &v, vector<Point> &p )
{
// collect point and value lists

	int	np = ws * hs, ok = true;

	for( int i = 0; i < np; ++i ) {

		if( ras[i] ) {

			int	iy = i / ws,
				ix = i - ws * iy;

			v.push_back( ras[i] );
			p.push_back( Point( ix, iy ) );
		}
	}

	if( !(np = p.size()) ) {
		ok = false;
		goto exit;
	}

// get points origin and translate to zero

	DBox	bb;

	BBoxFromPoints( bb, p );
	Opts = Point( bb.L, bb.B );

	for( int i = 0; i < np; ++i ) {

		p[i].x -= Opts.x;
		p[i].y -= Opts.y;
	}

// normalize values

	{
		double	sd = Normalize( v );

		if( !sd || !isfinite( sd ) )
			ok = false;
	}

exit:
	KillRas();

	return ok;
}

/* --------------------------------------------------------------- */
/* WriteMeta ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::WriteMeta()
{
	fprintf( flog,
	"*%c: z scl [ws,hs] [x0,y0]\n", clbl );

	fprintf( flog,
	"%d %d [%d,%d] [%g,%g]\n",
	TS.vtil[is0].z, gArgs.abscl, ws, hs, x0, y0 );
}

/* --------------------------------------------------------------- */
/* FindPairs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void FindPairs(
	vector<double>		&Asum,
	vector<TAffine>		&vT,
	vector<Pair>		&P,
	const CSuperscape	&A,
	const CSuperscape	&B )
{
	int		it = vT.size() - 1;
	TAffine	Tm = vT[it];

	for( int ka = 0; ka < gDat.ntil; ++ka ) {

		int	ia = A.vID[ka];

		if( Asum[ka] >= kTileAnchorHi )
			continue;

		TAffine Ta = Tm * TS.vtil[ia].T;

		for( int ib = B.is0; ib < B.isN; ++ib ) {

			TAffine	Tab;
			Tab.FromAToB( Ta, TS.vtil[ib].T );

			double	area = TS.ABOlap( ia, ib, &Tab );

			if( area >= kPairOlapLo ) {

				P.push_back( Pair( ia, ib, it ) );
				Asum[ka] += area;
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* ThisBZ -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ThisBZ(
	vector<double>		&Asum,
	vector<TAffine>		&vT,
	vector<Pair>		&P,
	const CSuperscape	&A,
	ThmRec				&thm,
	int					&next_isN )
{
	clock_t		t0 = StartTiming();
	CSuperscape	B;
	CThmScan	S;
	CorRec		best;

	fprintf( flog, "\n--- Start B layer ----\n" );

	B.SetLabel( 'B' );
	next_isN = B.FindLayerIndices( next_isN );

	if( next_isN == -1 ) {
		fprintf( flog, "$$$ Exhausted B layers $$$\n" );
		return;
	}

	if( !B.MakeRasB( A.B ) ) {
		fprintf( flog, "No B tiles for z=%d.\n", TS.vtil[B.is0].z );
		return;
	}

	B.DrawRas();

	if( !B.MakePoints( thm.bv, thm.bp ) ) {
		fprintf( flog, "No B points for z=%d.\n", TS.vtil[B.is0].z );
		return;
	}

	B.WriteMeta();
	t0 = StopTiming( flog, "MakeRasB", t0 );

	thm.ftc.clear();
	thm.reqArea	= int(kPairOlapLo * A.ws * A.hs);
	thm.olap1D	= 4;
	thm.scl		= 1;

	int	Ox	= int(A.x0 - B.x0),
		Oy	= int(A.y0 - B.y0),
		Rx	= int((1.0 - gArgs.xyconf) * A.ws),
		Ry	= int((1.0 - gArgs.xyconf) * A.hs);

	S.Initialize( flog, best );
	S.SetRThresh( gArgs.abcorr );
	S.SetNbMaxHt( 0.99 );
	S.SetSweepConstXY( false );
	S.SetSweepPretweak( true );
	S.SetUseCorrR( true );
	S.SetDisc( Ox, Oy, Rx, Ry );

	if( gArgs.abdbg ) {

		S.Pretweaks( 0, gArgs.abctr, thm );
		dbgCor = true;
		S.RFromAngle( best, gArgs.abctr, thm );
	}
	else {

		S.Pretweaks( 0, 0, thm );

		if( S.DenovoBestAngle( best, 0, 4, .2, thm ) ) {

			Point	Aorigin = A.Opts;

			best.T.Apply_R_Part( Aorigin );

			best.X += B.Opts.x - Aorigin.x;
			best.Y += B.Opts.y - Aorigin.y;

			best.T.SetXY( best.X, best.Y );
		}
		else {
			// return a block-block transform that
			// converts to identity montage-montage
			best.T.NUSetOne();
			best.T.SetXY( A.x0 - B.x0, A.y0 - B.y0 );
		}

		fprintf( flog, "*T: [0,1,2,3,4,5] (block-block)\n" );
		fprintf( flog, "[%g,%g,%g,%g,%g,%g]\n",
		best.T.t[0], best.T.t[1], best.T.t[2],
		best.T.t[3], best.T.t[4], best.T.t[5] );
	}

	t0 = StopTiming( flog, "Corr", t0 );

	if( gArgs.abdbg )
		return;

// Build: montage -> montage transform

	TAffine	s, t;

	// A montage -> A block image
	s.NUSetScl( gArgs.inv_abscl );
	s.AddXY( -A.x0, -A.y0 );

	// A block image -> B block image
	t = best.T * s;

	// B block image -> B montage
	t.AddXY( B.x0, B.y0 );
	s.NUSetScl( gArgs.abscl );
	best.T = s * t;

// Append to list

	vT.push_back( best.T );

// Accumulate pairs

	FindPairs( Asum, vT, P, A, B );
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P).
//
static void WriteMakeFile(
	const vector<TAffine>	&vT,
	const vector<Pair>		&P )
{
	FILE	*f;
	int		np = P.size();

	if( !np ) {
		fprintf( flog, "FAIL: No tiles matched.\n" );
		return;
	}

// open the file

	f = FileOpenOrDie( "make.down", "w", flog );

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
		TAffine			T;

		fprintf( f,
		"%d/%d.%d.map.tif:\n",
		A.id, B.z, B.id );

		T = vT[P[i].it] * A.T;
		T.FromAToB( T, B.T );

		fprintf( f,
		"\tptest %d/%d@%d/%d -Tab=%g,%g,%g,%g,%g,%g%s ${EXTRA}\n\n",
		A.z, A.id, B.z, B.id,
		T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
		option_nf );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteThumbFiles ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteThumbFiles( const vector<Pair> &P )
{
	int	lastZ = -1, np = P.size();

	for( int i = 0; i < np; ++i ) {

		int	z = TS.vtil[P[i].b].z;

		if( z != lastZ ) {

			char	name[128];
			sprintf( name, "ThmPair_%d_@_%d.txt", gDat.za, z );
			FILE	*f = FileOpenOrDie( name, "w", flog );
			WriteThmPairHdr( f );
			fclose( f );

			lastZ = z;
		}
	}
}

/* --------------------------------------------------------------- */
/* LayerLoop ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Here's the strategy:
// We paint a scape for the whole A-block and match it to the
// nearest whole B-block below, and then we pair off the tiles.
//
// An A-tile is paired to a B-tile if their overlap is >= 0.02.
// Each A-tile keeps a running sum of these pairwise overlaps.
// The summing of area is deliberately coarse: we just sum it
// up over different B-tiles without worrying about orientation.
// An A-tile is very well anchored if its area sum is > 0.80 and
// is adequately anchored if sum > 0.30.
//
// Next we survey the results of the pair-ups by counting how many
// A-tiles are adequately anchored. If the count is below 0.90 of
// the total number of A-block tiles, we repeat the above process
// by matching the A-block to the next lower B-block in the stack,
// doing the tile-wise pairing, and accumulating overlap areas into
// the existing A-tile sums.
//
// That is repeated until the block is adequately anchored
// or we have already used the lowest allowed B-layer (zmin).
//
static void LayerLoop()
{
	clock_t			t0 = StartTiming();
	vector<double>	Asum( gDat.ntil, 0.0 );
	vector<TAffine>	vT;
	vector<Pair>	P;
	CSuperscape		A;
	ThmRec			thm;
	int				next_isN;

	fprintf( flog, "\n--- Start A layer ----\n" );

	A.SetLabel( 'A' );
	next_isN = A.FindLayerIndices( TS.vtil.size() );
	A.vID_From_sID();
	A.CalcBBox();
	A.MakeRasA();
	A.DrawRas();
	A.MakePoints( thm.av, thm.ap );
	A.WriteMeta();
	t0 = StopTiming( flog, "MakeRasA", t0 );

	for(;;) {

		// align B-block and pair tiles

		int	psize = P.size();

		ThisBZ( Asum, vT, P, A, thm, next_isN );

		// any changes to survey?

		if( next_isN == -1 )
			break;

		if( P.size() <= psize )
			continue;

		// survey the coverage

		int	ngood = 0;

		for( int i = 0; i < gDat.ntil; ++i ) {

			if( Asum[i] >= kTileAnchorLo )
				++ngood;
		}

		fprintf( flog, "Block coverage %.2f\n",
			(double)ngood/gDat.ntil );

		if( ngood >= kBlockAnchorHi * gDat.ntil )
			break;
	}

	fprintf( flog, "\n--- Write Files ----\n" );

	WriteMakeFile( vT, P );
	WriteThumbFiles( P );
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

/* --------------- */
/* Read block data */
/* --------------- */

	gDat.ReadFile();

/* ---------------- */
/* Read source data */
/* ---------------- */

	TS.FillFromTrakEM2( gDat.xmlfile, gDat.zmin, gDat.za );

	fprintf( flog, "Got %d images.\n", TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.GetTileDims( gW, gH );

	t0 = StopTiming( flog, "ReadFile", t0 );

/* ------------- */
/* Sort by layer */
/* ------------- */

	TS.SortAll_z();

/* ----- */
/* Stuff */
/* ----- */

	LayerLoop();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


