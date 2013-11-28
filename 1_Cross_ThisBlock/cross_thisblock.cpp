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

#include	<algorithm>
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

class BlkZ {
public:
	TAffine	T;
	double	R;
	int		Z,
			used;
public:
	BlkZ( TAffine T, double R, int Z )
	: T(T), R(R), Z(Z), used(0) {};
};

class Pair {
public:
	double	area;
	int		id, iz;
public:
	Pair()	{};	// for resize()
	Pair( double area, int id, int iz )
	: area(area), id(id), iz(iz) {};
};

/* --------------------------------------------------------------- */
/* Superscape ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSuperscape {

public:
	vector<int>	vID;
	DBox	bb;			// oriented bounding box
	Point	Opts;		// origin of aligned point list
	double	x0, y0;		// scape corner in oriented system
	uint8	*ras;		// scape pixels
	uint32	ws, hs;		// scape dims
	int		is0, isN;	// layer index range
	int		clbl, rsvd;	// 'A' or 'B'

public:
	CSuperscape() : ras(NULL) {};

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
	bool MakeRasB( const DBox &Abb );

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
			blkmincorr,
			blknomcorr,
			abctr,
			xyconf;		// search radius = (1-conf)(blockwide)
	int		abscl,
			ablgord,
			absdev,
			dbgz;
	bool	abdbg,
			NoFolds;

public:
	CArgs_scp()
	{
		blkmincorr	= 0.45;
		blknomcorr	= 0.50;
		abctr		= 0.0;
		xyconf		= 0.75;
		abscl		= 50;
		ablgord		= 1;	// 1  probably good for Davi EM
		absdev		= 42;	// 42 useful for Davi EM
		dbgz		= -1;
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
		else if( GetArg( &blkmincorr, "-blkmincorr=%lf", argv[i] ) )
			;
		else if( GetArg( &blknomcorr, "-blknomcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &abctr, "-abctr=%lf", argv[i] ) )
			;
		else if( GetArg( &xyconf, "-xyconf=%lf", argv[i] ) ) {

			if( xyconf < 0.0 || xyconf > 1.0 )
				xyconf = 0.5;
		}
		else if( GetArg( &dbgz, "-abdbg=%d", argv[i] ) )
			abdbg = true;
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
	bb.L =  BIGD;
	bb.R = -BIGD;
	bb.B =  BIGD;
	bb.T = -BIGD;

	for( int i = 0; i < gDat.ntil; ++i ) {

		vector<Point>	cnr( 4 );

		cnr[0] = Point( 0.0 , 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point( 0.0 , gH-1 );

		TS.vtil[vID[i]].T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			bb.L = fmin( bb.L, cnr[k].x );
			bb.R = fmax( bb.R, cnr[k].x );
			bb.B = fmin( bb.B, cnr[k].y );
			bb.T = fmax( bb.T, cnr[k].y );
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

bool CSuperscape::MakeRasB( const DBox &Abb )
{
	vector<int>	vid;
	int			W2 = gW/2, H2 = gH/2;

	for( int i = is0; i < isN; ++i ) {

		const CUTile& U = TS.vtil[i];

		Point	p( W2, H2 );
		U.T.Transform( p );

		if( p.x >= Abb.L && p.x <= Abb.R &&
			p.y >= Abb.B && p.y <= Abb.T ) {

			vid.push_back( i );
		}
	}

	if( vid.size() < 0.05 * gDat.ntil ) {

		fprintf( flog, "Low B tile count [%d] for z=%d.\n",
		vid.size(), TS.vtil[is0].z );

		return false;
	}

	ras = TS.Scape( ws, hs, x0, y0,
			vid, gArgs.inv_abscl, 1, 0,
			gArgs.ablgord, gArgs.absdev );

	if( !ras ) {

		fprintf( flog, "Empty B scape for z=%d.\n",
		TS.vtil[is0].z );

		return false;
	}

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

	v.clear();
	p.clear();

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

	DBox	ptsbb;

	BBoxFromPoints( ptsbb, p );
	Opts = Point( ptsbb.L, ptsbb.B );

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
	const vector<BlkZ>		&vZ,
	vector<vector<Pair> >	&P,
	const CSuperscape		&A,
	const CSuperscape		&B )
{
	int		iz = vZ.size() - 1;
	TAffine	Tm = vZ[iz].T;

	for( int ia = 0; ia < gDat.ntil; ++ia ) {

		int				aid = A.vID[ia];
		vector<Pair>	&p  = P[ia];
		TAffine			Ta  = Tm * TS.vtil[aid].T;

		for( int bid = B.is0; bid < B.isN; ++bid ) {

			TAffine	Tab;
			Tab.FromAToB( Ta, TS.vtil[bid].T );

			double	area = TS.ABOlap( aid, bid, &Tab );

			if( area >= kPairOlapLo )
				p.push_back( Pair( area, bid, iz ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* ThisBZ -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if changes made.
//
static bool ThisBZ(
	vector<BlkZ>			&vZ,
	vector<vector<Pair> >	&P,
	const CSuperscape		&A,
	ThmRec					&thm,
	int						&next_isN )
{
	clock_t		t0 = StartTiming();
	CSuperscape	B;
	CThmScan	S;
	CorRec		best;

	fprintf( flog, "\n--- Start B layer ----\n" );

	B.SetLabel( 'B' );
	next_isN = B.FindLayerIndices( next_isN );

	if( gArgs.abdbg ) {
		while( next_isN != -1 && TS.vtil[B.is0].z != gArgs.dbgz )
			next_isN = B.FindLayerIndices( next_isN );
	}

	if( next_isN == -1 ) {
		fprintf( flog, "$$$ Exhausted B layers $$$\n" );
		return false;
	}

	if( !B.MakeRasB( A.bb ) )
		return false;

	B.DrawRas();

	if( !B.MakePoints( thm.bv, thm.bp ) ) {
		fprintf( flog, "No B points for z=%d.\n", TS.vtil[B.is0].z );
		return false;
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
	S.SetRThresh( gArgs.blkmincorr );
	S.SetNbMaxHt( 0.99 );
	S.SetSweepConstXY( false );
	S.SetSweepPretweak( true );
	S.SetSweepNThreads( 3 );
	S.SetUseCorrR( true );
	S.SetDisc( Ox, Oy, Rx, Ry );

	if( gArgs.abdbg ) {

		S.Pretweaks( 0, gArgs.abctr, thm );
		dbgCor = true;
		S.RFromAngle( best, gArgs.abctr, thm );
	}
	else {

		S.Pretweaks( 0, 0, thm );

		if( !S.DenovoBestAngle( best, 0, 4, .2, thm, false ) ) {

			fprintf( flog, "Low corr [%g] for z=%d.\n",
			best.R, TS.vtil[B.is0].z );

			return false;
		}

		Point	Aorigin = A.Opts;

		best.T.Apply_R_Part( Aorigin );

		best.X += B.Opts.x - Aorigin.x;
		best.Y += B.Opts.y - Aorigin.y;

		best.T.SetXY( best.X, best.Y );

		fprintf( flog, "*T: [0,1,2,3,4,5] (block-block)\n" );
		fprintf( flog, "[%g,%g,%g,%g,%g,%g]\n",
		best.T.t[0], best.T.t[1], best.T.t[2],
		best.T.t[3], best.T.t[4], best.T.t[5] );
	}

	t0 = StopTiming( flog, "Corr", t0 );

	if( gArgs.abdbg )
		return false;

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

	vZ.push_back( BlkZ( best.T, best.R, TS.vtil[B.is0].z ) );

// Accumulate pairs

	FindPairs( vZ, P, A, B );

	return true;
}

/* --------------------------------------------------------------- */
/* OrderByR ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static const vector<BlkZ>	*_vZ;

static bool OrderByR( const Pair &I, const Pair &J )
{
	return ((*_vZ)[I.iz].R > (*_vZ)[J.iz].R);
}

/* --------------------------------------------------------------- */
/* BlockCoverage ------------------------------------------------- */
/* --------------------------------------------------------------- */

// For each tile, see how much area we can get using only
// the highest quality tissue (returned prime coverage),
// and including the lower quality tissue (scdry filled in).
//
static double BlockCoverage(
	double					&scdry,
	const vector<BlkZ>		&vZ,
	vector<vector<Pair> >	&P )
{
	double	prime = 0;

	scdry	= 0;
	_vZ		= &vZ;

	for( int ia = 0; ia < gDat.ntil; ++ia ) {

		vector<Pair>	&p = P[ia];
		double			sum = 0.0;
		int				iblast = -1, nb = p.size();

		sort( p.begin(), p.end(), OrderByR );

		// first use only best (prime)

		for( int ib = 0; ib < nb; ++ib ) {

			if( vZ[p[ib].iz].R >= gArgs.blknomcorr )
				sum += p[iblast = ib].area;
			else
				break;
		}

		if( sum >= kTileAnchorLo ) {
			++prime;
			++scdry;
		}
		else {	// try adding lesser quality (scdry)

			for( int ib = iblast + 1; ib < nb; ++ib )
				sum += p[ib].area;

			if( sum >= kTileAnchorLo )
				++scdry;
		}
	}

	scdry /= gDat.ntil;

	return prime / gDat.ntil;
}

/* --------------------------------------------------------------- */
/* ZSeen --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return the fraction of maximal layer range that we've seen.
//
static double ZSeen( const vector<BlkZ> &vZ )
{
	double	nseen = vZ.size();

// do at least two
	if( nseen < 2 )
		return 0.0;

// actual fraction
	const vector<CUTile>	&vU = TS.vtil;

	return (nseen / (vU[vU.size()-1].z - vU[0].z));
}

/* --------------------------------------------------------------- */
/* KeepBest ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void KeepBest(
	vector<BlkZ>			&vZ,
	vector<vector<Pair> >	&P )
{
	for( int ia = 0; ia < gDat.ntil; ++ia ) {

		vector<Pair>	&p = P[ia];
		double			sum = 0.0;
		int				nb = p.size(),
						ng = 0;	// count good b-tiles

		// cover tile with highest quality

		do {

			if( ng >= nb )
				break;

			int	curz = p[ng].iz;

			for( int ib = ng; ib < nb; ++ib ) {

				if( p[ib].iz != curz )
					break;

				sum += p[ib].area;
				++ng;

				vZ[curz].used = 1;
			}

		} while( sum < kTileAnchorHi );

		// kill the rest

		if( ng < nb )
			p.resize( ng );
	}
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P).
//
static void WriteMakeFile(
	const CSuperscape			&A,
	const vector<BlkZ>			&vZ,
	const vector<vector<Pair> >	&P )
{
	FILE	*f;
	int		nz = vZ.size(), nused = 0;

	for( int i = 0; i < nz; ++i )
		nused += vZ[i].used;

	if( !nused ) {
		fprintf( flog, "FAIL: No tiles matched.\n" );
		return;
	}

// open the file

	f = FileOpenOrDie( "make.down", "w", flog );

// write 'all' targets line

	fprintf( f, "all: " );

	for( int ia = 0; ia < gDat.ntil; ++ia ) {

		const CUTile&		a  = TS.vtil[A.vID[ia]];
		const vector<Pair>	&p = P[ia];
		int					nb = p.size();

		for( int ib = 0; ib < nb; ++ib ) {

			const CUTile&	b = TS.vtil[p[ib].id];

			fprintf( f, "%d/%d.%d.map.tif ", a.id, b.z, b.id );
		}
	}

	fprintf( f, "\n\n" );

// Write each 'target: dependencies' line
//		and each 'rule' line

	const char	*option_nf = (gArgs.NoFolds ? " -nf" : "");

	for( int ia = 0; ia < gDat.ntil; ++ia ) {

		const CUTile&		a  = TS.vtil[A.vID[ia]];
		const vector<Pair>	&p = P[ia];
		int					nb = p.size();

		for( int ib = 0; ib < nb; ++ib ) {

			const CUTile&	b = TS.vtil[p[ib].id];
			TAffine			T;

			fprintf( f,
			"%d/%d.%d.map.tif:\n",
			a.id, b.z, b.id );

			T = vZ[p[ib].iz].T * a.T;
			T.FromAToB( T, b.T );

			fprintf( f,
			"\tptest %d/%d@%d/%d"
			" -Tab=%g,%g,%g,%g,%g,%g%s ${EXTRA}\n\n",
			a.z, a.id, b.z, b.id,
			T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
			option_nf );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteThumbFiles ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteThumbFiles( const vector<BlkZ> &vZ )
{
	int	nz = vZ.size();

	for( int iz = 0; iz < nz; ++iz ) {

		if( !vZ[iz].used )
			continue;

		char	name[128];
		sprintf( name, "ThmPair_%d_@_%d.txt", gDat.za, vZ[iz].Z );
		FILE	*f = FileOpenOrDie( name, "w", flog );
		WriteThmPairHdr( f );
		fclose( f );
	}
}

/* --------------------------------------------------------------- */
/* LayerLoop ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Here's the strategy:
// We paint a scape for the whole A-block and match it to the
// scape for the nearest whole B-block below. If the correlation
// (R-block) is at least blkmincorr then we use this B-block and
// pair off the tiles. Otherwise, the match is too uncertain and
// we ignore this B-block. Remember that the block match transform
// will set a search disc for the tile-tile jobs, and this will not
// be questioned at the tile level. Therefore, we need to feed very
// reliable block-block matches to the tile jobs.
//
// To do tile pairing, for each A-tile we keep a list of 'pair'
// records with these fields {
// - id of B-tile that overlaps it (iff overlap >= 0.02).
// - the fractional overlap area.
// - link (iz) to data about this B-block, especially R-block.
// }
//
// Next we determine if we have enough high quality data matched
// to our A-tiles, or we need to examine the next B-block. We have
// enough data if at least 0.90 of the A-tiles have good matches,
// or if we have already used the lowest allowed B-layer (zmin).
//
// The test of good matches for an A-tile goes like this:
// - Sort its pair records, descending R order.
// - Sum the fractional areas of records with R >= blknomcorr.
// - Call this A-tile matched if its sum is >= 0.30.
// Summing is admittedly coarse, we just sum fractional overlap
// areas without regard to orientation, and it's good enough for
// a coverage assessment.
//
// There's an additional clause to stop looking at more B-blocks
// a bit earlier, and this really comes into play if few of the
// blocks in this region are as good as blknomcorr. Let's call
// the coverage obtained using only great tissue (> blknomcorr)
// 'prime coverage' and the coverage we can get from all blocks
// done so far 'secondary coverage'. If we've already looked at
// half of the prospective blocks and the secondary coverage is
// already adequate (0.90 A-tiles matched) then let's give up.
//
// Once we've finished accumulating pair records, we proceed
// with writing the make.down job files. We do one A-tile at
// a time. Its list of pair records is already sorted by R. As
// we walk the pair list we keep a running sum of area. If the
// current area sum is < 0.80, we sum in all of the tiles with
// the current iz index (all come from the same layer) and we
// write jobs for those matches. If the sum is still below 0.80
// coverage we repeat for the next set of tiles with a poorer
// quality, and so on until coverage is >= 0.80 or the records
// are exhausted.
//
static void LayerLoop()
{
	clock_t					t0 = StartTiming();
	vector<BlkZ>			vZ;
	vector<vector<Pair> >	P( gDat.ntil );
	CSuperscape				A;
	ThmRec					thm;
	int						next_isN;

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

	double	prime = 0.0, scdry = 0.0;

	for(;;) {

		// align B-block and pair tiles

		int	changed = ThisBZ( vZ, P, A, thm, next_isN );

		// any changes to survey?

		if( gArgs.abdbg )
			return;

		if( next_isN == -1 )
			break;

		if( !changed )
			continue;

		// survey the coverage

		prime = BlockCoverage( scdry, vZ, P );

		fprintf( flog, "Block coverage prime %.2f scdry %.2f  Z %d\n",
		prime, scdry, vZ[vZ.size()-1].Z );

		if( (prime >= kBlockAnchorHi) ||
			(scdry >= kBlockAnchorHi && ZSeen( vZ ) >= 0.50) ) {

			break;
		}
	}

	fprintf( flog,
	"Final coverage prime %.2f scdry %.2f ntiles %d\n",
	prime, scdry, gDat.ntil );

	fprintf( flog, "\n--- Write Files ----\n" );

	KeepBest( vZ, P );
	WriteMakeFile( A, vZ, P );
	WriteThumbFiles( vZ );
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

	if( gArgs.abdbg ) {

		if( gArgs.dbgz == -1 )
			gArgs.dbgz = gDat.za - 1;

		gDat.zmin = gArgs.dbgz;
	}

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


