//
// Modeled after scapeops, align down subblocks.
//


#include	"Cmdline.h"
#include	"File.h"
#include	"CTileSet.h"
#include	"Scape.h"
#include	"CThmScan.h"
#include	"Geometry.h"
#include	"Maths.h"
#include	"ImageIO.h"
#include	"Timer.h"
#include	"Debug.h"

#include	<set>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

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
	void vID_From_sID();
	void CalcBBox();

	bool MakeRasA();
	bool MakeRasB( const DBox &A );

	void WriteMeta( char clbl, int z, int scl );

	void MakePoints( vector<double> &v, vector<Point> &p );
};

/* --------------------------------------------------------------- */
/* CBlockDat ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CBlockDat {

public:
	char		xmlfile[2048];
	int			za, zb;
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
			abctr;
	char	*pat;
	int		abscl,
			absdev;
	bool	abdbg,
			NoFolds;

public:
	CArgs_scp()
	{
		abcorr		= 0.20;
		abctr		= 0.0;
		pat			= "/N";
		abscl		= 200;
		absdev		= 0;	// 12 useful for Davi EM
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

		if( GetArgStr( pat, "-p=", argv[i] ) )
			;
		else if( GetArg( &abscl, "-abscl=%d", argv[i] ) )
			inv_abscl = 1.0/abscl;
		else if( GetArg( &absdev, "-absdev=%d", argv[i] ) )
			;
		else if( GetArg( &abcorr, "-abcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &abctr, "-abctr=%lf", argv[i] ) )
			;
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
			sscanf( LS.line, "ZaZb=%d,%d", &za, &zb );
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

void CSuperscape::FindLayerIndices( int z )
{
	TS.GetLayerLimits( is0 = 0, isN );

	while( isN != -1 && TS.vtil[is0].z != z )
		TS.GetLayerLimits( is0 = isN, isN );
}

/* --------------------------------------------------------------- */
/* vID_From_sID -------------------------------------------------- */
/* --------------------------------------------------------------- */

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
	vector<ScpTile>	S;

	for( int i = 0; i < gDat.ntil; ++i ) {

		const CUTile& U = TS.vtil[vID[i]];

		S.push_back( ScpTile( U.name, U.T ) );
	}

	ras = Scape( ws, hs, x0, y0, S, gW, gH,
			gArgs.inv_abscl, 1, 0, gArgs.absdev, flog );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasB( const DBox &A )
{
	vector<ScpTile>	S;
	int				W2 = gW/2, H2 = gH/2;

	for( int i = is0; i < isN; ++i ) {

		const CUTile& U = TS.vtil[i];

		Point	p( W2, H2 );
		U.T.Transform( p );

		if( p.x >= A.L && p.x <= A.R && p.y >= A.B && p.y <= A.T )
			S.push_back( ScpTile( U.name, U.T ) );
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
	"*%c: z scl [ws,hs] [x0,y0]\n", clbl );

	fprintf( flog,
	"%d %d [%d,%d] [%g,%g]\n",
	z, scl, ws, hs, x0, y0 );
}

/* --------------------------------------------------------------- */
/* MakePoints ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CSuperscape::MakePoints( vector<double> &v, vector<Point> &p )
{
// collect point and value lists

	int		np = ws * hs;

	for( int i = 0; i < np; ++i ) {

		if( ras[i] ) {

			int	iy = i / ws,
				ix = i - ws * iy;

			v.push_back( ras[i] );
			p.push_back( Point( ix, iy ) );
		}
	}

// get points origin and translate to zero

	DBox	bb;

	BBoxFromPoints( bb, p );
	Opts = Point( bb.L, bb.B );
	np = p.size();

	for( int i = 0; i < np; ++i ) {

		p[i].x -= Opts.x;
		p[i].y -= Opts.y;
	}

// normalize values

	double	sd = Normalize( v );

	if( !sd || !isfinite( sd ) ) {
		fprintf( flog, "FAIL: Block has stdev: %f\n", sd );
		exit( 42 );
	}

	KillRas();
}

/* --------------------------------------------------------------- */
/* FindPairs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void FindPairs(
	vector<Pair>		&P,
	const CSuperscape	&A,
	const CSuperscape	&B,
	const TForm			&Tm )
{
	int	na = A.vID.size();

	for( int ka = 0; ka < na; ++ka ) {

		TForm	Ta;
		int		ia = A.vID[ka];

		MultiplyTrans( Ta, Tm, TS.vtil[ia].T );

		for( int ib = B.is0; ib < B.isN; ++ib ) {

			TForm	Tab;

			AToBTrans( Tab, Ta, TS.vtil[ib].T );

			if( TS.ABOlap( ia, ib, &Tab ) >= 0.02 )
				P.push_back( Pair( ia, ib ) );
		}
	}
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P).
//
static void WriteMakeFile(
	const vector<Pair>	&P,
	const TForm			&Tm )
{
	FILE	*f;
	int		np = P.size();

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
		TForm			T;

		fprintf( f,
		"%d/%d.%d.map.tif:\n",
		A.id, B.z, B.id );

		MultiplyTrans( T, Tm, A.T );
		AToBTrans( T, TForm( T ), B.T );

		fprintf( f,
		"\tptest %d/%d@%d/%d -Tab=%g,%g,%g,%g,%g,%g%s ${EXTRA}\n\n",
		A.z, A.id, B.z, B.id,
		T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
		option_nf );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ScapeStuff ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScapeStuff()
{
	clock_t		t0 = StartTiming();
	CSuperscape	A, B;
	ThmRec		thm;
	CThmScan	S;
	CorRec		best;

	A.FindLayerIndices( gDat.za );
	A.vID_From_sID();
	A.CalcBBox();
#if 1
	A.MakeRasA();
	A.DrawRas( "Aras.png" );
#else
	A.Load( "Aras.png", flog );
#endif
	t0 = StopTiming( flog, "MakeRasA", t0 );

	B.FindLayerIndices( gDat.zb );
#if 1
	if( !B.MakeRasB( A.B ) ) {
		fprintf( flog, "FAIL: No B tiles in block.\n" );
		return;
	}
	B.DrawRas( "Bras.png" );
#else
	B.Load( "Bras.png", flog );
#endif
	t0 = StopTiming( flog, "MakeRasB", t0 );

	A.MakePoints( thm.av, thm.ap );
	A.WriteMeta( 'A', gDat.za, gArgs.abscl );

	B.MakePoints( thm.bv, thm.bp );
	B.WriteMeta( 'B', gDat.zb, gArgs.abscl );

	thm.ftc.clear();
	thm.reqArea	= 0;
	thm.olap1D	= 0;
	thm.scl		= 1;

	int	Ox	= int(A.x0 - B.x0),
		Oy	= int(A.y0 - B.y0),
		Rx	= int(0.33 * A.ws),
		Ry	= int(0.33 * A.hs);

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

			best.T.Apply_R_Part( A.Opts );

			best.X += B.Opts.x - A.Opts.x;
			best.Y += B.Opts.y - A.Opts.y;

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

	TForm	s, t;

	// A montage -> A block image
	s.NUSetScl( gArgs.inv_abscl );
	s.AddXY( -A.x0, -A.y0 );

	// A block image -> B block image
	MultiplyTrans( t, best.T, s );

	// B block image -> B montage
	t.AddXY( B.x0, B.y0 );
	s.NUSetScl( gArgs.abscl );
	MultiplyTrans( best.T, s, t );

// List pairs

	vector<Pair>	P;

	FindPairs( P, A, B, best.T );

// Write make file for block

	WriteMakeFile( P, best.T );
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

/* --------------- */
/* Read block data */
/* --------------- */

	gDat.ReadFile();

/* ---------------- */
/* Read source data */
/* ---------------- */

	TS.FillFromTrakEM2( gDat.xmlfile, gDat.zb, gDat.za );

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

	ScapeStuff();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


