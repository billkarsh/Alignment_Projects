//
// scapeops
//
// Given inpath -idb=idbpath...
//
// Perform montage drawing and/or strip aligning as follows:
//
//	If drawing a montage...
//
//	-mb -zb=%d -scl=%d
//
//		[-lgord=%d] [-sdev=%d] [-resmask]
//
// If aligning strips...
//
//	-ab -za=%d -zb=%d -scl=%d -abwide=%d -stpcorr=%lf
//
//		[-lgord=%d] [-sdev=%d] [-resmask] [-abdbg] [-abctr=%lf]
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
#include	"Memory.h"
#include	"Debug.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Superscape ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSuperscape {

public:
	DBox	B;			// oriented bounding box
	Point	Opts;		// origin of aligned point list
	double	x0, y0;		// scape corner in oriented system
	uint8	*ras;		// scape pixels
	uint32	ws, hs;		// scape dims
	int		is0, isN,	// layer index range
			Bxc, Byc,	// oriented layer center
			Bxw, Byh,	// oriented layer span
			deg;		// rotate this much to orient

public:
	CSuperscape() : ras(NULL) {};

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
	double	inv_scl,
			stpcorr,
			abctr;
	char	*inpath;
	string	idb;
	int		za,
			zb,
			scl,
			lgord,
			sdev,
			abwide;
	bool	ismb,
			isab,
			resmask,
			abdbg;

public:
	CArgs_scp()
	{
		stpcorr		= 0.02;
		abctr		= 0.0;
		inpath		= NULL;
		za			= -1;
		zb			= -1;
		scl			= 50;
		lgord		= 1;	// 1  probably good for Davi EM
		sdev		= 42;	// 42 useful for Davi EM
		abwide		= 15;
		ismb		= false;
		isab		= false;
		abdbg		= false;

		inv_scl	= 1.0/scl;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_scp	gArgs;
static CTileSet		TS;
static CSuperscape	*gA, *gB;
static CThmScan		*gS;
static FILE*		flog	= NULL;
static int			gW		= 0,	// universal pic dims
					gH		= 0;






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

	const char	*pchar;

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			inpath = argv[i];
		else if( GetArgStr( pchar, "-idb=", argv[i] ) )
			idb = pchar;
		else if( GetArg( &za, "-za=%d", argv[i] ) )
			;
		else if( GetArg( &zb, "-zb=%d", argv[i] ) )
			;
		else if( GetArg( &scl, "-scl=%d", argv[i] ) )
			inv_scl = 1.0/scl;
		else if( GetArg( &lgord, "-lgord=%d", argv[i] ) )
			;
		else if( GetArg( &sdev, "-sdev=%d", argv[i] ) )
			;
		else if( GetArg( &abwide, "-abwide=%d", argv[i] ) )
			;
		else if( GetArg( &stpcorr, "-stpcorr=%lf", argv[i] ) )
			;
		else if( GetArg( &abctr, "-abctr=%lf", argv[i] ) )
			;
		else if( IsArg( "-mb", argv[i] ) )
			ismb = true;
		else if( IsArg( "-ab", argv[i] ) )
			isab = true;
		else if( IsArg( "-resmask", argv[i] ) )
			resmask = true;
		else if( IsArg( "-abdbg", argv[i] ) )
			abdbg = true;
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

	deg = TightestBBox( B, C );

	if( B.R - B.L > B.T - B.B ) {

		// Make tall: rotate 90 degrees if wider than tall

		deg = (deg < 0 ? deg + 90 : deg - 90);
		R.NUSetRot( deg*PI/180 );
		R.Apply_R_Part( C );
		BBoxFromPoints( B, C );
	}

	R.NUSetRot( deg*PI/180 );

	for( int i = is0; i < isN; ++i ) {

		TAffine&	T = TS.vtil[i].T;

		T = R * T;
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
	vector<int>	vid( isN - is0 );

	for( int i = is0; i < isN; ++i )
		vid[i - is0] = i;

	ras = TS.Scape( ws, hs, x0, y0,
			vid, 1.0/gArgs.scl, 1, 0,
			gArgs.lgord, gArgs.sdev,
			gArgs.resmask );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasV ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasV()
{
// Collect strip tiles

	vector<int>	vid;
	int			w1, w2, h1, h2;

	w1 = (gArgs.abwide * gW)/2;
	w2 = Bxc + w1;
	w1 = Bxc - w1;

	h1 = int(Byh * 0.45);
	h2 = Byc + h1;
	h1 = Byc - h1;

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		if( U.T.t[2] + gW > w1 && U.T.t[2] < w2 &&
			U.T.t[5] + gH > h1 && U.T.t[5] < h2 ) {

			vid.push_back( i );
		}
	}

	ras = TS.Scape( ws, hs, x0, y0,
			vid, gArgs.inv_scl, 1, 0,
			gArgs.lgord, gArgs.sdev,
			gArgs.resmask );

	return (ras != NULL);
}

/* --------------------------------------------------------------- */
/* MakeRasH ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CSuperscape::MakeRasH()
{
// Collect strip tiles

	vector<int>	vid;
	int			w1, w2, h1, h2;

	w1 = int(Bxw * 0.45);
	w2 = Bxc + w1;
	w1 = Bxc - w1;

	h1 = (gArgs.abwide * gH)/2;
	h2 = Byc + h1;
	h1 = Byc - h1;

	for( int i = is0; i < isN; ++i ) {

		const CUTile&	U = TS.vtil[i];

		if( U.T.t[2] + gW > w1 && U.T.t[2] < w2 &&
			U.T.t[5] + gH > h1 && U.T.t[5] < h2 ) {

			vid.push_back( i );
		}
	}

	ras = TS.Scape( ws, hs, x0, y0,
			vid, gArgs.inv_scl, 1, 0,
			gArgs.lgord, gArgs.sdev,
			gArgs.resmask );

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

	if( !(np = p.size()) ) {
		fprintf( flog, "FAIL: Block has no non-zero pixels.\n" );
		exit( 42 );
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

	if( !Normalize( v ) ) {
		fprintf( flog, "FAIL: Scape stdev = 0.\n" );
		exit( 42 );
	}

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
/* NewAngProc ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// At any given angle the corners of A and B strip coincide
// and are at (0,0) in B-system. (Ox,Oy) brings center of A
// coincident with center of B, so, Oxy = Bc - Rot(Ac).
//
static void NewAngProc(
	int		&Ox,
	int		&Oy,
	int		&Rx,
	int		&Ry,
	double	deg )
{
	double	r  = deg*PI/180,
			c  = cos( r ),
			s  = sin( r );

	Ox = int(gB->ws - (c*gA->ws - s*gA->hs))/2;
	Oy = int(gB->hs - (s*gA->ws + c*gA->hs))/2;
	Rx = int(gA->ws * 0.25);
	Ry = int(gB->hs * 0.25);
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
//	TAffine	Rbi, Ra, t;
//
// Scale back up
//	best.T.MulXY( gArgs.scl );
//	A.x0 *= gArgs.scl;
//	A.y0 *= gArgs.scl;
//	B.x0 *= gArgs.scl;
//	B.y0 *= gArgs.scl;
//
// A-montage -> A-oriented
//	Ra.NUSetRot( A.deg*PI/180 );
//
// A-oriented -> A-scape
//	Ra.AddXY( -A.x0, -A.y0 );
//
// A-scape -> B-scape
//	t = best.T * Ra;
//
// B-scape -> B-oriented
//	t.AddXY( B.x0, B.y0 );
//
// B-oriented -> B-montage
//	Rbi.NUSetRot( -B.deg*PI/180 );
//	best.T = Rbi * t;
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
		B.WriteMeta( 'M', gArgs.zb, gArgs.scl );
		t0 = StopTiming( flog, "MakeMontage", t0 );
	}

	if( !gArgs.isab )
		return;

	A.FindLayerIndices( gArgs.za );
	A.OrientLayer();

	MakeStripRasters( A, B );
	t0 = StopTiming( flog, "MakeStrips", t0 );

	A.MakePoints( thm.av, thm.ap );
	A.WriteMeta( 'A', gArgs.za, gArgs.scl );

	B.MakePoints( thm.bv, thm.bp );
	B.WriteMeta( 'B', gArgs.zb, gArgs.scl );

	thm.ftc.clear();
	thm.reqArea	= int(gW * gH * gArgs.inv_scl * gArgs.inv_scl);
	thm.olap1D	= int(gW * gArgs.inv_scl * 0.5);
	thm.scl		= 1;

	S.Initialize( flog, best );
	S.SetRThresh( gArgs.stpcorr );
	S.SetNbMaxHt( 0.99 );
	S.SetSweepConstXY( false );
	S.SetSweepPretweak( true );
	S.SetSweepNThreads( 8 );
	S.SetUseCorrR( true );
	S.SetNewAngProc( NewAngProc );

	gA = &A;
	gB = &B;
	gS = &S;

	if( gArgs.abdbg ) {

		//S.RFromAngle( best, gArgs.abctr, thm );
		//S.Pretweaks( best.R, gArgs.abctr, thm );
		dbgCor = true;
		S.RFromAngle( best, gArgs.abctr, thm );
	}
	else {

		S.DenovoBestAngle( best, 0, 180, 5, thm, true );

		best.T.Apply_R_Part( A.Opts );

		best.X += B.Opts.x - A.Opts.x;
		best.Y += B.Opts.y - A.Opts.y;

		best.T.SetXY( best.X, best.Y );

		fprintf( flog, "*T: [0,1,2,3,4,5] (strip-strip)\n" );
		fprintf( flog, "[%f,%f,%f,%f,%f,%f]\n",
		best.T.t[0], best.T.t[1], best.T.t[2],
		best.T.t[3], best.T.t[4], best.T.t[5] );
	}

	t0 = StopTiming( flog, "Corr", t0 );
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

/* ---------------- */
/* Read source data */
/* ---------------- */

	if( gArgs.zb >= 0 && gArgs.za < 0 )
		gArgs.za = gArgs.zb;

	TS.FillFromRgns( gArgs.inpath, gArgs.idb, gArgs.zb, gArgs.za );

	fprintf( flog, "Got %d images.\n", (int)TS.vtil.size() );

	if( !TS.vtil.size() )
		goto exit;

	TS.SetTileDimsFromImageFile();
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
	VMStats( flog );
	fclose( flog );

	return 0;
}


