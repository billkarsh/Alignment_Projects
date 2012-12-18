//
// scapeops
//
// Given file, id pattern:
//
// xmlfile -p=%s	# e.g. -p=_N_
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
#include	"File.h"
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

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

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
			*pat;
	int		za,
			zb,
			abwide,
			mbscl,
			abscl,
			mbsdev,
			absdev;
	bool	ismb,
			isab,
			abdbg;

public:
	CArgs_scp()
	{
		abcorr		= 0.20;
		abctr		= 0.0;
		infile		= NULL;
		pat			= "/N";
		za			= -1;
		zb			= -1;
		abwide		= 5;
		mbscl		= 200;
		abscl		= 200;
		mbsdev		= 0;	// 12 useful for Davi EM
		absdev		= 0;	// 12 useful for Davi EM
		ismb		= false;
		isab		= false;
		abdbg		= false;

		inv_abscl	= 1.0/abscl;
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

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &za, "-za=%d", argv[i] ) )
			;
		else if( GetArg( &zb, "-zb=%d", argv[i] ) )
			;
		else if( GetArgStr( pat, "-p=", argv[i] ) )
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
		else if( IsArg( "-mb", argv[i] ) )
			ismb = true;
		else if( IsArg( "-ab", argv[i] ) )
			isab = true;
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

		if( U.T.t[2] + gW > w1 && U.T.t[2] < w2 &&
			U.T.t[5] + gH > h1 && U.T.t[5] < h2 ) {

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

		if( U.T.t[2] + gW > w1 && U.T.t[2] < w2 &&
			U.T.t[5] + gH > h1 && U.T.t[5] < h2 ) {

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

	double	sd = Normalize( v );

	if( !sd || !isfinite( sd ) ) {

		fprintf( flog,
		"FAIL: Image intersection region has stdev: %f\n", sd );

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
static void NewAngProc( double deg )
{
	double	r  = deg*PI/180,
			c  = cos( r ),
			s  = sin( r );
	int		Ox = int(gB->ws - (c*gA->ws - s*gA->hs))/2,
			Oy = int(gB->hs - (s*gA->ws + c*gA->hs))/2,
			Rx = int(gA->ws * 0.25),
			Ry = int(gB->hs * 0.25);

	gS->SetDisc( Ox, Oy, Rx, Ry );
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
// B-scape -> B-oriented
//	t.AddXY( B.x0, B.y0 );
//
// B-oriented -> B-montage
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
		t0 = StopTiming( flog, "MakeMontage", t0 );
	}

	if( !gArgs.isab )
		return;

	A.FindLayerIndices( gArgs.za );
	A.OrientLayer();

	MakeStripRasters( A, B );
	t0 = StopTiming( flog, "MakeStrips", t0 );

	A.MakePoints( thm.av, thm.ap );
	A.WriteMeta( 'A', gArgs.za, gArgs.abscl );

	B.MakePoints( thm.bv, thm.bp );
	B.WriteMeta( 'B', gArgs.zb, gArgs.abscl );

	thm.ftc.clear();
	thm.reqArea	= int(gW * gH * gArgs.inv_abscl * gArgs.inv_abscl);
	thm.olap1D	= int(gW * gArgs.inv_abscl * 0.20);
	thm.scl		= 1;

	S.Initialize( flog, best );
	S.SetRThresh( gArgs.abcorr );
	S.SetNbMaxHt( 0.99 );
	S.SetSweepConstXY( false );
	S.SetSweepPretweak( false );
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

		S.DenovoBestAngle( best, 0, 175, 5, thm );
		best.T.SetXY( best.X, best.Y );

		fprintf( flog, "*T: [0,1,2,3,4,5] (strip-strip)\n" );
		fprintf( flog, "[%g,%g,%g,%g,%g,%g]\n",
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


