

#include	"lsq_Globals.h"
#include	"lsq_Magnitude.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"TAffine.h"
#include	"THmgphy.h"
#include	"Timer.h"

#include	<string.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#define TOPN	10
#define	TOPL	(TOPN-1)

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

namespace magnitude {

class EI {
// Error and point using local indexing
public:
	double	e;
	int		iz, i;	// local tform indices
public:
	EI() : e(-1), iz(0), i(-1) {};
};

static bool Sort_EI_inc( const EI &A, const EI &B )
{
	if( A.e == -1 )
		return false;
	else if( B.e == -1 )
		return true;
	return A.e < B.e;
}

static bool Sort_EI_dec( const EI &A, const EI &B )
{
	return A.e > B.e;
}

class EG {
// Error and point using global indexing
public:
	double	e;
	int		z, i, r;	// global indices
public:
	EG() : e(-1) {};

	void FromEI( const EI &rhs );
};

static bool Sort_EG_inc( const EG &A, const EG &B )
{
	if( A.e == -1 )
		return false;
	else if( B.e == -1 )
		return true;
	return A.e < B.e;
}

static bool Sort_EG_dec( const EG &A, const EG &B )
{
	return A.e > B.e;
}

class Stat {
// statistics summaries using local point indexing
private:
	vector<EI>::iterator	eis0, eib0;
public:
	vector<EI>	eis, eib;	// topn
	double		sum;		// sum
	int			n;			// count
	EI			cur;		// current
public:
	// accumulate layerwise data
	void Init( int iz );
	void Add();

	// combine
	void Init_Smy( const Stat &rhs );
	void Add_Smy( const Stat &rhs );
	void SmyMyLayers();

	// report
	double RMS() const	{return (n ? sqrt( sum/(2*n) ) : -1);};
	void Topn( FILE *f, int SorB ) const;
};

class StatG {
// Stat type using cross-worker global indexing
private:
	typedef struct {
		double	sum;
		int		n;
		EG		egs[TOPN],
				egb[TOPN];
	} MPIBUF;
private:
	vector<EG>::iterator	egs0, egb0;
public:
	vector<EG>	egs, egb;	// topn
	double		sum;		// sum
	int			n;			// count
public:
	// Scatter, gather
	void FromStat( const Stat &rhs );
	void Send();
	void Recv( int iw );
	void Add( const StatG &rhs );

	// report
	double RMS() const	{return (n ? sqrt( sum/(2*n) ) : -1);};
	void Topn( FILE *f, int SorB ) const;
};

}	// namespace magnitude

using namespace magnitude;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const XArray	*gX;
static vector<Stat>	vS;
static int			nthr;






/* --------------------------------------------------------------- */
/* EG::FromEI ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void EG::FromEI( const EI &rhs )
{
	if( (e = rhs.e) > -1 )
		RealZIDR( z, i, r, rhs.iz, rhs.i );
}

/* --------------------------------------------------------------- */
/* Stat::Init ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Init( int iz )
{
	eis.resize( TOPN );
	eib.resize( TOPN );
	eis0	= eis.begin();
	eib0	= eib.begin();
	sum		= 0.0;
	n		= 0;
	cur.iz	= iz;
}

/* --------------------------------------------------------------- */
/* Stat::Add ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Add()
{
	sum += cur.e;
	++n;

// smallest
	{
		EI&	ei = eis[TOPL];

		if( ei.e == -1 || cur.e < ei.e ) {
			ei = cur;
			sort( eis0, eis0 + TOPN, Sort_EI_inc );
		}
	}

// biggest
	{
		EI&	ei = eib[TOPL];

		if( cur.e > ei.e ) {
			ei = cur;
			sort( eib0, eib0 + TOPN, Sort_EI_dec );
		}
	}
}

/* --------------------------------------------------------------- */
/* Stat::Init_Smy ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Use one Stat item to initialize cummulative summary.
//
void Stat::Init_Smy( const Stat &rhs )
{
	eis.resize( 2*TOPN );
	eib.resize( 2*TOPN );
	eis0 = eis.begin();
	eib0 = eib.begin();
	memcpy( &eis[0], &rhs.eis[0], TOPN * sizeof(EI) );
	memcpy( &eib[0], &rhs.eib[0], TOPN * sizeof(EI) );
	sum	= rhs.sum;
	n	= rhs.n;
}

/* --------------------------------------------------------------- */
/* Stat::Add_Smy ------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Add_Smy( const Stat &rhs )
{
	memcpy( &eis[TOPN], &rhs.eis[0], TOPN * sizeof(EI) );
	sort( eis0, eis0 + 2*TOPN, Sort_EI_inc );
	memcpy( &eib[TOPN], &rhs.eib[0], TOPN * sizeof(EI) );
	sort( eib0, eib0 + 2*TOPN, Sort_EI_dec );
	sum += rhs.sum;
	n	+= rhs.n;
}

/* --------------------------------------------------------------- */
/* Stat::SmyMyLayers --------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::SmyMyLayers()
{
	int	nL = zihi - zilo + 1;

	Init_Smy( vS[0] );

	for( int iL = 1; iL < nL; ++iL )
		Add_Smy( vS[iL] );
}

/* --------------------------------------------------------------- */
/* Stat::Topn ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Topn( FILE *f, int SorB ) const
{
// print the errors...

	const vector<EI>&	vei = (SorB == 'S' ? eis : eib);
	int					ne  = 0;

	for( int i = 0; i < TOPN; ++i ) {

		double	e = vei[i].e;

		if( e > -1 ) {
			fprintf( f, "\t%.3f", sqrt( e/2 ) );
			++ne;
		}
	}

// ...and their labels

	if( ne ) {

		for( int ie = 0; ie < ne; ++ie ) {

			int	z, i, r;
			RealZIDR( z, i, r, vei[ie].iz, vei[ie].i );
			fprintf( f, "\t%d.%d-%d", z, i, r );
		}
	}
	else
		fprintf( f, "-1\n" );

	fprintf( f, "\n" );
}

/* --------------------------------------------------------------- */
/* StatG::FromStat ----------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::FromStat( const Stat &rhs )
{
	egs.resize( 2*TOPN );
	egb.resize( 2*TOPN );
	egs0	= egs.begin();
	egb0	= egb.begin();
	sum		= rhs.sum;
	n		= rhs.n;

	for( int i = 0; i < TOPN; ++i ) {
		egs[i].FromEI( rhs.eis[i] );
		egb[i].FromEI( rhs.eib[i] );
	}
}

/* --------------------------------------------------------------- */
/* StatG::Send --------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Send()
{
	MPIBUF	B;

	B.sum = sum;
	B.n   = n;
	memcpy( B.egs, &egs[0], TOPN * sizeof(EG) );
	memcpy( B.egb, &egb[0], TOPN * sizeof(EG) );

	MPISend( &B, sizeof(MPIBUF), 0, wkid );
}

/* --------------------------------------------------------------- */
/* StatG::Recv ----------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Recv( int iw )
{
	MPIBUF	B;
	MPIRecv( &B, sizeof(MPIBUF), iw, iw );

	egs.resize( 2*TOPN );
	egs0 = egs.begin();
	egb.resize( 2*TOPN );
	egb0 = egb.begin();

	sum = B.sum;
	n   = B.n;
	memcpy( &egs[0], B.egs, TOPN * sizeof(EG) );
	memcpy( &egb[0], B.egb, TOPN * sizeof(EG) );
}

/* --------------------------------------------------------------- */
/* StatG::Add ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Add( const StatG &rhs )
{
	memcpy( &egs[TOPN], &rhs.egs[0], TOPN * sizeof(EG) );
	sort( egs0, egs0 + 2*TOPN, Sort_EG_inc );
	memcpy( &egb[TOPN], &rhs.egb[0], TOPN * sizeof(EG) );
	sort( egb0, egb0 + 2*TOPN, Sort_EG_dec );
	sum += rhs.sum;
	n	+= rhs.n;
}

/* --------------------------------------------------------------- */
/* StatG::Topn --------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Topn( FILE *f, int SorB ) const
{
// print the errors...

	const vector<EG>&	veg = (SorB == 'S' ? egs : egb);
	int					ne  = 0;

	for( int i = 0; i < TOPN; ++i ) {

		double	e = veg[i].e;

		if( e > -1 ) {
			fprintf( f, "\t%.3f", sqrt( e/2 ) );
			++ne;
		}
	}

// ...and their labels

	if( ne ) {

		for( int i = 0; i < ne; ++i ) {

			const EG&	G = veg[i];

			fprintf( f, "\t%d.%d-%d", G.z, G.i, G.r );
		}
	}
	else
		fprintf( f, "-1\n" );

	fprintf( f, "\n" );
}

/* --------------------------------------------------------------- */
/* _Mag ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Mag( void* ithr )
{
	int	ns = zihi - zilo + 1;

// For each layer...

	for( int is = (long)ithr; is < ns; is += nthr ) {

		int						iz	= is + zilo;
		Stat&					S	= vS[is];
		const Rgns&				R	= vR[iz];
		const vector<double>&	x	= gX->X[iz];

		S.Init( iz );

		// For each rgn...

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( !FLAG_ISUSED( R.flag[ir] ) )
				continue;

			Point	p0, p1( 1, 1 );

			if( gX->NE == 6 ) {
				const TAffine&	T = X_AS_AFF( x, ir );
				T.Transform( p0 );
				T.Transform( p1 );
			}
			else {
				const THmgphy&	T = X_AS_HMY( x, ir );
				T.Transform( p0 );
				T.Transform( p1 );
			}

			S.cur.i = ir;
			S.cur.e = p1.DistSqr( p0 );
			S.Add();
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* CalcLayerwiseMags --------------------------------------------- */
/* --------------------------------------------------------------- */

static void CalcLayerwiseMags( const XArray &X )
{
	gX = &X;

	int	ns = zihi - zilo + 1;

	vS.resize( ns );

	nthr = maxthreads;

	if( nthr > ns )
		nthr = ns;

	if( !EZThreads( _Mag, nthr, 1, "_Mag" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteLocalFiles ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteLocalFiles()
{
	char	buf[64];
	FILE	*f;

// Smalls

	if( !wkid )
		strcpy( buf, "MagSmall.txt" );
	else
		sprintf( buf, "MagTemp/topn_S_%d.txt", wkid );

	f = FileOpenOrDie( buf, "w" );

	// header

	if( !wkid )
		fprintf( f, "Z\tRMS\tTOPN\n" );

	for( int iz = zilo; iz <= zihi; ++iz ) {

		const Stat&	S = vS[iz - zilo];

		fprintf( f, "%d\t%.3f", vR[iz].z, S.RMS() );
		S.Topn( f, 'S' );
	}

	fclose( f );

// Bigs

	if( !wkid )
		strcpy( buf, "MagBig.txt" );
	else
		sprintf( buf, "MagTemp/topn_B_%d.txt", wkid );

	f = FileOpenOrDie( buf, "w" );

	// header

	if( !wkid )
		fprintf( f, "Z\tRMS\tTOPN\n" );

	for( int iz = zilo; iz <= zihi; ++iz ) {

		const Stat&	S = vS[iz - zilo];

		fprintf( f, "%d\t%.3f", vR[iz].z, S.RMS() );
		S.Topn( f, 'B' );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* LogLocalSmy --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LogLocalSmy( Stat &Sw, const char* title )
{
	printf( "%s:\n", title );

	printf( "Small RMS %.3f TopN", Sw.RMS() );
	Sw.Topn( stdout, 'S' );

	printf( "\nBig RMS %.3f TopN", Sw.RMS() );
	Sw.Topn( stdout, 'B' );
}

/* --------------------------------------------------------------- */
/* LogGlobalSmy -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LogGlobalSmy( StatG &S0 )
{
	printf( "\nAll workers:\n" );

	printf( "Small RMS %.3f TopN", S0.RMS() );
	S0.Topn( stdout, 'S' );

	printf( "\nBig RMS %.3f TopN", S0.RMS() );
	S0.Topn( stdout, 'B' );
}

/* --------------------------------------------------------------- */
/* Consolidate --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Consolidate()
{
	Stat	Sw;

	Sw.SmyMyLayers();

	if( nwks <= 1 )
		LogLocalSmy( Sw, "All workers" );
	else if( wkid > 0 ) {

		StatG	Sg;
		Sg.FromStat( Sw );
		Sg.Send();

		LogLocalSmy( Sw, "This worker" );
	}
	else {

		StatG	S0;
		S0.FromStat( Sw );

		for( int iw = 1; iw < nwks; ++iw ) {

			// accumulate worker stats

			StatG	Si;
			Si.Recv( iw );
			S0.Add( Si );

			// append temp files

			char	buf[128];

			sprintf( buf,
			"cat MagTemp/topn_S_%d.txt >> MagSmall.txt", iw );
			system( buf );

			sprintf( buf,
			"cat MagTemp/topn_B_%d.txt >> MagBig.txt", iw );
			system( buf );
		}

		system( "rm -rf MagTemp" );
		LogLocalSmy( Sw, "This worker" );
		LogGlobalSmy( S0 );
	}
}

/* --------------------------------------------------------------- */
/* Magnitude ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Reports on tforms size distortions:
//
// - Files 'MagBig.txt' and 'MagSmall.txt' tabulate, for
//		each layer, the RMS and the topn largest distortions.
//
// - Logs summarize RMS and topn over {bigs, smalls}, shown
//		by each worker and over all workers.
//
void Magnitude( const XArray &X )
{
	printf( "\n---- Magnitudes ----\n" );

	clock_t	t0 = StartTiming();

	if( nwks > 1 )
		DskCreateDir( "MagTemp", stdout );

	CalcLayerwiseMags( X );
	WriteLocalFiles();
	Consolidate();
	vS.clear();

	StopTiming( stdout, "Mags", t0 );
}


