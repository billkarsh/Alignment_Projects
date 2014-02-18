

#include	"lsq_Error.h"
#include	"lsq_Globals.h"
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

namespace error {

class FileErr {
// Buffered file writing
private:
	static const int bufsize = 2048;
	FILE*			f;
	vector<float>	ve;
	int				n;
public:
	FileErr( int SorD, int z );
	virtual ~FileErr();
	void Add( double err );
};

class EI {
// Error and point using local indexing
public:
	double	e;
	int		i;	// index into vC[]
public:
	EI() : e(-1), i(-1) {};

	bool operator < ( const EI &rhs ) const
		{return e > rhs.e;};
};

class EG {
// Error and point using global indexing
public:
	double	e;
	int		z1, i1, r1,	// global idb values
			z2, i2, r2;
public:
	EG() : e(-1) {};

	void FromEI( const EI &rhs );

	bool operator < ( const EG &rhs ) const
		{return e > rhs.e;};
};

class Stat {
// statistics summaries using local point indexing
private:
	vector<EI>::iterator	eis0, eid0;
public:
	vector<EI>	eis, eid;	// topn
	double		sms, smd;	// sum
	int			ns,  nd;	// count
	EI			cur;		// current
public:
	// accumulate layerwise data
	void Init();
	void AddS();
	void AddD();

	// combine
	void Init_Smy( const Stat &rhs );
	void Add_Smy( const Stat &rhs );
	void SmyMyLayers();
	void Total( const Stat &rhs );

	// report
	double RMSS() const	{return (ns ? sqrt( sms/ns ) : -1);};
	double RMSD() const	{return (nd ? sqrt( smd/nd ) : -1);};
	void Topn( FILE *f, int SorD ) const;
};

class StatG {
// Stat type using cross-worker global indexing
private:
	typedef struct {
		double	sms, smd;
		int		ns,  nd;
		EG		egs[TOPN],
				egd[TOPN];
	} MPIBUF;
private:
	vector<EG>::iterator	egs0, egd0;
public:
	vector<EG>	egs, egd;	// topn
	double		sms, smd;	// sum
	int			ns,  nd;	// count
public:
	// Scatter, gather
	void FromStat( const Stat &rhs );
	void Send();
	void Recv( int iw );
	void Add( const StatG &rhs );

	// combine
	void Total( const StatG &rhs );

	// report
	double RMSS() const	{return (ns ? sqrt( sms/ns ) : -1);};
	double RMSD() const	{return (nd ? sqrt( smd/nd ) : -1);};
	void Topn( FILE *f, int SorD ) const;
};

}	// namespace error

using namespace error;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const XArray	*gX;
static vector<Stat>	vS;
static double		fnlrms = -1,
					fnlmax = -1;
static int			nthr;






/* --------------------------------------------------------------- */
/* FileErr::FileErr ---------------------------------------------- */
/* --------------------------------------------------------------- */

FileErr::FileErr( int SorD, int z ) : f(NULL), n(0)
{
	if( SorD == 'S' || zolo != zohi ) {

		char	buf[32];
		sprintf( buf, "Error/Err_%c_%d.bin", SorD, z );
		f = FileOpenOrDie( buf, "wb" );

		ve.resize( bufsize );
	}
}

/* --------------------------------------------------------------- */
/* FileErr::~FileErr --------------------------------------------- */
/* --------------------------------------------------------------- */

FileErr::~FileErr()
{
	if( f ) {

		if( n )
			fwrite( &ve[0], sizeof(float), n, f );

		fclose( f );
	}
}

/* --------------------------------------------------------------- */
/* FileErr::Add -------------------------------------------------- */
/* --------------------------------------------------------------- */

void FileErr::Add( double err )
{
	ve[n++] = (float)sqrt( err );

	if( n >= bufsize ) {
		fwrite( &ve[0], sizeof(float), n, f );
		n = 0;
	}
}

/* --------------------------------------------------------------- */
/* EG::FromEI ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void EG::FromEI( const EI &rhs )
{
	const CorrPnt&	C = vC[rhs.i];

	if( (e = rhs.e) > -1 ) {
		RealZIDR( z1, i1, r1, C.z1, C.i1 );
		RealZIDR( z2, i2, r2, C.z2, C.i2 );
	}
}

/* --------------------------------------------------------------- */
/* Stat::Init ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Init()
{
	eis.resize( TOPN );
	eis0	= eis.begin();
	sms		= 0.0;
	ns		= 0;

	if( zolo != zohi ) {
		eid.resize( TOPN );
		eid0	= eid.begin();
		smd		= 0.0;
		nd		= 0;
	}
}

/* --------------------------------------------------------------- */
/* Stat::AddS ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::AddS()
{
	EI&		ei = eis[TOPL];

	sms += cur.e;
	++ns;

	if( cur.e > ei.e ) {
		ei = cur;
		sort( eis0, eis0 + TOPN );
	}
}

/* --------------------------------------------------------------- */
/* Stat::AddD ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::AddD()
{
	EI&		ei = eid[TOPL];

	smd += cur.e;
	++nd;

	if( cur.e > ei.e ) {
		ei = cur;
		sort( eid0, eid0 + TOPN );
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
	eis0 = eis.begin();
	memcpy( &eis[0], &rhs.eis[0], TOPN * sizeof(EI) );
	sms	= rhs.sms;
	ns	= rhs.ns;

	if( zolo != zohi ) {
		eid.resize( 2*TOPN );
		eid0 = eid.begin();
		memcpy( &eid[0], &rhs.eid[0], TOPN * sizeof(EI) );
		smd	= rhs.smd;
		nd	= rhs.nd;
	}
}

/* --------------------------------------------------------------- */
/* Stat::Add_Smy ------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Add_Smy( const Stat &rhs )
{
	sms += rhs.sms;
	ns	+= rhs.ns;

	memcpy( &eis[TOPN], &rhs.eis[0], TOPN * sizeof(EI) );
	sort( eis0, eis0 + 2*TOPN );

	if( zolo != zohi ) {

		smd += rhs.smd;
		nd	+= rhs.nd;

		memcpy( &eid[TOPN], &rhs.eid[0], TOPN * sizeof(EI) );
		sort( eid0, eid0 + 2*TOPN );
	}
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
/* Stat::Total --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set this same as combination of rhs same and down data.
//
void Stat::Total( const Stat &rhs )
{
	eis.resize( 2*TOPN );
	eis0	= eis.begin();
	sms		= rhs.sms + rhs.smd;
	ns		= rhs.ns  + rhs.nd;

	memcpy( &eis[0],    &rhs.eis[0], TOPN * sizeof(EI) );
	memcpy( &eis[TOPN], &rhs.eid[0], TOPN * sizeof(EI) );
	sort( eis0, eis0 + 2*TOPN );
}

/* --------------------------------------------------------------- */
/* Stat::Topn ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Stat::Topn( FILE *f, int SorD ) const
{
// print the errors...

	const vector<EI>&	vei = (SorD == 'S' ? eis : eid);
	int					ne  = 0;

	for( int i = 0; i < TOPN; ++i ) {

		double	e = vei[i].e;

		if( e > -1 ) {
			fprintf( f, "\t%.2f", sqrt( e ) );
			++ne;
		}
	}

// ...and their labels

	if( ne ) {

		for( int i = 0; i < ne; ++i ) {

			const CorrPnt&	C = vC[vei[i].i];
			int				z1, i1, r1,
							z2, i2, r2;

			RealZIDR( z1, i1, r1, C.z1, C.i1 );
			RealZIDR( z2, i2, r2, C.z2, C.i2 );

			fprintf( f, "\t%d.%d-%d^%d.%d-%d",
				z1, i1, r1, z2, i2, r2 );
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
	egs0	= egs.begin();
	sms		= rhs.sms;
	ns		= rhs.ns;

	for( int i = 0; i < TOPN; ++i )
		egs[i].FromEI( rhs.eis[i] );

	egd.resize( 2*TOPN );
	smd		= rhs.smd;
	nd		= rhs.nd;

	if( zolo != zohi ) {

		egd0 = egd.begin();

		for( int i = 0; i < TOPN; ++i )
			egd[i].FromEI( rhs.eid[i] );
	}
}

/* --------------------------------------------------------------- */
/* StatG::Send --------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Send()
{
	MPIBUF	B;

	B.sms = sms;
	B.smd = smd;
	B.ns  = ns;
	B.nd  = nd;
	memcpy( B.egs, &egs[0], TOPN * sizeof(EG) );
	memcpy( B.egd, &egd[0], TOPN * sizeof(EG) );

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

	egd.resize( 2*TOPN );
	egd0 = egd.begin();

	sms = B.sms;
	smd = B.smd;
	ns  = B.ns;
	nd  = B.nd;
	memcpy( &egs[0], B.egs, TOPN * sizeof(EG) );
	memcpy( &egd[0], B.egd, TOPN * sizeof(EG) );
}

/* --------------------------------------------------------------- */
/* StatG::Add ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Add( const StatG &rhs )
{
	sms += rhs.sms;
	ns	+= rhs.ns;

	memcpy( &egs[TOPN], &rhs.egs[0], TOPN * sizeof(EG) );
	sort( egs0, egs0 + 2*TOPN );

	if( zolo != zohi ) {

		smd += rhs.smd;
		nd	+= rhs.nd;

		memcpy( &egd[TOPN], &rhs.egd[0], TOPN * sizeof(EG) );
		sort( egd0, egd0 + 2*TOPN );
	}
}

/* --------------------------------------------------------------- */
/* StatG::Total -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set this same as combination of rhs same and down data.
//
void StatG::Total( const StatG &rhs )
{
	egs.resize( 2*TOPN );
	egs0	= egs.begin();
	sms		= rhs.sms + rhs.smd;
	ns		= rhs.ns  + rhs.nd;

	memcpy( &egs[0],    &rhs.egs[0], TOPN * sizeof(EG) );
	memcpy( &egs[TOPN], &rhs.egd[0], TOPN * sizeof(EG) );
	sort( egs0, egs0 + 2*TOPN );
}

/* --------------------------------------------------------------- */
/* StatG::Topn --------------------------------------------------- */
/* --------------------------------------------------------------- */

void StatG::Topn( FILE *f, int SorD ) const
{
// print the errors...

	const vector<EG>&	veg = (SorD == 'S' ? egs : egd);
	int					ne  = 0;

	for( int i = 0; i < TOPN; ++i ) {

		double	e = veg[i].e;

		if( e > -1 ) {
			fprintf( f, "\t%.2f", sqrt( e ) );
			++ne;
		}
	}

// ...and their labels

	if( ne ) {

		for( int i = 0; i < ne; ++i ) {

			const EG&	G = veg[i];

			fprintf( f, "\t%d.%d-%d^%d.%d-%d",
				G.z1, G.i1, G.r1, G.z2, G.i2, G.r2 );
		}
	}
	else
		fprintf( f, "-1\n" );

	fprintf( f, "\n" );
}

/* --------------------------------------------------------------- */
/* _ErrorA ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _ErrorA( void* ithr )
{
	int	ns = zihi - zilo + 1;

// For each layer...

	for( int is = (long)ithr; is < ns; is += nthr ) {

		int						iz	= is + zilo;
		Stat&					S	= vS[is];
		const Rgns&				Ra	= vR[iz];
		const vector<double>&	xa	= gX->X[iz];
		FileErr					FS( 'S', Ra.z ),
								FD( 'D', Ra.z );

		S.Init();

		// For each rgn...

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( !FLAG_ISUSED( Ra.flag[ir] ) )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const TAffine*		Ta = &X_AS_AFF( xa, ir );
			const TAffine*		Tb;
			int					lastbi,
								lastbz	= -1,
								np		= P.size();

			// For each of its points...

			for( int ip = 0; ip < np; ++ip ) {

				CorrPnt&	C = vC[S.cur.i = P[ip]];

				if( !C.used )
					continue;

				if( C.z1 == C.z2 ) {

					// no double counting
					if( C.i1 != ir )
						continue;

					if( C.z2 != lastbz ) {
						lastbz = C.z2;
						lastbi = -1;
					}

					if( C.i2 != lastbi ) {

						if( !FLAG_ISUSED( Ra.flag[C.i2] ) )
							continue;

						Tb = &X_AS_AFF( xa, C.i2 );
						lastbi = C.i2;
					}

					Point	pa = C.p1,
							pb = C.p2;

					Ta->Transform( pa );
					Tb->Transform( pb );

					S.cur.e = pb.DistSqr( pa );
					S.AddS();
					FS.Add( S.cur.e );
				}
				else if( C.z1 == zolo )
					continue;
				else if( FLAG_ISCUTD( Ra.flag[ir] ) )
					continue;
				else if( C.z1 == iz ) {

					if( C.z2 != lastbz ) {
						lastbz = C.z2;
						lastbi = -1;
					}

					if( C.i2 != lastbi ) {

						if( !FLAG_ISUNCT( vR[C.z2].flag[C.i2] ) )
							continue;

						Tb = &X_AS_AFF( gX->X[C.z2], C.i2 );
						lastbi = C.i2;
					}

					Point	pa = C.p1,
							pb = C.p2;

					Ta->Transform( pa );
					Tb->Transform( pb );

					S.cur.e = pb.DistSqr( pa );
					S.AddD();
					FD.Add( S.cur.e );
				}
			}
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* _ErrorH ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _ErrorH( void* ithr )
{
	int	ns = zihi - zilo + 1;

// For each layer...

	for( int is = (long)ithr; is < ns; is += nthr ) {

		int						iz	= is + zilo;
		Stat&					S	= vS[is];
		const Rgns&				Ra	= vR[iz];
		const vector<double>&	xa	= gX->X[iz];
		FileErr					FS( 'S', Ra.z ),
								FD( 'D', Ra.z );

		S.Init();

		// For each rgn...

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( !FLAG_ISUSED( Ra.flag[ir] ) )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const THmgphy*		Ta = &X_AS_HMY( xa, ir );
			const THmgphy*		Tb;
			int					lastbi,
								lastbz	= -1,
								np		= P.size();

			// For each of its points...

			for( int ip = 0; ip < np; ++ip ) {

				CorrPnt&	C = vC[S.cur.i = P[ip]];

				if( !C.used )
					continue;

				if( C.z1 == C.z2 ) {

					// no double counting
					if( C.i1 != ir )
						continue;

					if( C.z2 != lastbz ) {
						lastbz = C.z2;
						lastbi = -1;
					}

					if( C.i2 != lastbi ) {

						if( !FLAG_ISUSED( Ra.flag[C.i2] ) )
							continue;

						Tb = &X_AS_HMY( xa, C.i2 );
						lastbi = C.i2;
					}

					Point	pa = C.p1,
							pb = C.p2;

					Ta->Transform( pa );
					Tb->Transform( pb );

					S.cur.e = pb.DistSqr( pa );
					S.AddS();
					FS.Add( S.cur.e );
				}
				else if( C.z1 == zolo )
					continue;
				else if( FLAG_ISCUTD( Ra.flag[ir] ) )
					continue;
				else if( C.z1 == iz ) {

					if( C.z2 != lastbz ) {
						lastbz = C.z2;
						lastbi = -1;
					}

					if( C.i2 != lastbi ) {

						if( !FLAG_ISUNCT( vR[C.z2].flag[C.i2] ) )
							continue;

						Tb = &X_AS_HMY( gX->X[C.z2], C.i2 );
						lastbi = C.i2;
					}

					Point	pa = C.p1,
							pb = C.p2;

					Ta->Transform( pa );
					Tb->Transform( pb );

					S.cur.e = pb.DistSqr( pa );
					S.AddD();
					FD.Add( S.cur.e );
				}
			}
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* CalcLayerwiseError -------------------------------------------- */
/* --------------------------------------------------------------- */

static void CalcLayerwiseError( const XArray &X )
{
	gX = &X;

	int	ns = zihi - zilo + 1;

	vS.resize( ns );

	nthr = maxthreads;

	if( nthr > ns )
		nthr = ns;

// Select affine or hmgphy

	EZThreadproc	proc;
	const char		*sproc;

	if( X.NE == 6 ) {
		proc	= _ErrorA;
		sproc	= "_ErrorA";
	}
	else {
		proc	= _ErrorH;
		sproc	= "_ErrorH";
	}

	if( !EZThreads( proc, nthr, 1, sproc ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteLocalFiles ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteLocalFiles()
{
	char	buf[64];
	FILE	*f;

// Sames

	if( !wkid )
		strcpy( buf, "ErrSame.txt" );
	else
		sprintf( buf, "ErrTemp/topn_S_%d.txt", wkid );

	f = FileOpenOrDie( buf, "w" );

	// header

	if( !wkid )
		fprintf( f, "Z\tRMS\tTOPN\n" );

	for( int iz = zilo; iz <= zihi; ++iz ) {

		const Stat&	S = vS[iz - zilo];

		fprintf( f, "%d\t%.2f", vR[iz].z, S.RMSS() );
		S.Topn( f, 'S' );
	}

	fclose( f );

// Downs

	if( zolo == zohi )
		return;

	if( !wkid )
		strcpy( buf, "ErrDown.txt" );
	else
		sprintf( buf, "ErrTemp/topn_D_%d.txt", wkid );

	f = FileOpenOrDie( buf, "w" );

	// header

	if( !wkid )
		fprintf( f, "Z\tRMS\tTOPN\n" );

	for( int iz = zilo + (zilo == zolo); iz <= zihi; ++iz ) {

		const Stat&	S = vS[iz - zilo];

		fprintf( f, "%d\t%.2f", vR[iz].z, S.RMSD() );
		S.Topn( f, 'D' );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* LogLocalSmy --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LogLocalSmy( Stat &Sw, const char* title )
{
	printf( "%s:\n", title );

	printf( "Same RMS %.2f TopN", Sw.RMSS() );
	Sw.Topn( stdout, 'S' );

	fnlrms	= Sw.RMSS();
	fnlmax	= Sw.eis[0].e;

	if( zolo != zohi ) {

		printf( "\nDown RMS %.2f TopN", Sw.RMSD() );
		Sw.Topn( stdout, 'D' );

		Stat	St;
		St.Total( Sw );
		printf( "\nAll RMS %.2f TopN", St.RMSS() );
		St.Topn( stdout, 'S' );

		fnlrms	= St.RMSS();
		fnlmax	= St.eis[0].e;
	}
}

/* --------------------------------------------------------------- */
/* LogGlobalSmy -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LogGlobalSmy( StatG &S0 )
{
	printf( "\nAll workers:\n" );

	printf( "Same RMS %.2f TopN", S0.RMSS() );
	S0.Topn( stdout, 'S' );

	printf( "\nDown RMS %.2f TopN", S0.RMSD() );
	S0.Topn( stdout, 'D' );

	StatG	St;
	St.Total( S0 );
	printf( "\nAll RMS %.2f TopN", St.RMSS() );
	St.Topn( stdout, 'S' );

	fnlrms	= St.RMSS();
	fnlmax	= St.egs[0].e;
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
			"cat ErrTemp/topn_S_%d.txt >> ErrSame.txt", iw );
			system( buf );

			sprintf( buf,
			"cat ErrTemp/topn_D_%d.txt >> ErrDown.txt", iw );
			system( buf );
		}

		system( "rm -rf ErrTemp" );
		LogLocalSmy( Sw, "This worker" );
		LogGlobalSmy( S0 );
	}
}

/* --------------------------------------------------------------- */
/* Error --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Using the vC inliers (used = true) calculate several metrics
// from the residual errors:
//
// - Files 'ErrSame.txt' and 'ErrDown.txt' tabulate, for
//		each layer, the RMS and the topn largest errors.
//
// - Logs summarize RMS and topn over {just sames, downs,
//		all}, shown by each worker and over all workers.
//
// - Folder 'Error' with files-by-layer 'Err_S_i.bin' and
//		'Err_D_i.bin' with packed |err| values as floats.
//		These are histogrammed using separate eview tool.
//
void Error( const XArray &X )
{
	printf( "\n---- Error statistics ----\n" );

	clock_t	t0 = StartTiming();

	DskCreateDir( "Error", stdout );

	if( nwks > 1 )
		DskCreateDir( "ErrTemp", stdout );

	CalcLayerwiseError( X );
	WriteLocalFiles();
	Consolidate();
	vS.clear();

	StopTiming( stdout, "Errors", t0 );
}

/* --------------------------------------------------------------- */
/* GetFinalError ------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetFinalError( double &erms, double &emax )
{
	erms	= fnlrms;
	emax	= (fnlmax >= 0 ? sqrt( fnlmax ) : -1);
}


