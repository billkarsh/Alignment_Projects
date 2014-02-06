

#include	"lsq_Globals.h"
#include	"lsq_Split.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"

#include	<stack>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// colors per layer
#define	PERZ	10000

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static Split*		ME;
static const char	*gpath;
static int			nthr, saveclr;






/* --------------------------------------------------------------- */
/* Resize -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Allocate space for K vector.
//
void Split::Resize()
{
	int	nz = zohi - zolo + 1;

	K.resize( nz );

	for( int iz = 0; iz < nz; ++iz )
		K[iz].resize( vR[iz].nr, 0 );
}

/* --------------------------------------------------------------- */
/* _ColorMontages ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void* _ColorMontages( void* ithr )
{
	int	nz = zihi - zilo + 1;

	for( int iz = zilo + (long)ithr; iz < nz; iz += nthr ) {

		const Rgns&		R = vR[iz];
		vector<int>&	k = ME->K[iz];
		int				clr = (R.z ? PERZ * R.z : 1) - 1;

		for( int ir = 0; ir < R.nr; ++ir ) {

			// Valid tile?
			if( !FLAG_ISUSED( R.flag[ir] ) )
				continue;

			// Already assigned?
			if( k[ir] )
				continue;

			// New seed; new stack; new color
			stack<int>	s;
			s.push( ir );
			++clr;

			while( !s.empty() ) {

				// Process this guy

				int	jr = s.top();
				s.pop();

				if( k[jr] )
					continue;

				k[jr] = clr;

				// Push neighbors that are:
				// same-layer, virgin, valid.

				const vector<int>&	P  = R.pts[jr];
				int					np = P.size(), prev = -1;

				for( int ip = 0; ip < np; ++ip ) {

					const CorrPnt&	C = vC[P[ip]];

					if( !C.used )
						continue;

					if( C.z1 != iz || C.z2 != iz )
						continue;

					int other = (C.i1 == jr ? C.i2 : C.i1);

					if( other == prev )
						continue;

					prev = other;

					if( !k[other] && FLAG_ISUSED( R.flag[other] ) )
						s.push( other );
				}
			}
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* ColorMontages ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Walk pts to colorize islands per montage.
//
void Split::ColorMontages()
{
	int	nz = zihi - zilo + 1;

	ME		= (Split*)this;
	nthr	= maxthreads;

	if( nthr > nz )
		nthr = nz;

	if( !EZThreads( _ColorMontages, nthr, 1, "_ColorMontages" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* Send ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::Send( int zlo, int zhi, int toLorR )
{
	void	*buf;
	int		bytes,
			wdst;

	if( toLorR == 'L' ) {

		wdst = wkid - 1;

		if( wdst < 0 )
			return true;
	}
	else {
		wdst = wkid + 1;

		if( wdst >= nwks )
			return true;
	}

	for( int iz = zlo; iz <= zhi; ++iz ) {

		buf		= (void*)&K[iz][0];
		bytes	= sizeof(int) * vR[iz].nr;

		MPISend( buf, bytes, wdst, iz - zlo );
	}

	return true;
}

/* --------------------------------------------------------------- */
/* Recv ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::Recv( int zlo, int zhi, int fmLorR )
{
	void	*buf;
	int		bytes,
			wsrc;

	if( fmLorR == 'L' ) {

		wsrc = wkid - 1;

		if( wsrc < 0 )
			return true;
	}
	else {
		wsrc = wkid + 1;

		if( wsrc >= nwks )
			return true;
	}

	for( int iz = zlo; iz <= zhi; ++iz ) {

		buf		= (void*)&K[iz][0];
		bytes	= sizeof(int) * vR[iz].nr;

		MPIRecv( buf, bytes, wsrc, iz - zlo );
	}

	return true;
}

/* --------------------------------------------------------------- */
/* Updt ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::Updt()
{
	if( nwks <= 1 )
		return true;

// To avoid deadlocks, we arrange that at any given time,
// any node is just sending or just receiving. Easily done
// by odd/even role assignment.

// Odd send

	if( wkid & 1 ) {

		Send( zLlo, zLhi, 'L' );
		Send( zRlo, zRhi, 'R' );
	}
	else {

		Recv( zihi + 1, zohi, 'R' );
		Recv( zolo, zilo - 1, 'L' );
	}

// Even send

	if( !(wkid & 1) ) {

		Send( zLlo, zLhi, 'L' );
		Send( zRlo, zRhi, 'R' );
	}
	else {

		Recv( zihi + 1, zohi, 'R' );
		Recv( zolo, zilo - 1, 'L' );
	}

	return true;
}

/* --------------------------------------------------------------- */
/* ReportCount --------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::ReportCount( const map<int,int>& m )
{
	map<int,int>::const_iterator	it, ie = m.end();

	for( it = m.begin(); it != ie; ++it )
		printf( "Color %9d Count %9d\n", it->first, it->second );
}

/* --------------------------------------------------------------- */
/* CountColors --------------------------------------------------- */
/* --------------------------------------------------------------- */

class mpiclr {
// share color map with workers
public:
	int	first, second;
public:
	mpiclr() {};
	mpiclr( map<int,int>::iterator& it )
	: first(it->first), second(it->second) {};
};

void Split::CountColors( map<int,int>& m )
{
// Cache iterator and end()

	map<int,int>::iterator	it, ie = m.end();
	int						clast = -1;

// Calc my colors

	for( int iz = zilo; iz <= zihi; ++iz ) {

		vector<int>&	k = K[iz];
		int				nr = vR[iz].nr;

		for( int ir = 0; ir < nr; ++ir ) {

			int	clr = k[ir];

			if( !clr )
				continue;

			if( clr != clast )
				it = m.find( clast = clr );

			if( it != ie )
				++it->second;
			else {	// new entry: reset cache
				m[clr]	= 1;
				ie		= m.end();
				clast	= -1;
			}
		}
	}

// Report my subcounts

	if( nwks > 1 )
		printf( "This worker:\n" );
	else
		printf( "All workers:\n" );

	ReportCount( m );

// Send my colors to master; then get global colors

	if( wkid > 0 ) {

		// send count
		int	nm = m.size();
		MPISend( &nm, sizeof(int), 0, wkid );

		// send elems
		for( it = m.begin(); it != ie; ++it ) {
			mpiclr	c( it );
			MPISend( &c, sizeof(mpiclr), 0, wkid );
		}

		// get count
		m.clear();
		MPIRecv( &nm, sizeof(int), 0, wkid );

		// get elems
		for( int i = 0; i < nm; ++i ) {
			mpiclr	c;
			MPIRecv( &c, sizeof(mpiclr), 0, wkid );
			m[c.first] = c.second;
		}
	}
	else if( nwks > 1 ) {

		int	nm;

		// get each worker
		for( int iw = 1; iw < nwks; ++iw ) {

			MPIRecv( &nm, sizeof(int), iw, iw );

			for( int i = 0; i < nm; ++i ) {

				mpiclr	c;
				MPIRecv( &c, sizeof(mpiclr), iw, iw );

				it = m.find( c.first );

				if( it != m.end() )
					it->second += c.second;
				else
					m[c.first] = c.second;
			}
		}

		// send each worker
		nm = m.size();
		ie = m.end();

		for( int iw = 1; iw < nwks; ++iw ) {

			MPISend( &nm, sizeof(int), iw, iw );

			for( it = m.begin(); it != ie; ++it ) {
				mpiclr	c( it );
				MPISend( &c, sizeof(mpiclr), iw, iw );
			}
		}

		printf( "\nAll workers:\n" );
		ReportCount( m );
	}
}

/* --------------------------------------------------------------- */
/* _Save --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SaveXBin( const vector<double> &x, int z )
{
	char	buf[64];
	FILE	*f;

	sprintf( buf, "%s/X_%c_%d.bin",
		gpath, (ME->X.NE == 6 ? 'A' : 'H'), z );

	f = FileOpenOrDie( buf, "wb" );
	fwrite( &x[0], sizeof(double), x.size(), f );
	fclose( f );
}


static void SaveFBin( const vector<uint8> &f, int z )
{
	char	buf[64];
	FILE	*q;

	sprintf( buf, "%s/F_%d.bin", gpath, z );
	q = FileOpenOrDie( buf, "wb" );
	fwrite( &f[0], sizeof(uint8), f.size(), q );
	fclose( q );
}


static void* _Save( void* ithr )
{
	int	nz = zihi - zilo + 1;

	for( int iz = zilo + (long)ithr; iz < nz; iz += nthr ) {

		const Rgns&		R = vR[iz];
		vector<int>&	k = ME->K[iz];
		vector<uint8>	f = R.flag;

		// Flag adjustment:
		// Basically, we preserve existing flags,
		// but if it's not the right color, mark
		// it as 'transform missing'.
		//
		// Saveclr is set by caller:
		// >0 => match specific color
		// -1 => consolidate remaining colors
		//  0 => not implemented.

		if( saveclr >= 0 ) {

			for( int j = 0; j < R.nr; ++j ) {

				if( k[j] == saveclr )
					k[j] = -1;	// mark it removed
				else if( k[j] )	// kill all others
					f[j] = fmRead;
			}
		}
		else {	// save all positive colors together

			for( int j = 0; j < R.nr; ++j ) {

				if( k[j] < 0 )	// previously removed
					f[j] = fmRead;
			}
		}

		SaveXBin( ME->X.X[iz], R.z );
		SaveFBin( f, R.z );
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* Save ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Save tforms and <<modified>> flags...
//
void Split::Save()
{
	char	buf[64];
	int		nz = zihi - zilo + 1;

	ME		= (Split*)this;
	gpath	= buf;
	nthr	= maxthreads;

	if( nthr > nz )
		nthr = nz;

	if( saveclr > 0 ) {
		sprintf( buf, "Splits/X_%c_BIN_CLR%d",
		(X.NE == 6 ? 'A' : 'H'), saveclr );
	}
	else {
		sprintf( buf, "Splits/X_%c_BIN_REM",
		(X.NE == 6 ? 'A' : 'H') );
	}

	DskCreateDir( buf, stdout );

	if( !EZThreads( _Save, nthr, 1, "_SplitSave" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* BreakOut ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void Split::BreakOut( const map<int,int>& m )
{
	int	nm = m.size();

	if( nm <= 1 ) {
		printf(
		"No splits created;"
		" (single color or all below splitmin).\n" );
		return;
	}

// Splits folder

	DskCreateDir( "Splits", stdout );

// First separately write each color over threshold.
// The Save operation will replace that color by (-1)
// removing it from the K vectors.

	map<int,int>::const_iterator	it, ie = m.end();

	for( it = m.begin(); it != ie; ++it ) {

		if( it->second >= splitmin ) {

			saveclr = it->first;
			Save();
			--nm;
		}
	}

// Now consolidate remaining 'REM' fragments.

	if( nm > 0 ) {
		saveclr = -1;
		Save();
	}
}

/* --------------------------------------------------------------- */
/* Run ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::Run()
{
	printf( "\n---- Split ----\n" );

	clock_t	t0 = StartTiming();

	Resize();
	ColorMontages();
	Updt();

	map<int,int>	m;
	CountColors( m );
//	BreakOut( m );

	StopTiming( stdout, "Split", t0 );
}


