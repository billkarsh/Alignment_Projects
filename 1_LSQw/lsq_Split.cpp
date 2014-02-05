

#include	"lsq_Globals.h"
#include	"lsq_Split.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

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

void Split::Resize()
{
	int	nz = zohi - zolo + 1;

	C.resize( nz );

	for( int iz = 0; iz < nz; ++iz )
		C[iz].resize( vR[iz].nr, 0 );
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

//		vector<int>&	c = C[iz];
		vector<uint8>&	c = vR[iz].flag;
		int				nr = vR[iz].nr;

		for( int ir = 0; ir < nr; ++ir ) {

			int	clr = c[ir];

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
		vector<int>&	c = ME->C[iz];
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

				if( c[j] == saveclr )
					c[j] = -1;	// mark it removed
				else if( c[j] > 0 )
					f[j] = fmRead;
			}
		}
		else {	// save all positive colors together

			for( int j = 0; j < R.nr; ++j ) {

				if( c[j] < 0 )	// previously removed
					f[j] = fmRead;
			}
		}

		SaveXBin( ME->X.X[iz], R.z );
		SaveFBin( R.flag, R.z );
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

	if( !EZThreads( _Save, nthr, 1, "Save" ) )
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
// removing it from the C vectors.

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

	map<int,int>	m;
	CountColors( m );
	BreakOut( m );

	StopTiming( stdout, "Split", t0 );
}


