

#include	"lsq_Globals.h"
#include	"lsq_Split.h"
#include	"lsq_MPI.h"

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
/* Run ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::Run()
{
	printf( "\n---- Split ----\n" );

	clock_t	t0 = StartTiming();

	Resize();

	map<int,int>	m;
	CountColors( m );

	StopTiming( stdout, "Split", t0 );
}


