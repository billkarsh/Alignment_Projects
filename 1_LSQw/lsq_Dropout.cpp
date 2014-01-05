

#include	"lsq_Dropout.h"
#include	"lsq_Globals.h"
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

static vector<Dropout>	vD;
static int				nthr;






/* --------------------------------------------------------------- */
/* _Scan --------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Scan( void* ithr )
{
	int	nd = zihi - zilo + 1;

// For each layer...

	for( int id = (long)ithr; id < nd; id += nthr ) {

		int			iz	= id + zilo;
		const Rgns&	R	= vR[iz];
		Dropout&	D	= vD[id];
		FILE		*q	= NULL;

		// For each rgn...

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( R.flag[ir] & mNewBad ) {

				if( !q ) {
					DskCreateDir( "Dropouts", stdout );
					char	buf[64];
					sprintf( buf, "Dropouts/drop_%d.txt", R.z );
					q = FileOpenOrDie( buf, "w" );
				}

				int	z, i, r;
				RealZIDR( z, i, r, iz, ir );

				if( R.flag[ir] & mOnIter ) {
					fprintf( q, "I %d.%d:%d\n", z, i, r );
					++D.iter;
				}
				else {
					fprintf( q, "P %d.%d:%d\n", z, i, r );
					++D.pnts;
				}
			}
		}

		if( q )
			fclose( q );
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* ScanEachLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScanEachLayer()
{
	int	nd = zihi - zilo + 1;

	vD.resize( nd );

	nthr = (zolo != zohi ? 16 : 1);

	if( nthr > nd )
		nthr = nd;

	if( !EZThreads( _Scan, nthr, 1, "_ScanEachLayer" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* GatherCounts -------------------------------------------------- */
/* --------------------------------------------------------------- */

void Dropout::GatherCounts()
{
	int	nd = zihi - zilo + 1;

	for( int id = 0; id < nd; ++id )
		Add( vD[id] );

	printf(
	"Worker %03d: ITER-DROPS %8ld PNTS-DROPS %8ld\n",
	wkid, iter, pnts );

	if( wkid > 0 )
		MPISend( this, sizeof(Dropout), 0, wkid );
	else {
		for( int iw = 1; iw < nwks; ++iw ) {

			Dropout	D;
			MPIRecv( &D, sizeof(Dropout), iw, iw );
			Add( D );

			printf(
			"Worker %03d: ITER-DROPS %8ld PNTS-DROPS %8ld\n",
			iw, D.iter, D.pnts );
		}

		printf(
		"--------------------------------------------------------\n"
		"     Total: ITER-DROPS %8ld PNTS-DROPS %8ld\n\n",
		iter, pnts );
	}
}

/* --------------------------------------------------------------- */
/* Scan ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Dropout::Scan()
{
	printf( "\n---- Dropouts ----\n" );

	clock_t	t0 = StartTiming();

	ScanEachLayer();
	GatherCounts();
	vD.clear();

	StopTiming( stdout, "Drops ", t0 );
}


