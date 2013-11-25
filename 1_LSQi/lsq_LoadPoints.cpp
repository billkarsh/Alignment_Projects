

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"

#include	"File.h"
#include	"Timer.h"

#include	<pthread.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CLoadPoints		*ME;
static pthread_mutex_t	mutex_vC = PTHREAD_MUTEX_INITIALIZER;






/* --------------------------------------------------------------- */
/* InitTablesToMaximum ------------------------------------------- */
/* --------------------------------------------------------------- */

static void InitTablesToMaximum()
{
	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		int	z = vL[iL].z;

		mZ[z] = iL;
		vR.push_back( Rgns( z ) );
	}
}

/* --------------------------------------------------------------- */
/* _Loader ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Loader( void* ithr )
{
	const int blksz = 1000;

	int	nmax = blksz;
	vector<CorrPnt>	vc( nmax );

	for( int j = (long)ithr; j < ME->njob; j += ME->nthr ) {

		CLoadPoints::CJob	&J = ME->vJ[j];

		char	buf[2048];
		sprintf( buf, "%s/%d/%c%d_%d/pts.%s",
		ME->tempdir, J.z, J.SorD, J.x, J.y,
		(J.SorD == 'S' ? "same" : "down") );

		FILE	*f = fopen( buf, "r" );

		if( f ) {

			int	n = 0;

			for(;;) {

				if( n >= nmax )
					vc.resize( nmax += blksz );

				if( vc[n].FromFile( f ) )
					++n;
				else
					break;
			}

			fclose( f );

			pthread_mutex_lock( &mutex_vC );
			vC.insert( vC.end(), vc.begin(), vc.begin() + n );
			pthread_mutex_unlock( &mutex_vC );
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* AppendJobs ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::AppendJobs(
	int			z,
	int			SorD,
	int			xhi,
	int			yhi )
{
	for( int y = 0; y <= yhi; ++y ) {

		for( int x = 0; x <= xhi; ++x )
			vJ.push_back( CJob( z, SorD, x, y ) );
	}
}

/* --------------------------------------------------------------- */
/* Remap --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Presorting before remapping was expected to improve
// cache performance, but in practice it doesn't matter.
// --- Certainly the points are very well sorted on disk.
//
static bool Sort_vC_inc( const CorrPnt& A, const CorrPnt& B )
{
	if( A.z1 < B.z1 )
		return true;
	if( A.z1 > B.z1 )
		return false;
	if( A.i1 < B.i1 )
		return true;
	if( A.i1 > B.i1 )
		return false;
	if( A.r1 < B.r1 )
		return true;
	if( A.r1 > B.r1 )
		return false;

	if( A.z2 < B.z2 )
		return true;
	if( A.z2 > B.z2 )
		return false;
	if( A.i2 < B.i2 )
		return true;
	if( A.i2 > B.i2 )
		return false;

	return A.r2 < B.r2;
}


void CLoadPoints::Remap()
{
	clock_t	t0 = StartTiming();

	int	nc = vC.size();

	for( int i = 0; i < nc; ++i ) {

		CorrPnt&	C = vC[i];

		MapZPair( C.z1, C.z2, C.z1, C.z2 );

		Rgns&	R1 = vR[C.z1];
		Rgns&	R2 = vR[C.z2];

		C.i1 = R1.Map( C.i1, C.r1 );
		C.i2 = R2.Map( C.i2, C.r2 );

		C.used = false;
	}

	StopTiming( stdout, "Mapping", t0 );
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::Load( const char *tempdir, bool isstack )
{
	printf( "\n---- Loading points ----\n" );

	clock_t	t0 = StartTiming();

	InitTablesToMaximum();

// Create list of reader file specs

	ME				= this;
	this->tempdir	= tempdir;

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		const Layer&	L = vL[iL];

		AppendJobs( L.z, 'S', L.sx, L.sy );

		if( isstack )
			AppendJobs( L.z, 'D', L.dx, L.dy );
	}

	njob = vJ.size();

// Create reader threads to scan points
// I will be thread zero.

	nthr = (isstack ? 16 : 2);

	if( nthr > njob )
		nthr = njob;

	vector<pthread_t>	vthr( nthr );

	if( nthr > 1 ) {

		pthread_attr_t	attr;
		pthread_attr_init( &attr );
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

		for( int i = 1; i < nthr; ++i ) {

			int	ret =
			pthread_create( &vthr[i], &attr, _Loader, (void*)i );

			if( ret ) {
				printf(
				"Error %d starting _Loader thread %d\n", ret, i );
				for( int j = 1; j < i; ++j )
					pthread_cancel( vthr[j] );
				exit( 42 );
			}
		}

		pthread_attr_destroy( &attr );
	}

// Do my own work

	_Loader( 0 );

// Join/wait my coworkers

	if( nthr > 1 ) {

		void*	ret;

		for( int i = 1; i < nthr; ++i )
			pthread_join( vthr[i], &ret );
	}

	vJ.clear();	// in case caller doesn't delete loader

	printf( "Loaded %d point pairs.\n", vC.size() );

	StopTiming( stdout, "Loading", t0 );

// Postprocessing

	Remap();
}


