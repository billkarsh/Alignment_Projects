

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"

#include	"File.h"
#include	"Timer.h"

#include	<pthread.h>
#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static int	nthr = 2;

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CLoadPoints	*ME;






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
	for( int ij = (long)ithr; ij < ME->nj; ij += nthr ) {

		CLoadPoints::CJob	&J = ME->vJ[ij];

		char	buf[2048];
		sprintf( buf, "%s/%d/%c%d_%d/pts.%s",
		ME->tempdir, J.z, J.SorD, J.x, J.y,
		(J.SorD == 'S' ? "same" : "down") );

		FILE	*f = fopen( buf, "r" );

		if( f ) {

			CLoadPoints::PntPair	P;

			while( 10 == fscanf( f, "CPOINT2"
				" %d.%d:%d %lf %lf"
				" %d.%d:%d %lf %lf\n",
				&P.z1, &P.d1, &P.r1, &P.x1, &P.y1,
				&P.z2, &P.d2, &P.r2, &P.x2, &P.y2 ) ) {

				J.vP.push_back( P );
			}

			fclose( f );
		}

		J.done = true;
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* LoadPoints ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::AppendJobs(
	int			z,
	int			SorD,
	int			xhi,
	int			yhi )
{
	if( xhi < 0 )
		return;

	for( int y = 0; y <= yhi; ++y ) {

		for( int x = 0; x <= xhi; ++x )
			vJ.push_back( CJob( z, SorD, x, y ) );
	}
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::Load( const char *tempdir, bool isstack )
{
	printf( "\n---- Loading points ----\n" );

	clock_t	t0 = StartTiming();

	ME				= this;
	this->tempdir	= tempdir;

	InitTablesToMaximum();

// define reader job params

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		const Layer&	L = vL[iL];

		AppendJobs( L.z, 'S', L.sx, L.sy );

		if( isstack )
			AppendJobs( L.z, 'D', L.dx, L.dy );
	}

	nj = vJ.size();

// make readers

	nthr = (isstack ? 16 : 2);

	if( nthr > nj )
		nthr = nj;

	vector<pthread_t>	vthr( nthr );

	pthread_attr_t	attr;
	pthread_attr_init( &attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	for( int i = 0; i < nthr; ++i ) {

		int	ret =
		pthread_create( &vthr[i], &attr, _Loader, (void*)i );

		if( ret ) {
			printf( "Error %d starting _Loader thread %d\n", ret, i );
			for( int j = 0; j < i; ++j )
				pthread_cancel( vthr[j] );
			exit( 42 );
		}
	}

	pthread_attr_destroy( &attr );

// Chase the readers and process the points

//	for(;;) {
//
//		for( int ij = 0; ij < nj; ++ij ) {
//
//			if( !vJ[ij].done ) {
//				sleep( 2 );
//				goto still_reading;
//			}
//		}
//
//		break;
//
//still_reading:;
//	}

// Join/wait my readers

	void*	ret;

	for( int i = 0; i < nthr; ++i )
		pthread_join( vthr[i], &ret );

	vJ.clear();

	StopTiming( stdout, "Load", t0 );
}


