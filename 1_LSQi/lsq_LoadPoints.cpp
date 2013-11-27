

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"

#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"

#include	<limits.h>


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
static pthread_mutex_t	mutex_fpnts = PTHREAD_MUTEX_INITIALIZER;






/* --------------------------------------------------------------- */
/* NameBinary ---------------------------------------------------- */
/* --------------------------------------------------------------- */

char* CLoadPoints::NameBinary( char *buf )
{
	sprintf( buf, "pnts_%d_%d_%d.dat", wkid, zolo, zohi );
	return buf;
}

/* --------------------------------------------------------------- */
/* IsBinary ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CLoadPoints::IsBinary()
{
	char	buf[2048];
	return DskExists( NameBinary( buf ) );
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
/* _Gather ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Gather( void* ithr )
{
	const int		ngrow = 1000;
	int				nmax = ngrow;
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
					vc.resize( nmax += ngrow );

				if( vc[n].FromFile( f ) )
					++n;
				else
					break;
			}

			fclose( f );

			pthread_mutex_lock( &mutex_fpnts );
			fwrite( &vc[0], sizeof(CorrPnt), n, ME->fpnts );
			pthread_mutex_unlock( &mutex_fpnts );
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* MakeBinary ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::MakeBinary()
{
	clock_t	t0 = StartTiming();

// Output binary points file

	char	buf[2048];
	fpnts = FileOpenOrDie( NameBinary( buf ), "wb" );

// Create list of input file specs

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		const Layer&	L = vL[iL];

		AppendJobs( L.z, 'S', L.sx, L.sy );

		if( zolo != zohi )
			AppendJobs( L.z, 'D', L.dx, L.dy );
	}

	njob = vJ.size();

// Create reader threads to scan points
// I will be thread zero.

	nthr = (zolo != zohi ? 16 : 2);

	if( nthr > njob )
		nthr = njob;

	vector<pthread_t>	vthr( nthr );

	if( nthr > 1 ) {

		pthread_attr_t	attr;
		pthread_attr_init( &attr );
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );
		pthread_attr_setstacksize( &attr, PTHREAD_STACK_MIN );

		for( int i = 1; i < nthr; ++i ) {

			int	ret =
			pthread_create( &vthr[i], &attr, _Gather, (void*)i );

			if( ret ) {
				printf(
				"Error %d starting _Gather thread %d\n", ret, i );
				for( int j = 1; j < i; ++j )
					pthread_cancel( vthr[j] );
				exit( 42 );
			}
		}

		pthread_attr_destroy( &attr );
	}

// Do my own work

	_Gather( 0 );

// Join/wait my coworkers

	if( nthr > 1 ) {

		for( int i = 1; i < nthr; ++i ) {
			pthread_join( vthr[i], NULL );
			pthread_detach( vthr[i] );
		}
	}

	fclose( fpnts );
	pthread_mutex_destroy( &mutex_fpnts );

	StopTiming( stdout, "WrBin", t0 );
}

/* --------------------------------------------------------------- */
/* LoadBinary ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::LoadBinary()
{
	clock_t	t0 = StartTiming();

	char	buf[2048];
	long	n = (long)DskBytes( NameBinary( buf ) ) / sizeof(CorrPnt);

	vC.resize( n );

	FILE	*f = FileOpenOrDie( buf, "rb" );
	fread( &vC[0], sizeof(CorrPnt), n, f );
	fclose( f );

	StopTiming( stdout, "RdBin", t0 );
}

/* --------------------------------------------------------------- */
/* Remap --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::Load(
	const char	*tempdir,
	int			wkid,
	int			zolo,
	int			zohi )
{
	printf( "\n---- Loading points ----\n" );

	clock_t	t0 = StartTiming();

	InitTablesToMaximum();

	ME				= this;
	this->tempdir	= tempdir;
	this->wkid		= wkid;
	this->zolo		= zolo;
	this->zohi		= zohi;

	if( !IsBinary() )
		MakeBinary();

	LoadBinary();

	RemapIndices();

	StopTiming( stdout, "Total", t0 );

	printf( "Loaded %d point pairs.\n", vC.size() );
}


