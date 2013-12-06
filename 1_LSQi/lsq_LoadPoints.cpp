

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"

#include	"Disk.h"
#include	"File.h"
#include	"EZThreads.h"
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

static CLoadPoints		*ME;
static pthread_mutex_t	mutex_fpnts = PTHREAD_MUTEX_INITIALIZER;






/* --------------------------------------------------------------- */
/* NameBinary ---------------------------------------------------- */
/* --------------------------------------------------------------- */

char* CLoadPoints::NameBinary( char *buf )
{
	sprintf( buf, "pnts_%d_%d_%d.dat", wkid, vR[zolo].z, vR[zohi].z );
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

	nthr = (zolo != zohi ? 16 : 2);

	if( nthr > njob )
		nthr = njob;

	if( !EZThreads( _Gather, nthr, 1, "_Gather" ) )
		exit( 42 );

	fclose( fpnts );
	pthread_mutex_destroy( &mutex_fpnts );
	vJ.clear();

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
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::Load( const char *tempdir, int wkid )
{
	printf( "\n---- Loading points ----\n" );

	clock_t	t0 = StartTiming();

	ME				= this;
	this->tempdir	= tempdir;
	this->wkid		= wkid;

	if( !IsBinary() )
		MakeBinary();

	LoadBinary();

	RemapIndices();

	StopTiming( stdout, "Total", t0 );

	printf( "Loaded %d point pairs.\n", vC.size() );
}


