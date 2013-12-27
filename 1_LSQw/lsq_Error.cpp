

#include	"lsq_Error.h"
#include	"lsq_Globals.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"TAffine.h"
#include	"THmgphy.h"
#include	"Timer.h"

#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class FileErr {
// Buffered file writing
private:
	static const int bufsize = 2048;
	FILE*			f;
	vector<float>	ve;
	int				n;
public:
	FileErr( int S_or_D, int z );
	virtual ~FileErr();
	void Add( double err );
};

class Stat {
public:
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const XArray	*gX;
static vector<Stat>	vS;
static int			nthr;






/* --------------------------------------------------------------- */
/* FileErr::FileErr ---------------------------------------------- */
/* --------------------------------------------------------------- */

FileErr::FileErr( int S_or_D, int z ) : f(NULL), n(0)
{
	if( S_or_D == 'S' || zolo != zohi ) {

		char	buf[32];
		sprintf( buf, "Error/Err_%c_%d.bin", S_or_D, z );
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
/* _Error -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Error( void* ithr )
{
	int	ns = vS.size();

// For each layer...

	for( int is = (long)ithr; is < ns; is += nthr ) {

		int						iz	= is + zilo;
		const Rgns&				Ra	= vR[iz];
		const vector<double>&	xa	= gX->X[iz];
		FileErr					FS( 'S', Ra.z ),
								FD( 'D', Ra.z );

		// For each rgn...

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( !Ra.used[ir] )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const TAffine*		Ta = &X_AS_AFF( xa, ir );
			const TAffine*		Tb;
			int					lastbi,
								lastbz	= -1,
								np		= P.size();

			// For each of its points...

			for( int ip = 0; ip < np; ++ip ) {

				CorrPnt&	C = vC[P[ip]];

				if( C.used <= 0 )
					continue;

				if( C.z1 == C.z2 ) {

					// Negate C.used to signify visited;
					// only needed for sames; thread-safe;
					// undone by CalcLayerwiseError().

					C.used = -1;

					Point	p1 = C.p1,
							p2 = C.p2;

					Ta->Transform( p1 );
					Ta->Transform( p2 );

					double	err = p2.DistSqr( p1 );

					FS.Add( err );
				}
				else if( zolo == zohi )
					continue;
				else if( C.z1 == iz ) {

					if( C.z2 != lastbz ) {
						lastbz = C.z2;
						lastbi = -1;
					}

					if( C.i2 != lastbi ) {

						if( !vR[C.z2].used[C.i2] )
							continue;

						Tb = &X_AS_AFF( gX->X[C.z2], C.i2 );
						lastbi = C.i2;
					}

					Point	pa = C.p1,
							pb = C.p2;

					Ta->Transform( pa );
					Tb->Transform( pb );

					double	err = pb.DistSqr( pa );

					FD.Add( err );
				}
				else
					continue;
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

	nthr = (zolo != zohi ? 16 : 1);

	if( nthr > ns )
		nthr = ns;

	if( !EZThreads( _Error, nthr, 1, "_Error" ) )
		exit( 42 );

// Restore C.used flags

	int	nc = vC.size();

	for( int ic = 0; ic < nc; ++ic ) {

		if( vC[ic].used < 0 )
			vC[ic].used = 1;
	}
}

/* --------------------------------------------------------------- */
/* Error --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Using the vC inliers (used = true) calculate several metrics
// from the residual errors:
//
// - Folder 'Error' with files-by-layer 'Err_S_i.bin' and
//		'Err_D_i.bin' with packed |err| values as floats.
//
void Error( const XArray &X )
{
	clock_t	t0 = StartTiming();

	DskCreateDir( "Error", stdout );
	CalcLayerwiseError( X );

	StopTiming( stdout, "Error", t0 );
}


