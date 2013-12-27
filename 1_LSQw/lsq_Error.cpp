

#include	"lsq_Error.h"
#include	"lsq_Globals.h"

#include	"EZThreads.h"
#include	"Disk.h"
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
/* _Error -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Error( void* ithr )
{
	int	ns = vS.size();

// For each layer...

	for( int is = (long)ithr; is < ns; is += nthr ) {

		int						iz	= is + zilo;
		const Rgns&				Ra	= vR[iz];
		const vector<double>&	x	= gX->X[iz];

		// For each rgn...

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( !Ra.used[ir] )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const TAffine*		Ta = &X_AS_AFF( x, ir );
			const TAffine*		Tb;
			int					lastb	= -1;
			int					np		= P.size();
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


