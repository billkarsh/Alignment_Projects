

#include	"lsq_Bounds.h"
#include	"lsq_Globals.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"TAffine.h"
#include	"THmgphy.h"
#include	"Timer.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const XArray	*gX;
static vector<DBox>	vB;
static int			nthr;






/* --------------------------------------------------------------- */
/* _Bounds ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Bounds( void* ithr )
{
	int				nb = zihi - zilo + 1;
	vector<Point>	cnr;
	Set4Corners( cnr, gW, gH );

// For each box...

	for( int ib = (long)ithr; ib < nb; ib += nthr ) {

		int						iz	= ib + zilo;
		const Rgns&				R	= vR[iz];
		const vector<double>&	x	= gX->X[iz];
		DBox&					B	= vB[ib];

		B.L = BIGD;
		B.R = -BIGD;
		B.B = BIGD;
		B.T = -BIGD;

		// For each rgn...

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( !FLAG_ISUSED( R.flag[ir] ) )
				continue;

			vector<Point>	c( 4 );
			memcpy( &c[0], &cnr[0], 4*2*sizeof(double) );

			if( gX->NE == 6 )
				X_AS_AFF( x, ir ).Transform( c );
			else
				X_AS_HMY( x, ir ).Transform( c );

			for( int k = 0; k < 4; ++k ) {
				B.L = fmin( B.L, c[k].x );
				B.R = fmax( B.R, c[k].x );
				B.B = fmin( B.B, c[k].y );
				B.T = fmax( B.T, c[k].y );
			}
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* CalcLayerwiseBoxes -------------------------------------------- */
/* --------------------------------------------------------------- */

// Set each vB[i] = bounds of gX->X[i].
//
// Only need range zi.
//
static void CalcLayerwiseBoxes( const XArray &X )
{
	gX = &X;

	int	nb = zihi - zilo + 1;

	vB.resize( nb );

	nthr = maxthreads;

	if( nthr > nb )
		nthr = nb;

	if( !EZThreads( _Bounds, nthr, 1, "_Bounds" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* GlobalBounds -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GlobalBounds( DBox &B )
{
// Calc my local overall box

	int	nb = zihi - zilo + 1;

	B = vB[0];

	for( int ib = 1; ib < nb; ++ib ) {

		const DBox&	b = vB[ib];

		B.L = fmin( B.L, b.L );
		B.R = fmax( B.R, b.R );
		B.B = fmin( B.B, b.B );
		B.T = fmax( B.T, b.T );
	}

// Send my box to master; then get global box from master

	if( wkid > 0 ) {
		MPISend( &B, sizeof(DBox), 0, wkid );
		MPIRecv( &B, sizeof(DBox), 0, wkid );
	}
	else if( nwks > 1 ) {

		for( int iw = 1; iw < nwks; ++iw ) {

			DBox	b;
			MPIRecv( &b, sizeof(DBox), iw, iw );
			B.L = fmin( B.L, b.L );
			B.R = fmax( B.R, b.R );
			B.B = fmin( B.B, b.B );
			B.T = fmax( B.T, b.T );
		}

		for( int iw = 1; iw < nwks; ++iw )
			MPISend( &B, sizeof(DBox), iw, iw );
	}
}

/* --------------------------------------------------------------- */
/* _Apply -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Apply( void* ithr )
{
	double	xorg	= -vB[0].L,
			yorg	= -vB[0].B;
	int		nL		= zihi - zilo + 1;
	THmgphy	M( 1,0,xorg, 0,1,yorg, 0,0 );

// For each layer...

	for( int iL = (long)ithr; iL < nL; iL += nthr ) {

		int						iz	= iL + zilo;
		const Rgns&				R	= vR[iz];
		const vector<double>&	x	= gX->X[iz];

		// For each rgn...

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( !FLAG_ISUSED( R.flag[ir] ) )
				continue;

			if( gX->NE == 6 )
				X_AS_AFF( x, ir ).AddXY( xorg, yorg );
			else {
				THmgphy&	T = X_AS_HMY( x, ir );
				T = M * T;
			}
		}
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* Apply --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Apply()
{
	int	nb = zihi - zilo + 1;

	nthr = maxthreads;

	if( nthr > nb )
		nthr = nb;

	if( !EZThreads( _Apply, nthr, 1, "_Apply" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Calculate B = global image bounds.
//
void Bounds( DBox &B, const XArray &X )
{
	clock_t	t0 = StartTiming();

	CalcLayerwiseBoxes( X );
	GlobalBounds( B );

// Adjust all transform origins

	vB[0] = B;
	Apply();
	vB.clear();

// Report pretty enclosing bounds

	B.R = ceil( B.R - B.L + 1 );
	B.T = ceil( B.T - B.B + 1 );
	B.L = 0;
	B.B = 0;

	if( !wkid ) {
		printf( "\nGlobal bounds: x=[0 %.2f] y=[0 %.2f].\n",
		B.R, B.T );
	}

	StopTiming( stdout, "Bounds", t0 );
}


