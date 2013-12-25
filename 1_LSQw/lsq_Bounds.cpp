

#include	"lsq_Bounds.h"
#include	"lsq_Globals.h"
#include	"lsq_Msg.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
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

static const XArray		*gX;
static vector<DBox>		vB;
static int				nthr;






/* --------------------------------------------------------------- */
/* Cmd_CalcAndShare ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void Cmd_CalcAndShare()
{
	if( !wkid && nwks > 1 )
		MsgSend( "bounds-share" );
}

/* --------------------------------------------------------------- */
/* _Bounds ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Bounds( void* ithr )
{
	int				nb = vB.size();
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

			if( !R.used[ir] )
				continue;

			vector<Point>	c( 4 );
			memcpy( &c[0], &cnr[0], 4*2*sizeof(double) );
			X_AS_AFF( x, ir ).Transform( c );

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

	nthr = (zolo != zohi ? 16 : 2);

	if( nthr > nb )
		nthr = nb;

	if( !EZThreads( _Bounds, nthr, 1, "_Bounds" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteMyBox ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMyBox( DBox &B )
{
// Everyone get the box for local zi range

	int	nb = vB.size();

	B = vB[0];

	for( int ib = 1; ib < nb; ++ib ) {

		const DBox&	b = vB[ib];

		B.L = fmin( B.L, b.L );
		B.R = fmax( B.R, b.R );
		B.B = fmin( B.B, b.B );
		B.T = fmax( B.T, b.T );
	}

// Everyone write if multiple

	if( nwks <= 1 )
		return;

	DskCreateDir( "Bounds", stdout );

	char	buf[32];
	sprintf( buf, "Bounds/%d.txt", wkid );
	FILE	*f = FileOpenOrDie( buf, "w" );
	fprintf( f, "%.21f %.21f %.21f %.21f\n", B.L, B.R, B.B, B.T );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* Ack_CalcAndShare ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void Ack_CalcAndShare()
{
	if( wkid > 0 ) {

		if( !MsgWaitMatch( "bounds-share" ) )
			exit( 32 );

		MsgAck( wkid );
	}
}

/* --------------------------------------------------------------- */
/* Cmd_Bounds ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Zero syncs by waiting for all acks, then issues "bounds".
// All others wait for "bounds".
//
static void Cmd_Bounds()
{
	if( nwks <= 1 )
		return;

	if( !wkid ) {

		if( !MsgWaitAck( nwks ) )
			exit( 32 );

		MsgSend( "bounds" );
	}
	else {
		if( !MsgWaitMatch( "bounds" ) )
			exit( 32 );
	}
}

/* --------------------------------------------------------------- */
/* GlobalBounds -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GlobalBounds( DBox &B )
{
	for( int iw = 0; iw < nwks; ++iw ) {

		if( iw == wkid )
			continue;

		char	buf[32];
		sprintf( buf, "Bounds/%d.txt", iw );
		FILE	*f = FileOpenOrDie( buf, "r" );

		DBox	b;

		if( 4 != fscanf( f, "%lf%lf%lf%lf",
				&b.L, &b.R, &b.B, &b.T ) ) {

			exit( 32 );
		}

		B.L = fmin( B.L, b.L );
		B.R = fmax( B.R, b.R );
		B.B = fmin( B.B, b.B );
		B.T = fmax( B.T, b.T );

		fclose( f );
	}

	B.R = ceil( B.R - B.L + 1 );
	B.T = ceil( B.T - B.B + 1 );
	B.L = 0;
	B.B = 0;
}

/* --------------------------------------------------------------- */
/* Ack_Bounds ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Ack_Bounds()
{
	if( nwks <= 1 )
		return;

	if( !wkid ) {

		if( !MsgWaitAck( nwks ) )
			exit( 32 );
	}
	else
		MsgAck( wkid );
}

/* --------------------------------------------------------------- */
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Calculate B = global image bounds.
//
void Bounds( DBox &B, const XArray &X )
{
	printf( "\n---- Global bounds ----\n" );

	clock_t	t0 = StartTiming();

	Cmd_CalcAndShare();
		CalcLayerwiseBoxes( X );
		WriteMyBox( B );
		vB.clear();
	Ack_CalcAndShare();

	Cmd_Bounds();
		GlobalBounds( B );
	Ack_Bounds();

	if( !wkid ) {
		printf( "Global image bounds: x=[%f %f] y=[%f %f].\n",
		B.L, B.R, B.B, B.T );
	}

	StopTiming( stdout, "Bounds", t0 );
}


