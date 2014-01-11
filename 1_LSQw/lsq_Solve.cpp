

#include	"lsq_Solve.h"
#include	"lsq_Globals.h"

#include	"EZThreads.h"
#include	"TAffine.h"
#include	"THmgphy.h"
#include	"Timer.h"

#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Wb is a sort of old and new solution mixing parameter that
// favorably affects the rate of convergence. It's applied as
// follows: Ordinarily, point pairs {a,b} enter the equations
// as Ta'(a) = Tb(b), where (') marks the new solution. Here,
// we calculate Ta'(a) = Wb x Tb(b) + (1-Wb) x Ta(a). The value
// is chosen empirically and both same and down equations seem
// to like roughly the same value.

// sqrtol is tolerance for a new tform's squareness deviation,
// as calculated by {THmgphy,TAffine}.Squareness().

static const double Wb		= 0.9;
static const double	sqrtol	= sin( 15 * PI/180 );

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Todo {
// Which 0-based {iz,ir} a thread should work on
public:
	int	iz, ir;
public:
	bool First( int ithr );
	bool Next();
	static bool UseThreads( int minrgns );
	static int RgnCount();
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray	*Xs, *Xd;
static int		editdelay,
				pass, nthr;






/* --------------------------------------------------------------- */
/* Todo::First --------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Todo::First( int ithr )
{
	ir = ithr;

	for( iz = zilo; iz <= zihi; ++iz ) {

		const Rgns&	R = vR[iz];

		for( ; ir < R.nr; ir += nthr ) {

			if( FLAG_ISUSED( R.flag[ir] ) )
				return true;
		}

		ir -= R.nr;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* Todo::Next ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Todo::Next()
{
	ir += nthr;

	for( ; iz <= zihi; ++iz ) {

		const Rgns&	R = vR[iz];

		for( ; ir < R.nr; ir += nthr ) {

			if( FLAG_ISUSED( R.flag[ir] ) )
				return true;
		}

		ir -= R.nr;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* Todo::RgnCount ------------------------------------------------ */
/* --------------------------------------------------------------- */

int Todo::RgnCount()
{
	int	N = 0;

	for( int iz = zilo; iz <= zihi; ++iz ) {

		const vector<uint8>&	f  = vR[iz].flag;
		int						nr = vR[iz].nr;

		for( int ir = 0; ir < nr; ++ir ) {

			if( FLAG_ISUSED( f[ir] ) )
				++N;
		}
	}

	printf( "NRgns: %d\n", N );

	return N;
}

/* --------------------------------------------------------------- */
/* _A2A ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _A2A( void* ithr )
{
	return NULL;
}

/* --------------------------------------------------------------- */
/* _A2H ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _A2H( void* ithr )
{
	return NULL;
}

/* --------------------------------------------------------------- */
/* _H2H ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _H2H( void* ithr )
{
	return NULL;
}

/* --------------------------------------------------------------- */
/* Do1Pass ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Do1Pass( EZThreadproc proc )
{
	if( !EZThreads( proc, nthr, 1, "Solveproc" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* Solve --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Solve for tforms:
// On each pass, all of the src tforms are held fixed and new
// dst tforms are computed (and distributed across workers).
// Next the src/dst roles are swapped (Xs/Xd pointer swap) and
// the process repeats.
//
void Solve( XArray &Xsrc, XArray &Xdst, int iters )
{
	clock_t	t0 = StartTiming();

/* ---------------------- */
/* Mode specific settings */
/* ---------------------- */

	EZThreadproc	proc;
	int				cS, cD;

	if( Xsrc.NE == 6 ) {

		cS = 'A';

		if( Xdst.NE == 6 ) {
			editdelay	= 200;
			proc		= _A2A;
			cD			= 'A';
		}
		else {
			editdelay	= iters + 1;	// disable
			proc		= _A2H;
			cD			= 'H';
		}
	}
	else {
		editdelay	= 4;
		proc		= _H2H;
		cS			= 'A';
		cD			= 'H';
	}

	printf( "Solve: %c to %c iters %d\n", cS, cD, iters );

/* -------- */
/* Set nthr */
/* -------- */

	const int minrgnsperthr = 100;

	int	nrgn = Todo::RgnCount();

	nthr = maxthreads;

	while( nthr > 1 && nrgn < minrgnsperthr * nthr )
		--nthr;

/* ------- */
/* Iterate */
/* ------- */

	Xs = &Xsrc;
	Xd = &Xdst;

	for( pass = 0; pass < iters; ++pass ) {

		Do1Pass( proc );
		Xd->Updt();							// spread the love
		XArray	*Xt = Xs; Xs = Xd; Xd = Xt;	// swap Xs<->Xd
	}

	StopTiming( stdout, "Solve", t0 );
}


