

#include	"lsq_Solve.h"
#include	"lsq_Globals.h"

#include	"EZThreads.h"
#include	"TAffine.h"
#include	"THmgphy.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

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
static int		pass, nthr;






/* --------------------------------------------------------------- */
/* Todo::First --------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Todo::First( int ithr )
{
	ir = ithr;

	for( iz = zilo; iz <= zihi; ++iz ) {

		const Rgns&	R = vR[iz];

		for( ; ir < R.nr; ir += nthr ) {

			if( !R.flag[ir] )
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

			if( !R.flag[ir] )
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

			if( !f[ir] )
				++N;
		}
	}

	printf( "Solving for %d regions.\n", N );

	return N;
}

/* --------------------------------------------------------------- */
/* Do1Pass ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Do1Pass()
{
	const int minrgnsperthr = 100;

	int	nrgn = Todo::RgnCount();

	nthr = maxthreads;

	while( nthr > 1 && nrgn < minrgnsperthr * nthr )
		--nthr;
}

/* --------------------------------------------------------------- */
/* Solve --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Solve for tforms, first holding src tforms fixed to calculate
// a new set of dst tforms, then swapping these roles on each
// subsequent pass (by swapping Xs<->Xd pointers).
//
void Solve( XArray &Xsrc, XArray &Xdst, int iters )
{
	clock_t	t0 = StartTiming();

	Xs = &Xsrc;
	Xd = &Xdst;

	for( pass = 0; pass < iters; ++pass ) {

		Do1Pass();
		Xd->Updt();							// spread the love
		XArray	*Xt = Xs; Xs = Xd; Xd = Xt;	// swap Xs<->Xd
	}

	StopTiming( stdout, "Solve", t0 );
}


