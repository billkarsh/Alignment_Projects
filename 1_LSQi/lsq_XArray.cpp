

#include	"lsq_Globals.h"
#include	"lsq_XArray.h"

#include	"EZThreads.h"
#include	"PipeFiles.h"
#include	"Timer.h"

#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray*	ME;
static int		nthr;






/* --------------------------------------------------------------- */
/* _AFromIDB ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _AFromIDB( void* ithr )
{
	for( int iz = (long)ithr; iz < vR.size(); iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		vector<Til2Img>	t2i;

		// Get rgn #1 tforms

		if( !IDBT2IGet_JustIDandT( t2i, idb, R.z ) )
			exit( 42 );

		x.resize( R.nr * 6 );

		int						nt = t2i.size();
		map<int,int>::iterator	en = R.m.end();

		// For each transform in IDB...

		for( int it = 0; it < nt; ++it ) {

			// Get its block start and limit {j0,jlim}

			const Til2Img&			T = t2i[it];
			map<int,int>::iterator	mi = R.m.find( T.id );
			int						j0, jlim;

			if( mi == en )
				continue;

			j0		= mi->second;
			jlim	= (++mi != en ? mi->second : R.nr);

			// Propagate rgn #1 tform to all block members

			for( int j = j0; j < jlim; ++j ) {

				if( R.pts[j].size() >= 3 ) {

					T.T.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* Load_AFromIDB ------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Load_AFromIDB()
{
	clock_t	t0 = StartTiming();

	int	nz = vR.size();

	X.resize( nz );

	ME		= this;
	nthr	= (zolo != zohi ? 16 : 1);

	if( nthr > nz )
		nthr = nz;

	if( !EZThreads( _AFromIDB, nthr, 1, "_AFromIDB" ) )
		exit( 42 );

	StopTiming( stdout, "AFromIDB", t0 );
}


