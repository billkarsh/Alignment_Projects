

#include	"lsq_Globals.h"
#include	"lsq_XArray.h"

#include	"PipeFiles.h"
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






/* --------------------------------------------------------------- */
/* Load_AFromIDB ------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Load_AFromIDB()
{
	clock_t	t0 = StartTiming();

	int	nz = vR.size();

	X.resize( nz );

// For each z...

	for( int iz = 0; iz < nz; ++iz ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = X[iz];

		x.resize( vR[iz].nr * 6 );

		// Set used flags and load tforms from IDB.
		// We use the map (m) to traverse {id,r}.

		map<int,int>::iterator	it, nx, en = R.m.end();

		for( it = R.m.begin(); it != en; ++it ) {

			const Til2Img*	t2i;
			int				id = it->first,
							j0 = it->second,
							jlim;

			nx		= it;
			jlim	= (++nx != en ? nx->second : R.nr);

			// Propagate rgn #1 tform to all slots

			if( !IDBT2ICacheNGet1( t2i, idb, R.z, id ) )
				t2i = NULL;

			for( int j = j0; j < jlim; ++j ) {

				if( R.pts[j].size() >= 3 && t2i ) {

					t2i->T.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}

	StopTiming( stdout, "AFromIDB", t0 );
}


