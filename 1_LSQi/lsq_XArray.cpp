

#include	"lsq_Globals.h"
#include	"lsq_XArray.h"

#include	"EZThreads.h"
#include	"File.h"
#include	"PipeFiles.h"
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

static XArray*		ME;
static const char	*gpath;
static int			nthr;






/* --------------------------------------------------------------- */
/* _AFromIDB ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _AFromIDB( void* ithr )
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
/* _AFromScfTxt -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CIDRA {
public:
	TAffine	A;
	int		id, r;
public:
	inline bool FromFile( FILE *f )
	{
		return 8 == fscanf( f,
			" %d %d %lf %lf %lf %lf %lf %lf\n",
			&id, &r,
			&A.t[0], &A.t[1], &A.t[2],
			&A.t[3], &A.t[4], &A.t[5] );
	};
};


static bool Read_vA( vector<CIDRA> &vA, int z )
{
	char	buf[2048];
	FILE	*f;
	CIDRA	A;
	bool	nf = true;	// default = no folds

	sprintf( buf, "%s/X_A_%d.txt", gpath, z );
	f = FileOpenOrDie( buf, "r" );

	while( A.FromFile( f ) ) {

		if( A.r > 1 )
			nf = false;

		vA.push_back( A );
	}

	fclose( f );

	return nf;
}


static void* _AFromScfTxt( void* ithr )
{
	for( int iz = (long)ithr; iz < vR.size(); iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		vector<CIDRA>	vA;
		int				nf = Read_vA( vA, R.z );

		x.resize( R.nr * 6 );

		int						na = vA.size();
		map<int,int>::iterator	en = R.m.end();

		if( nf ) {	// Propagate rgn #1 to all block members

			// For each transform in vA...

			for( int ia = 0; ia < na; ++ia ) {

				// Get its block start and limit {j0,jlim}

				const CIDRA&			A = vA[ia];
				map<int,int>::iterator	mi = R.m.find( A.id );
				int						j0, jlim;

				if( mi == en )
					continue;

				j0		= mi->second;
				jlim	= (++mi != en ? mi->second : R.nr);

				// Propagate rgn #1 tform to all block members

				for( int j = j0; j < jlim; ++j ) {

					if( R.pts[j].size() >= 3 ) {

						A.A.CopyOut( &x[j * 6] );
						R.used[j] = true;
					}
				}
			}
		}
		else {	// Move each A to its specified position

			// For each transform in vA...

			for( int ia = 0; ia < na; ++ia ) {

				// Get its location j

				const CIDRA&			A = vA[ia];
				map<int,int>::iterator	mi = R.m.find( A.id );
				int						j;

				if( mi == en )
					continue;

				j = mi->second + A.r - 1;

				if( R.pts[j].size() >= 3 ) {

					A.A.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* _AFromBin ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadXBin( vector<double> &x, int cAH, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/X_%c_%d.bin", gpath, cAH, z );
	f = FileOpenOrDie( buf, "rb" );
	fread( &x[0], sizeof(double), x.size(), f );
	fclose( f );
}


static void ReadUBin( vector<char> &u, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/U_%d.bin", gpath, z );
	f = FileOpenOrDie( buf, "rb" );
	fread( &u[0], sizeof(char), u.size(), f );
	fclose( f );
}


static void* _AFromBin( void* ithr )
{
	for( int iz = (long)ithr; iz < vR.size(); iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];

		x.resize( R.nr * 6 );
		ReadXBin( x, 'A', R.z );
		ReadUBin( R.used, R.z );

		for( int j = 0; j < R.nr; ++j )
			R.used[j] &= R.pts[j].size() >= 3;
	}
}

/* --------------------------------------------------------------- */
/* _HFromBin ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _HFromBin( void* ithr )
{
	for( int iz = (long)ithr; iz < vR.size(); iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];

		x.resize( R.nr * 8 );
		ReadXBin( x, 'H', R.z );
		ReadUBin( R.used, R.z );

		for( int j = 0; j < R.nr; ++j )
			R.used[j] &= R.pts[j].size() >= 4;
	}
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Load( const char *path )
{
	clock_t	t0 = StartTiming();

	int	nz = vR.size();

	X.resize( nz );

	ME		= this;
	gpath	= path;
	nthr	= (zolo != zohi ? 16 : 1);

	if( nthr > nz )
		nthr = nz;

	EZThreadproc	proc;
	const char		*sproc;

	if( !path ) {
		NE		= 6;
		proc	= _AFromIDB;
		sproc	= "AFromIDB";
	}
	else {

		const char *name = FileNamePtr( path );

		if( strstr( name, "X_A" ) ) {

			NE = 6;

			if( strstr( name, "X_A_TXT" ) ) {
				proc	= _AFromScfTxt;
				sproc	= "AFromScfTxt";
			}
			else if( strstr( name, "X_A_BIN" ) ) {
				proc	= _AFromBin;
				sproc	= "AFromBin";
			}
			else
				goto error;
		}
		else if( strstr( name, "X_H" ) ) {

			NE = 8;

			if( strstr( name, "X_H_BIN" ) ) {
				proc	= _HFromBin;
				sproc	= "HFromBin";
			}
			else
				goto error;
		}
		else {
error:
			printf( "Unknown prior type [%s].\n", name );
			exit( 42 );
		}
	}

	if( !EZThreads( proc, nthr, 1, sproc ) )
		exit( 42 );

	StopTiming( stdout, sproc, t0 );
}


