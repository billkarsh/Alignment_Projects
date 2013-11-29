

#include	"lsq_Globals.h"

#include	"File.h"
#include	"PipeFiles.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

string			idb;
vector<Layer>	vL;
map<int,int>	mZ;
vector<Rgns>	vR;
vector<CorrPnt>	vC;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Rgns::Rgns ---------------------------------------------------- */
/* --------------------------------------------------------------- */

Rgns::Rgns( int z ) : z(z), cached_id(-1)
{
	nr = IDBGetIDRgnMap( m, idb, z );
	pts.resize( nr );
	used.resize( nr, 0 );
}

/* --------------------------------------------------------------- */
/* Rgns::Map ----------------------------------------------------- */
/* --------------------------------------------------------------- */

int Rgns::Map( int id, int r )
{
	if( id != cached_id ) {
		cached_i  = m.find( id )->second + r - 1;
		cached_id = id;
	}

	return cached_i;
}

/* --------------------------------------------------------------- */
/* GetIDB -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetIDB( const char *tempdir )
{
	char	buf[2048];
	sprintf( buf, "%s/imageparams.txt", tempdir );

	FILE	*f = fopen( buf, "r" );

	if( f ) {

		CLineScan	LS;

		while( LS.Get( f ) > 0 ) {

			if( 1 == sscanf( LS.line, "IDBPATH %[^\n]", buf ) ) {

				idb = buf;
				goto close;
			}
		}

		printf( "Globals: imageparams.txt missing IDBPATH tag.\n" );

close:
		fclose( f );
	}
	else
		printf( "Globals: Can't open imageparams.txt.\n" );

	if( idb.empty() )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* InitTablesToMaximum ------------------------------------------- */
/* --------------------------------------------------------------- */

void InitTablesToMaximum()
{
	clock_t	t0 = StartTiming();

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		int	z = vL[iL].z;

		mZ[z] = iL;
		vR.push_back( Rgns( z ) );
	}

	StopTiming( stdout, "Table", t0 );
}

/* --------------------------------------------------------------- */
/* MapZPair ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MapZPair( int &ia, int &ib, int za, int zb )
{
	static int	cached_za = -1, cached_zb = -1,
				cached_ia,      cached_ib;

	if( za != cached_za ) {
		cached_ia = mZ.find( za )->second;
		cached_za = za;
	}

	ia = cached_ia;

	if( za == zb )
		ib = ia;
	else {

		if( zb != cached_zb ) {
			cached_ib = mZ.find( zb )->second;
			cached_zb = zb;
		}

		ib = cached_ib;
	}
}

/* --------------------------------------------------------------- */
/* Sort_vC_inc --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Presorting before remapping was expected to improve
// cache performance, but in practice it doesn't matter.
// --- Certainly the points are very well sorted on disk.
//
static bool Sort_vC_inc( const CorrPnt& A, const CorrPnt& B )
{
	if( A.z1 < B.z1 )
		return true;
	if( A.z1 > B.z1 )
		return false;
	if( A.i1 < B.i1 )
		return true;
	if( A.i1 > B.i1 )
		return false;
	if( A.r1 < B.r1 )
		return true;
	if( A.r1 > B.r1 )
		return false;

	if( A.z2 < B.z2 )
		return true;
	if( A.z2 > B.z2 )
		return false;
	if( A.i2 < B.i2 )
		return true;
	if( A.i2 > B.i2 )
		return false;

	return A.r2 < B.r2;
}

/* --------------------------------------------------------------- */
/* RemapIndices -------------------------------------------------- */
/* --------------------------------------------------------------- */

void RemapIndices()
{
	clock_t	t0 = StartTiming();

	int	nc = vC.size();

	for( int i = 0; i < nc; ++i ) {

		CorrPnt&	C = vC[i];

		MapZPair( C.z1, C.z2, C.z1, C.z2 );

		Rgns&	R1 = vR[C.z1];
		Rgns&	R2 = vR[C.z2];

		C.i1 = R1.Map( C.i1, C.r1 );
		C.i2 = R2.Map( C.i2, C.r2 );

		C.used = true;

		R1.pts[C.i1].push_back( i );
		R2.pts[C.i2].push_back( i );
	}

	StopTiming( stdout, "Remap", t0 );
}


