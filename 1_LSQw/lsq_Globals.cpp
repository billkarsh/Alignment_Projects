

#include	"lsq_Globals.h"
#include	"lsq_MPI.h"

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
int				zLlo, zLhi,	// LHS neib needs these
				zRlo, zRhi,	// RHS neib needs these
				zilo, zihi,
				zolo, zohi,	// index into vR
				gW,   gH;	// image dims
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
	flag.resize( nr, fAsRead );
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
/* -------------------------------------- */
/* Get IDB path from temp/imageparams.txt */
/* -------------------------------------- */

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

/* ----------------------- */
/* Get image dims from IDB */
/* ----------------------- */

	if( !IDBGetImageDims( gW, gH, idb ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* InitTables ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Maximally size indexing tables.
//
void InitTables( int argzilo, int argzihi )
{
	printf( "\n---- Initializing tables ----\n" );

	clock_t	t0 = StartTiming();

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		int	z = vL[iL].z;

		mZ[z] = iL;
		vR.push_back( Rgns( z ) );
	}

// Turn z-ranges into zero-based indices into vR

	if( wkid > 0 ) {
		zLlo = mZ.find( zLlo )->second;
		zLhi = mZ.find( zLhi )->second;
	}

	if( wkid < nwks - 1 ) {
		zRlo = mZ.find( zRlo )->second;
		zRhi = mZ.find( zRhi )->second;
	}

	MapZPair( zilo, zihi, argzilo, argzihi );
	zolo = 0;
	zohi = nL - 1;

	StopTiming( stdout, "Table", t0 );
}

/* --------------------------------------------------------------- */
/* MapZPair ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return true if both za,zb map to zo range.
//
bool MapZPair( int &ia, int &ib, int za, int zb )
{
	static map<int,int>::iterator	en = mZ.end();

	static int	cached_za = -1, cached_zb = -1,
				cached_ia,      cached_ib;

	map<int,int>::iterator	mi;

	if( za != cached_za ) {

		mi = mZ.find( za );

		if( mi == en )
			return false;

		cached_ia = mi->second;
		cached_za = za;
	}

	ia = cached_ia;

	if( za == zb )
		ib = ia;
	else {

		if( zb != cached_zb ) {

			mi = mZ.find( zb );

			if( mi == en )
				return false;

			cached_ib = mi->second;
			cached_zb = zb;
		}

		ib = cached_ib;
	}

	return true;
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

// After loading points and calling this function:
// - All data indexed by zero-based (iz,idx0).
// - Initial value of all vC[i].used set per connectivity.
// - Initial value of all vR[i][j].flag set fAsRead.
//
// - vC[i].used are modified later as outliers are detected
//		in solution iteration cycles.
// - vR[i][j].flag are modified later as starting solutions
//		are loaded, and as needed during solve cycles.
//
void RemapIndices()
{
	clock_t	t0 = StartTiming();

	int	nc = vC.size();

	for( int i = 0; i < nc; ++i ) {

		CorrPnt&	C = vC[i];

		if( MapZPair( C.z1, C.z2, C.z1, C.z2 ) ) {

			Rgns&	R1 = vR[C.z1];
			Rgns&	R2 = vR[C.z2];

			C.i1 = R1.Map( C.i1, C.r1 );
			C.i2 = R2.Map( C.i2, C.r2 );

			C.used = (C.z1 >= zilo && C.z1 <= zihi)
						|| (C.z2 >= zilo && C.z2 <= zihi);

			R1.pts[C.i1].push_back( i );
			R2.pts[C.i2].push_back( i );
		}
		else
			C.used = false;
	}

	StopTiming( stdout, "Remap", t0 );
}

/* --------------------------------------------------------------- */
/* RealZIDR ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Convert zero-based constraint indices {iz, idx0}
// back into real world indices {z, id, r}.
//
void RealZIDR( int &z, int &id, int &r, int iz, int idx )
{
	Rgns& R = vR[iz];

	z = R.z;

	map<int,int>::iterator	mi, mn, en = R.m.end();

	for( mi = R.m.begin(); ; mi = mn ) {

		mn = mi;

		if( ++mn == en || mn->second > idx ) {

			id = mi->first;
			r  = (idx - mi->second) + 1;
			return;
		}
	}
}


