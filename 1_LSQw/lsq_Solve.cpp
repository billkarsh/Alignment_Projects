

#include	"lsq_Solve.h"
#include	"lsq_Globals.h"

#include	"EZThreads.h"
#include	"LinEqu.h"
#include	"CRigid.h"
#include	"THmgphy.h"
#include	"Timer.h"

#include	<stdlib.h>
#include	<string.h>

#include	<algorithm>
using namespace std;


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
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	KILL( q )	vthr[(long)ithr].vkill.push_back( q )
#define	CUTD( q )	vthr[(long)ithr].vcutd.push_back( q )
#define	MARK( q )	vthr[(long)ithr].vmark.push_back( q )

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Todo {
// Which 0-based {iz,ir} a thread should work on
public:
	int	iz, ir;
public:
	Todo() {};
	Todo( int iz, int ir ) : iz(iz), ir(ir) {};
	bool First( int ithr );
	bool Next();
	static bool UseThreads( int minrgns );
	static int RgnCount();
};

class Thrdat {
// Each thread's edit tracking data
public:
	vector<Todo>	vmark,	// update pnt used flags
					vcutd,	// cut down pts
					vkill;	// can't rescue
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static double			Wr, Etol;
static XArray			*Xs, *Xd;
static vector<Thrdat>	vthr;
static int				regtype,
						editdelay,
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
/* SortPnts ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool SortPnts( int a, int b )
{
	const CorrPnt&	A = vC[a];
	const CorrPnt&	B = vC[b];

	if( A.z1 < B.z1 )
		return true;
	if( A.z1 > B.z1 )
		return false;
	if( A.i1 < B.i1 )
		return true;
	if( A.i1 > B.i1 )
		return false;

	if( A.z2 < B.z2 )
		return true;
	if( A.z2 > B.z2 )
		return false;
	if( A.i2 < B.i2 )
		return true;
	if( A.i2 > B.i2 )
		return false;

	if( A.p1.x < B.p1.x )
		return true;
	if( A.p1.x > B.p1.x )
		return false;
	if( A.p1.y < B.p1.y )
		return true;
	if( A.p1.y > B.p1.y )
		return false;

	if( A.p2.x < B.p2.x )
		return true;
	if( A.p2.x > B.p2.x )
		return false;

	return A.p2.y < B.p2.y;
}

/* --------------------------------------------------------------- */
/* ShortenList --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ShortenList( const Todo& Q, int ithr, int minpts )
{
	vector<int>		keep;
	vector<int>&	vp = vR[Q.iz].pts[Q.ir];
	int				np = vp.size(),
					nu;

	for( int ip = 0; ip < np; ++ip ) {

		if( vC[vp[ip]].used )
			keep.push_back( vp[ip] );
	}

	if( (nu = keep.size()) >= minpts ) {
		vp.resize( nu );
		memcpy( &vp[0], &keep[0], nu * sizeof(int) );
	}
	else
		KILL( Q );
}

/* --------------------------------------------------------------- */
/* Cut_A2A ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Cut_A2A( double *RHS, const Todo& Q, int ithr )
{
	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

// For rgn Q...

	CRigid				*rgd;
	const Rgns&			R  = vR[Q.iz];
	const vector<int>&	vp = R.pts[Q.ir];
	int					np = vp.size(),
						nu = 0;	// count pts used

	if( regtype == 'R' )
		rgd = new CRigid;
	else
		rgd = new CTrans;

	double		LHS[6*6];
	TAffine*	Ta = &X_AS_AFF( Xs->X[Q.iz], Q.ir );
	TAffine*	Tb;
	int			lastbi = -1;

	memset( RHS, 0, 6   * sizeof(double) );
	memset( LHS, 0, 6*6 * sizeof(double) );

	// For each of its points...

	for( int ip = 0; ip < np; ++ip ) {

		CorrPnt&	C = vC[vp[ip]];

		if( !C.used || C.z1 != C.z2 )
			continue;

		Point	A, B;

		// Which of {1,2} is the A-side?

		if( C.i1 == Q.ir ) {	// A is 1

			if( C.i2 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) )
					continue;

				Tb = &X_AS_AFF( Xs->X[C.z2], C.i2 );
				lastbi = C.i2;
			}

			Ta->Transform( A = C.p1 );
			Tb->Transform( B = C.p2 );

			if( pass > editdelay && A.DistSqr( B ) > Etol )
				continue;

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p1;

			rgd->Add( A, B );
		}
		else {	// A is 2

			if( C.i1 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) )
					continue;

				Tb = &X_AS_AFF( Xs->X[C.z1], C.i1 );
				lastbi = C.i1;
			}

			Ta->Transform( A = C.p2 );
			Tb->Transform( B = C.p1 );

			if( pass > editdelay && A.DistSqr( B ) > Etol )
				continue;

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p2;

			rgd->Add( A, B );
		}

		++nu;

		double	v[3] = { A.x, A.y, 1.0 };

		AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
		AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
	}

	if( nu < 3 || !Solve_Quick( LHS, RHS, 6 ) )
		KILL( Q );
	else {

		rgd->Regularize( RHS, 6, Wr );

		if( X_AS_AFF( RHS, 0 ).Squareness() > sqrtol )
			KILL( Q );
		else
			CUTD( Q );
	}

	delete rgd;
}

/* --------------------------------------------------------------- */
/* Cut_A2H ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Cut_A2H( double *RHS, const Todo& Q, int ithr )
{
	int	i1[5] = { 0, 1, 2, 6, 7 },
		i2[5] = { 3, 4, 5, 6, 7 };

// For rgn Q...

	CRigid				*rgd;
	const Rgns&			R  = vR[Q.iz];
	const vector<int>&	vp = R.pts[Q.ir];
	int					np = vp.size(),
						nu = 0;	// count pts used

	if( regtype == 'R' )
		rgd = new CRigid;
	else
		rgd = new CTrans;

	double		LHS[8*8];
	TAffine*	Ta = &X_AS_AFF( Xs->X[Q.iz], Q.ir );
	TAffine*	Tb;
	int			lastbi = -1;

	memset( RHS, 0, 8   * sizeof(double) );
	memset( LHS, 0, 8*8 * sizeof(double) );

	// For each of its points...

	for( int ip = 0; ip < np; ++ip ) {

		const CorrPnt&	C = vC[vp[ip]];

		if( !C.used || C.z1 != C.z2 )
			continue;

		Point	A, B;

		// Which of {1,2} is the A-side?

		if( C.i1 == Q.ir ) {	// A is 1

			if( C.i2 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) )
					continue;

				Tb = &X_AS_AFF( Xs->X[C.z2], C.i2 );
				lastbi = C.i2;
			}

			Ta->Transform( A = C.p1 );
			Tb->Transform( B = C.p2 );

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p1;

			rgd->Add( A, B );
		}
		else {	// A is 2

			if( C.i1 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) )
					continue;

				Tb = &X_AS_AFF( Xs->X[C.z1], C.i1 );
				lastbi = C.i1;
			}

			Ta->Transform( A = C.p2 );
			Tb->Transform( B = C.p1 );

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p2;

			rgd->Add( A, B );
		}

		++nu;

		double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

		AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

		v[3] = -A.x*B.y;
		v[4] = -A.y*B.y;

		AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
	}

	if( nu < 4 || !Solve_Quick( LHS, RHS, 8 ) )
		KILL( Q );
	else {

		rgd->Regularize( RHS, 8, Wr );

		if( X_AS_HMY( RHS, 0 ).Squareness() > sqrtol )
			KILL( Q );
		else
			CUTD( Q );
	}

	delete rgd;
}

/* --------------------------------------------------------------- */
/* Cut_H2H ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Cut_H2H( double *RHS, const Todo& Q, int ithr )
{
	int	i1[5] = { 0, 1, 2, 6, 7 },
		i2[5] = { 3, 4, 5, 6, 7 };

// For rgn Q...

	CRigid				*rgd;
	const Rgns&			R  = vR[Q.iz];
	const vector<int>&	vp = R.pts[Q.ir];
	int					np = vp.size(),
						nu = 0;	// count pts used

	if( regtype == 'R' )
		rgd = new CRigid;
	else
		rgd = new CTrans;

	double		LHS[8*8];
	THmgphy*	Ta = &X_AS_HMY( Xs->X[Q.iz], Q.ir );
	THmgphy*	Tb;
	int			lastbi = -1;

	memset( RHS, 0, 8   * sizeof(double) );
	memset( LHS, 0, 8*8 * sizeof(double) );

	// For each of its points...

	for( int ip = 0; ip < np; ++ip ) {

		CorrPnt&	C = vC[vp[ip]];

		if( !C.used || C.z1 != C.z2 )
			continue;

		Point	A, B;

		// Which of {1,2} is the A-side?

		if( C.i1 == Q.ir ) {	// A is 1

			if( C.i2 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) )
					continue;

				Tb = &X_AS_HMY( Xs->X[C.z2], C.i2 );
				lastbi = C.i2;
			}

			Ta->Transform( A = C.p1 );
			Tb->Transform( B = C.p2 );

			if( pass > editdelay && A.DistSqr( B ) > Etol )
				continue;

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p1;

			rgd->Add( A, B );
		}
		else {	// A is 2

			if( C.i1 != lastbi ) {

				if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) )
					continue;

				Tb = &X_AS_HMY( Xs->X[C.z1], C.i1 );
				lastbi = C.i1;
			}

			Ta->Transform( A = C.p2 );
			Tb->Transform( B = C.p1 );

			if( pass > editdelay && A.DistSqr( B ) > Etol )
				continue;

			B.x = Wb * B.x + (1 - Wb) * A.x;
			B.y = Wb * B.y + (1 - Wb) * A.y;
			A = C.p2;

			rgd->Add( A, B );
		}

		++nu;

		double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

		AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

		v[3] = -A.x*B.y;
		v[4] = -A.y*B.y;

		AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
	}

	if( nu < 4 || !Solve_Quick( LHS, RHS, 8 ) )
		KILL( Q );
	else {

		rgd->Regularize( RHS, 8, Wr );

		if( X_AS_HMY( RHS, 0 ).Squareness() > sqrtol )
			KILL( Q );
		else
			CUTD( Q );
	}

	delete rgd;
}

/* --------------------------------------------------------------- */
/* _A2A ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _A2A( void* ithr )
{
	Todo	Q;

	if( !Q.First( (long)ithr ) )
		return NULL;

	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

// For each of my rgns...

	do {

		CRigid			*rgd;
		Rgns&			R  = vR[Q.iz];
		vector<int>&	vp = R.pts[Q.ir];
		int				np = vp.size(),
						nu = 0;	// count pts used

		if( np < 3 ) {
			KILL( Q );
			continue;
		}

		if( regtype == 'R' )
			rgd = new CRigid;
		else
			rgd = new CTrans;

		double		*RHS = X_AS_AFF( Xd->X[Q.iz], Q.ir ).t;
		double		LHS[6*6];
		TAffine*	Ta = &X_AS_AFF( Xs->X[Q.iz], Q.ir );
		TAffine*	Tb;
		int			lastbi,
					lastbz	= -1;

		memset( RHS, 0, 6   * sizeof(double) );
		memset( LHS, 0, 6*6 * sizeof(double) );

		// Sort the points so that cummulative rounding
		// error tends to be same independent of nwks.

		if( !pass )
			sort( vp.begin(), vp.end(), SortPnts );

		// For each of its points...

		for( int ip = 0; ip < np; ++ip ) {

			CorrPnt&	C = vC[vp[ip]];

			if( !C.used )
				continue;

			Point	A, B;

			// Which of {1,2} is the A-side?

			if( C.z1 == Q.iz && C.i1 == Q.ir ) {	// A is 1

				if( C.z2 != lastbz ) {
					lastbz = C.z2;
					lastbi = -1;
				}

				if( C.i2 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get such points marked 'not used'.

						MARK( Todo( C.z2, C.i2 ) );
						continue;
					}

					Tb = &X_AS_AFF( Xs->X[C.z2], C.i2 );
					lastbi = C.i2;
				}

				Ta->Transform( A = C.p1 );
				Tb->Transform( B = C.p2 );

				if( pass > editdelay && A.DistSqr( B ) > Etol )
					continue;

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p1;

				rgd->Add( A, B );
			}
			else {	// A is 2

				if( C.z1 != lastbz ) {
					lastbz = C.z1;
					lastbi = -1;
				}

				if( C.i1 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get such points marked 'not used'.

						MARK( Todo( C.z1, C.i1 ) );
						continue;
					}

					Tb = &X_AS_AFF( Xs->X[C.z1], C.i1 );
					lastbi = C.i1;
				}

				Ta->Transform( A = C.p2 );
				Tb->Transform( B = C.p1 );

				if( pass > editdelay && A.DistSqr( B ) > Etol )
					continue;

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p2;

				rgd->Add( A, B );
			}

			++nu;

			double	v[3] = { A.x, A.y, 1.0 };

			AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
			AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
		}

		if( nu < 3 )
			KILL( Q );
		else if( !Solve_Quick( LHS, RHS, 6 ) )
			Cut_A2A( RHS, Q, (long)ithr );
		else {

			rgd->Regularize( RHS, 6, Wr );

			if( pass >= editdelay
				&& X_AS_AFF( RHS, 0 ).Squareness() > sqrtol ) {

				Cut_A2A( RHS, Q, (long)ithr );
			}
			else if( nu < np )
				ShortenList( Q, (long)ithr, 3 );
		}

		delete rgd;

	} while( Q.Next() );

	return NULL;
}

/* --------------------------------------------------------------- */
/* _A2H ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Solve for H using point-pairs {A,B}: B ~ HA.
//
// |Bx'|   |a b c|   |Ax|
// |By'| = |d e f| x |Ay| ; Bx=Bx'/Bw', By=By'/Bw'
// |Bw'|   |g h 1|   |1 |
//
// Expand:
// (gAx + hAy + 1)Bx = aAx + bAy + c
// (gAx + hAy + 1)By = dAx + eAy + f
//
// Gather:
// [Ax Ay 1 0 0 0 -AxBx -AyBx] [a b c d e f g h] = Bx
// [0 0 0 Ax Ay 1 -AxBy -AyBy] [a b c d e f g h] = By
//
static void* _A2H( void* ithr )
{
	Todo	Q;

	if( !Q.First( (long)ithr ) )
		return NULL;

	int	i1[5] = { 0, 1, 2, 6, 7 },
		i2[5] = { 3, 4, 5, 6, 7 };

// For each of my rgns...

	do {

		CRigid			*rgd;
		Rgns&			R  = vR[Q.iz];
		vector<int>&	vp = R.pts[Q.ir];
		int				np = vp.size(),
						nu = 0;	// count pts used

		if( np < 4 ) {
			KILL( Q );
			continue;
		}

		if( regtype == 'R' )
			rgd = new CRigid;
		else
			rgd = new CTrans;

		double		*RHS = X_AS_HMY( Xd->X[Q.iz], Q.ir ).t;
		double		LHS[8*8];
		TAffine*	Ta = &X_AS_AFF( Xs->X[Q.iz], Q.ir );
		TAffine*	Tb;
		int			lastbi,
					lastbz	= -1;

		memset( RHS, 0, 8   * sizeof(double) );
		memset( LHS, 0, 8*8 * sizeof(double) );

		// Sort the points so that cummulative rounding
		// error tends to be same independent of nwks.

		if( !pass )
			sort( vp.begin(), vp.end(), SortPnts );

		// For each of its points...

		for( int ip = 0; ip < np; ++ip ) {

			const CorrPnt&	C = vC[vp[ip]];

			if( !C.used )
				continue;

			Point	A, B;

			// Which of {1,2} is the A-side?

			if( C.z1 == Q.iz && C.i1 == Q.ir ) {	// A is 1

				if( C.z2 != lastbz ) {
					lastbz = C.z2;
					lastbi = -1;
				}

				if( C.i2 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get such points marked 'not used'.

						MARK( Todo( C.z2, C.i2 ) );
						continue;
					}

					Tb = &X_AS_AFF( Xs->X[C.z2], C.i2 );
					lastbi = C.i2;
				}

				Ta->Transform( A = C.p1 );
				Tb->Transform( B = C.p2 );

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p1;

				rgd->Add( A, B );
			}
			else {	// A is 2

				if( C.z1 != lastbz ) {
					lastbz = C.z1;
					lastbi = -1;
				}

				if( C.i1 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get those point marked 'not used'.

						MARK( Todo( C.z1, C.i1 ) );
						continue;
					}

					Tb = &X_AS_AFF( Xs->X[C.z1], C.i1 );
					lastbi = C.i1;
				}

				Ta->Transform( A = C.p2 );
				Tb->Transform( B = C.p1 );

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p2;

				rgd->Add( A, B );
			}

			++nu;

			double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

			AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

			v[3] = -A.x*B.y;
			v[4] = -A.y*B.y;

			AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
		}

		if( nu < 4 )
			KILL( Q );
		else if( !Solve_Quick( LHS, RHS, 8 ) )
			Cut_A2H( RHS, Q, (long)ithr );
		else {

			rgd->Regularize( RHS, 8, Wr );

			if( X_AS_HMY( RHS, 0 ).Squareness() > sqrtol ) {

				Cut_A2H( RHS, Q, (long)ithr );
			}
			else if( nu < np )
				ShortenList( Q, (long)ithr, 4 );
		}

		delete rgd;

	} while( Q.Next() );

	return NULL;
}

/* --------------------------------------------------------------- */
/* _H2H ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _H2H( void* ithr )
{
	Todo	Q;

	if( !Q.First( (long)ithr ) )
		return NULL;

	int	i1[5] = { 0, 1, 2, 6, 7 },
		i2[5] = { 3, 4, 5, 6, 7 };

// For each of my rgns...

	do {

		CRigid			*rgd;
		Rgns&			R  = vR[Q.iz];
		vector<int>&	vp = R.pts[Q.ir];
		int				np = vp.size(),
						nu = 0;	// count pts used

		if( np < 4 ) {
			KILL( Q );
			continue;
		}

		if( regtype == 'R' )
			rgd = new CRigid;
		else
			rgd = new CTrans;

		double		*RHS = X_AS_HMY( Xd->X[Q.iz], Q.ir ).t;
		double		LHS[8*8];
		THmgphy*	Ta = &X_AS_HMY( Xs->X[Q.iz], Q.ir );
		THmgphy*	Tb;
		int			lastbi,
					lastbz	= -1;

		memset( RHS, 0, 8   * sizeof(double) );
		memset( LHS, 0, 8*8 * sizeof(double) );

		// Sort the points so that cummulative rounding
		// error tends to be same independent of nwks.

		if( !pass )
			sort( vp.begin(), vp.end(), SortPnts );

		// For each of its points...

		for( int ip = 0; ip < np; ++ip ) {

			CorrPnt&	C = vC[vp[ip]];

			if( !C.used )
				continue;

			Point	A, B;

			// Which of {1,2} is the A-side?

			if( C.z1 == Q.iz && C.i1 == Q.ir ) {	// A is 1

				if( C.z2 != lastbz ) {
					lastbz = C.z2;
					lastbi = -1;
				}

				if( C.i2 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z2].flag[C.i2] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get those point marked 'not used'.

						MARK( Todo( C.z2, C.i2 ) );
						continue;
					}

					Tb = &X_AS_HMY( Xs->X[C.z2], C.i2 );
					lastbi = C.i2;
				}

				Ta->Transform( A = C.p1 );
				Tb->Transform( B = C.p2 );

				if( pass > editdelay && A.DistSqr( B ) > Etol )
					continue;

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p1;

				rgd->Add( A, B );
			}
			else {	// A is 2

				if( C.z1 != lastbz ) {
					lastbz = C.z1;
					lastbi = -1;
				}

				if( C.i1 != lastbi ) {

					if( !FLAG_ISUSED( vR[C.z1].flag[C.i1] ) ) {

						// Here's a pnt that's 'used' referencing
						// a rgn that's not. It's inefficient, so
						// we'll get those point marked 'not used'.

						MARK( Todo( C.z1, C.i1 ) );
						continue;
					}

					Tb = &X_AS_HMY( Xs->X[C.z1], C.i1 );
					lastbi = C.i1;
				}

				Ta->Transform( A = C.p2 );
				Tb->Transform( B = C.p1 );

				if( pass > editdelay && A.DistSqr( B ) > Etol )
					continue;

				B.x = Wb * B.x + (1 - Wb) * A.x;
				B.y = Wb * B.y + (1 - Wb) * A.y;
				A = C.p2;

				rgd->Add( A, B );
			}

			++nu;

			double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

			AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

			v[3] = -A.x*B.y;
			v[4] = -A.y*B.y;

			AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
		}

		if( nu < 4 )
			KILL( Q );
		else if( !Solve_Quick( LHS, RHS, 8 ) )
			Cut_H2H( RHS, Q, (long)ithr );
		else {

			rgd->Regularize( RHS, 8, Wr );

			if( pass >= editdelay
				&& X_AS_HMY( RHS, 0 ).Squareness() > sqrtol ) {

				Cut_H2H( RHS, Q, (long)ithr );
			}
			else if( nu < np )
				ShortenList( Q, (long)ithr, 4 );
		}

		delete rgd;

	} while( Q.Next() );

	return NULL;
}

/* --------------------------------------------------------------- */
/* UpdateFlags --------------------------------------------------- */
/* --------------------------------------------------------------- */

// A solving pass may encompass the following kinds edit
// operations on tracking and flag data:
//
// - (1) Mark a rgn flag as cut or killed.
// - (2) Change the used field of a CorrPnt.
// - (3) Shorten a rgn's point list.
//
// A solve pass has a multithreaded 'solve' phase followed by
// a single-threaded 'update' phase and that structure shapes
// how edits are managed.

// (1) In the solve phase we may not edit any of the shared input
// data, hence, we do not directly alter any flags or used fields.
// However, recommendations for rgn cuts or kills are enqueued for
// the update phase. In the update phase we actually do change the
// rgn flags, and in consequence we change the usage status of the
// points referenced by those rgns.
//
// (2) At this time we are not setting the used fields of points
// individually, as by some outlier identification scheme. Rather,
// it is the rgn flags that determine if a point is used. So the
// used fields get changed solely as described in (1). If lsq is
// restarted, there may be points that are initially 'used' but
// that reference rgns that are dead. The 'mark' facility updates
// the used fields for those cases.
//
// (3) As an efficiency measure, each rgn's list of points can be
// shortened to remove those that are not used. These lists are
// private to each thread so are edited in the solve phase.
//
static void UpdateFlags()
{
	int	ncut = 0, nkil = 0;
	int	nt = vthr.size();

// For each thread's lists...

	for( int it = 0; it < nt; ++it ) {

		vector<Todo>&	vmrk = vthr[it].vmark;
		vector<Todo>&	vcut = vthr[it].vcutd;
		vector<Todo>&	vkil = vthr[it].vkill;
		int				ne;	// n edits

		// Process marks

		ne = vmrk.size();

		for( int ie = 0; ie < ne; ++ie ) {

			const Todo&		e  = vmrk[ie];
			Rgns&			R  = vR[e.iz];
			vector<int>&	vp = R.pts[e.ir];
			int				np = vp.size();

			// mark all its pnts
			for( int ip = 0; ip < np; ++ip )
				vC[vp[ip]].used = false;
		}

		// Process cuts

		ne    = vcut.size();
		ncut += ne;

		for( int ie = 0; ie < ne; ++ie ) {

			const Todo&		e  = vcut[ie];
			Rgns&			R  = vR[e.iz];
			vector<int>&	vp = R.pts[e.ir];
			int				np = vp.size();

			// mark rgn cutd
			FLAG_ADDCUTD( R.flag[e.ir] );

			// mark its cross-pnts
			for( int ip = 0; ip < np; ++ip ) {

				CorrPnt& C = vC[vp[ip]];

				if( C.z1 != C.z2 )
					C.used = false;
			}
		}

		// Process kills

		ne    = vkil.size();
		nkil += ne;

		for( int ie = 0; ie < ne; ++ie ) {

			const Todo&		e  = vkil[ie];
			Rgns&			R  = vR[e.iz];
			vector<int>&	vp = R.pts[e.ir];
			int				np = vp.size();

			// mark rgn killed
			FLAG_ADDKILL( R.flag[e.ir] );

			// mark all its pnts
			for( int ip = 0; ip < np; ++ip )
				vC[vp[ip]].used = false;
		}
	}

// Report activity this pass

	if( ncut || nkil ) {
		printf( "Pass %d: tiles [cutd, killed] = [%d, %d].\n",
		pass, ncut, nkil );
	}
}

/* --------------------------------------------------------------- */
/* Do1Pass ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Do1Pass( EZThreadproc proc )
{
// multithreaded phase

	vthr.clear();
	vthr.resize( nthr );

	if( !EZThreads( proc, nthr, 1, "Solveproc" ) )
		exit( 42 );

// single-threaded phase

	UpdateFlags();
	vthr.clear();

// synchronize

	Xd->Updt();
}

/* --------------------------------------------------------------- */
/* SetSolveParams ------------------------------------------------ */
/* --------------------------------------------------------------- */

// When solving for affines, homographies or other transforms
// that include scaling, the objective function (minimization
// of square error) is satisfied by driving the scale toward
// zero. To counter that possibility, we regularize the system
// by solving both for the desired transform and for a model
// transform that does not permit scale change, say, a rigid
// transform. These two are combined in a weighted average:
// (1-Wr)*Aff + Wr*Rgd.
//
// errtol is the largest permitted point error.
//
void SetSolveParams( int type, double inWr, double inEtol )
{
	Etol	= inEtol * inEtol;
	Wr		= inWr;
	regtype	= type;
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
			editdelay	= 0;
			proc		= _A2H;
			cD			= 'H';
		}
	}
	else {
		editdelay	= 4;
		proc		= _H2H;
		cS			= 'H';
		cD			= 'H';
	}

	printf( "Solve: %c to %c (Wr %c, %g Etol %g iters %d)\n",
	cS, cD, regtype, Wr, sqrt( Etol ), iters );

/* -------- */
/* Set nthr */
/* -------- */

// Balance work/thread vs. thread overhead.
// Revisit if workload altered.

	const int minrgnsperthr = 8;

	int	nrgn = Todo::RgnCount();

	nthr = maxthreads;

	while( nthr > 1 && nrgn < minrgnsperthr * nthr )
		--nthr;

	printf( "NThrd: %d\n", nthr );

/* ------- */
/* Iterate */
/* ------- */

	Xs = &Xsrc;
	Xd = &Xdst;

	for( pass = 0; pass < iters; ++pass ) {

		Do1Pass( proc );

		// swap Xs<->Xd
		XArray	*Xt = Xs; Xs = Xd; Xd = Xt;

		if( !((long)DeltaSeconds( t0 ) % 300) )
			fflush( stdout );
	}

	StopTiming( stdout, "Solve", t0 );
}


