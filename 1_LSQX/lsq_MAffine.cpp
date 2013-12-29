

#include	"lsq_MTrans.h"
#include	"lsq_MAffine.h"

#include	"EZThreads.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"Timer.h"

#include	<math.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const double	SQRTOL		= sin( 15 * PI/180 );
const int		EDITDELAY	= 200;

/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SetPointPairs(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
{
	int	nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		double	fz =
		(vRgn[C.r1].z == vRgn[C.r2].z ? same_strength : 1);

		double	x1 = C.p1.x * fz / sc,
				y1 = C.p1.y * fz / sc,
				x2 = C.p2.x * fz / sc,
				y2 = C.p2.y * fz / sc;
		int		j  = vRgn[C.r1].itr * NX,
				k  = vRgn[C.r2].itr * NX;

		// A1(p1) - A2(p2) = 0

		double	v[6]  = { x1,  y1,  fz, -x2, -y2, -fz};
		int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
		int		i2[6] = {j+3, j+4, j+5, k+3, k+4, k+5};

		AddConstraint( LHS, RHS, 6, i1, v, 0.0 );
		AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
	}
}

/* --------------------------------------------------------------- */
/* SetIdentityTForm ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Explicitly set some TForm to Identity.
// @@@ Does it matter which one we use?
//
void MAffine::SetIdentityTForm(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				itr )
{
	double	stiff	= 1.0;

	double	one	= stiff;
	int		j	= itr * NX;

	AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;

// Don't do this with unite layer

	unite_layer = -1;

// Report which tile we set

	int	nr = vRgn.size();

	for( int k = 0; k < nr; ++k ) {

		if( vRgn[k].itr == itr ) {

			printf( "Ref region z=%d, id=%d\n",
			vRgn[k].z, vRgn[k].id );
			break;
		}
	}
}

/* --------------------------------------------------------------- */
/* SetUniteLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set one layer-full of TForms to those from a previous
// solution output file gArgs.unt_file.
//
void MAffine::SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
{
/* --------- */
/* Once only */
/* --------- */

	if( unite_layer < 0 )
		return;

/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

	map<MZIDR,TAffine>	M;

	LoadTAffineTbl_RngZ( M, unite_layer, unite_layer, unt_file );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

	double	stiff	= 0.001;

	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.z != unite_layer || R.itr < 0 )
			continue;

		map<MZIDR,TAffine>::iterator	it;

		it = M.find( MZIDR( R.z, R.id, R.rgn ) );

		if( it == M.end() )
			continue;

		double	one	= stiff,
				*t	= it->second.t;
		int		j	= R.itr * NX;

		AddConstraint( LHS, RHS, 1, &j, &one, one*t[0] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[1] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[2] / sc );	j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[3] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[4] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[5] / sc );	j++;
	}

/* --------- */
/* Once only */
/* --------- */

	unite_layer = -1;
}

/* --------------------------------------------------------------- */
/* SolveWithSquareness ------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SolveWithSquareness(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTr )
{
/* -------------------------- */
/* Add squareness constraints */
/* -------------------------- */

	double	stiff = square_strength;

	for( int i = 0; i < nTr; ++i ) {

		int	j = i * NX;

		// equal cosines
		{
			double	V[2] = {stiff, -stiff};
			int		I[2] = {j, j+4};

			AddConstraint( LHS, RHS, 2, I, V, 0.0 );
		}

		// opposite sines
		{
			double	V[2] = {stiff, stiff};
			int		I[2] = {j+1, j+3};

			AddConstraint( LHS, RHS, 2, I, V, 0.0 );
		}
	}

/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

	WriteSolveRead( X, LHS, RHS, "A-Square", nproc, false );
	PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* SolveWithUnitMag ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Effectively, we want to constrain the cosines and sines
// so that c^2 + s^2 = 1. We can't make constraints that are
// non-linear in the variables X[], but we can construct an
// approximation using the {c,s = X[]} of the previous fit:
// c*x + s*y = 1. To reduce sensitivity to the sizes of the
// previous fit c,s, we normalize them by m = sqrt(c^2 + s^2).
//
void MAffine::SolveWithUnitMag(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTR )
{
	double	stiff = scale_strength;

	for( int i = 0; i < nTR; ++i ) {

		int		j = i * NX;
		double	c = X[j];
		double	s = X[j+3];
		double	m = sqrt( c*c + s*s );

		// c*x/m + s*y/m = 1

		double	V[2] = {c * stiff, s * stiff};
		int		I[2] = {j, j+3};

		AddConstraint( LHS, RHS, 2, I, V, m * stiff );
	}

	WriteSolveRead( X, LHS, RHS, "A-Unimag", nproc, false );
	printf( "\t\t\t\t" );
	PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* RescaleAll ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::RescaleAll(
	vector<double>	&X,
	double			sc )
{
	int	nr	= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		X[itr+2] *= sc;
		X[itr+5] *= sc;
	}
}

/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::RotateAll(
	vector<double>	&X,
	double			degcw )
{
	TAffine	T, R;
	int		nr	= vRgn.size();

	R.SetCWRot( degcw, Point(0,0) );

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		TAffine	t( &X[itr] );

		T = R * t;
		T.CopyOut( &X[itr] );
	}
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::NewOriginAll(
	vector<double>	&X,
	double			xorg,
	double			yorg )
{
	int	nr	= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		X[itr+2] -= xorg;
		X[itr+5] -= yorg;
	}
}

/* --------------------------------------------------------------- */
/* DevFromTrans -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to see how much translation terms of each affine
// have moved from the trans-only starting values. We just list
// all tiles with dev > XXX, but do nothing with that for now.
//
void MAffine::DevFromTrans(
	const vector<double>	&T,
	const vector<double>	&X )
{
	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[i];

		if( I.itr < 0 )
			continue;

		const double	*J = &T[I.itr * 2],
						*K = &X[I.itr * NX];
		double			dx = J[0] - K[2],
						dy = J[1] - K[5];

		if( (dx = sqrt( dx*dx + dy*dy  )) > 200 ) {
			printf( "Dev: %d.%d:%d dr= %d\n",
			I.z, I.id, I.rgn, int(dx) );
		}
	}
}

/* --------------------------------------------------------------- */
/* DevFromPrior -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to see how much translation terms of each affine
// have moved from their prior starting values. We just list
// all tiles with dev > XXX, but do nothing with that for now.
//
void MAffine::DevFromPrior(
	const vector<double>	&A,
	const vector<double>	&X )
{
	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[i];

		if( I.itr < 0 )
			continue;

		const double	*J = &A[I.itr * NX],
						*K = &X[I.itr * NX];
		double			dx = J[2] - K[2],
						dy = J[5] - K[5];

		if( (dx = sqrt( dx*dx + dy*dy  )) > 200 ) {
			printf( "Dev: %d.%d:%d dr= %d\n",
			I.z, I.id, I.rgn, int(dx) );
		}
	}
}

/* --------------------------------------------------------------- */
/* LoadAffTable -------------------------------------------------- */
/* --------------------------------------------------------------- */

// z0 and nz resp. get first layer and num layers in scaffold.
//
void MAffine::LoadAffTable(
	vector<double>	&X,
	int				&z0,
	int				&nz,
	int				nTr )
{
	X.resize( nTr * NX );
	z0	= -1;
	nz	= 0;

// Load table

	printf( "Aff: Loading existing table.\n" );

	map<MZIDR,TAffine>	M;
	set<int>			Z;

	LoadTAffineTbl_AllZ( M, Z, priorafftbl );

// Data range

	if( !Z.size() ) {
		printf( "Aff: No layers in scaffold.\n" );
		exit( 42 );
	}

	z0 = *Z.begin();
	nz = *Z.rbegin();

	printf( "Aff: Z range [%d %d].\n", z0, nz );
	nz = nz - z0 + 1;

// Fill into X

	printf( "Aff: Mapping prior solutions.\n" );

	int	nr = vRgn.size(), nmapped = 0;

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		map<MZIDR,TAffine>::iterator	it;

		it = M.find( MZIDR( R.z, R.id, R.rgn ) );

		if( it != M.end() ) {

			memcpy( &X[R.itr*NX], it->second.t, NX*sizeof(double) );
			++nmapped;
		}
		else {
			// mark as no solution
			X[R.itr * NX] = 999.0;
			printf( "No prior for %d.%d:%d\n", R.z, R.id, R.rgn );
		}
	}

// Done

	printf( "Aff: Mapped %d affines.\n\n", nmapped );
	fflush( stdout );
}

/* --------------------------------------------------------------- */
/* UntwistScaffold ----------------------------------------------- */
/* --------------------------------------------------------------- */

class RgdSums {
public:
	double	Xa, Ya, Xb, Yb, XaXb, YaYb, XaYb, YaXb;
	int		za, zb, N;
};


// The standard way to obtain an external scaffold for use with
// the -prior option is to use XMLGetTF on HiRes.xml; the result
// of aligning low-res strips from very good montages. However,
// better angle calculations can be made later, after the down
// correspondence points are formed. Adjustments are calculated
// here as rigid transforms T{theta, kx, ky) from layer a to b.
//
// Let c=cos(theta), s=sin(theta), Sum over all point-pairs:
//
// E = Sum[Xb - cXa + sYa - kx]^2 + [Yb - sXa - cYa - ky]^2
//
// The params u = {theta,kx,ky) are determined by dE/du = 0.
// If we use notation [arg] => Sum[argi] over point-pairs,
//
// kx = ([Xb] - c[Xa] + s[Ya]) / N
// ky = ([Yb] - s[Xa] - c[Ya]) / N
//
//				[YaXb] - [XaYb] + ([Xa][Yb] - [Ya][Xb]) / N
// tan(theta) = --------------------------------------------
//				([Xa][Xb] + [Ya][Yb]) / N - [XaXb] - [YaYb]
//
void MAffine::UntwistScaffold(
	vector<double>	&X,
	int				z0,
	int				nz )
{
/* ------------------------------------------ */
/* Size sum vector for all layers in scaffold */
/* ------------------------------------------ */

	if( nz < 2 )
		return;

	vector<RgdSums>	vS( nz );

	memset( &vS[0], 0, nz*sizeof(RgdSums) );

/* ---------------------------- */
/* Do sums over all point pairs */
/* ---------------------------- */

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		const RGN	&Ra = vRgn[C.r1],
					&Rb = vRgn[C.r2];

		if( Ra.z != Rb.z + 1 )
			continue;

		if( X[Ra.itr * NX] == 999.0 )
			continue;

		if( X[Rb.itr * NX] == 999.0 )
			continue;

		Point	pa = C.p1,
				pb = C.p2;

		L2GPoint( pa, X, Ra.itr );
		L2GPoint( pb, X, Rb.itr );

		RgdSums	&S = vS[Ra.z - z0];

		S.Xa += pa.x;
		S.Ya += pa.y;
		S.Xb += pb.x;
		S.Yb += pb.y;

		S.XaXb += pa.x * pb.x;
		S.YaYb += pa.y * pb.y;
		S.XaYb += pa.x * pb.y;
		S.YaXb += pa.y * pb.x;

		S.za = Ra.z;
		S.zb = Rb.z;
		++S.N;
	}

/* --------------------------------------- */
/* Map cumulative rigids by affected layer */
/* --------------------------------------- */

	map<int,TAffine>	M;
	TAffine				Tprev;	// identity tform to start

	for( int i = 1; i < nz; ++i ) {

		RgdSums	&S = vS[i];

		if( S.N < 2 )
			continue;

		double
		theta = atan(
			(S.YaXb - S.XaYb + (S.Xa*S.Yb - S.Ya*S.Xb)/S.N) /
			((S.Xa*S.Xb + S.Ya*S.Yb)/S.N - S.XaXb - S.YaYb)
		),
		c  = cos( theta ),
		s  = sin( theta ),
		kx = (S.Xb - c*S.Xa + s*S.Ya) / S.N,
		ky = (S.Yb - s*S.Xa - c*S.Ya) / S.N;

		TAffine	T( c, -s, kx, s, c, ky );

		// propagate up the stack
		T = T * Tprev;
		Tprev = T;

		M[S.za] = T;
	}

/* ----- */
/* Apply */
/* ----- */

	TAffine	Tcache;
	int		zcache = -1, nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[(*zs)[i].i];

		if( R.itr < 0 )
			continue;

		if( X[R.itr * NX] == 999.0 )
			continue;

		TAffine	T( &X[R.itr * NX] );

		if( R.z != zcache ) {

			zcache = R.z;

			map<int,TAffine>::iterator	it = M.find( zcache );

			if( it != M.end() )
				Tcache = it->second;
			else
				Tcache.NUSetOne();
		}

		T = Tcache * T;
		T.CopyOut( &X[R.itr * NX] );
	}
}

/* --------------------------------------------------------------- */
/* AffineFromFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::AffineFromFile( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Apply previous final results at highest level, once only

	SetUniteLayer( LHS, RHS, sc );

// Standard starting point

	SetPointPairs( LHS, RHS, sc );

// Load the Affines A

	vector<double>	A;
	int				z0, nz;

	LoadAffTable( A, z0, nz, nTr );
	UntwistScaffold( A, z0, nz );

// Relatively weighted: A(pi) = A(pj)

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = A(p2)
		if( A[vRgn[C.r2].itr * NX] != 999.0 ) {

			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * scaf_strength / sc,
					y1 = C.p1.y * scaf_strength / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x * scaf_strength / sc;
			y2 = g2.y * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		if( A[vRgn[C.r1].itr * NX] != 999.0 ) {

			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * scaf_strength / sc,
					y1 = C.p2.y * scaf_strength / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x * scaf_strength / sc;
			y2 = g2.y * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-FrmFil", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );

//	DevFromPrior( A, X );

	fflush( stdout );
}

/* --------------------------------------------------------------- */
/* AffineFromFile2 ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::AffineFromFile2( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Apply previous final results at highest level, once only

	SetUniteLayer( LHS, RHS, sc );

// Standard starting point

	SetPointPairs( LHS, RHS, sc );

// Get the Affines A

	vector<double>	A;
	int				z0 = (*zs).begin()->z,
					nz = (*zs).rbegin()->z - z0 + 1;

	AffineFromFile( A, nTr );
	UntwistScaffold( A, z0, nz );

// Relatively weighted: A(pi) = A(pj)

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = A(p2)
		if( A[vRgn[C.r2].itr * NX] != 999.0 ) {

			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * scaf_strength / sc,
					y1 = C.p1.y * scaf_strength / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x * scaf_strength / sc;
			y2 = g2.y * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		if( A[vRgn[C.r1].itr * NX] != 999.0 ) {

			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * scaf_strength / sc,
					y1 = C.p2.y * scaf_strength / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x * scaf_strength / sc;
			y2 = g2.y * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-FrmFil2", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );

//	DevFromPrior( A, X );

	fflush( stdout );
}

/* --------------------------------------------------------------- */
/* AffineFromTransWt --------------------------------------------- */
/* --------------------------------------------------------------- */

// Preferred way to get montages:
// - Relate affines to each other A1(p1) = A2(p2).
// - Add scaffold relations A(pi) = T(pj) at reduced strength.
//
void MAffine::AffineFromTransWt( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Standard starting point

	SetPointPairs( LHS, RHS, sc );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;

	M.SetModelParams( gW, gH, -1, -1, -1, -1,
		nproc, -1, NULL, NULL, zs );
	M.SolveSystem( T, nTr );

// Relatively weighted A(pi) = T(pj)

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * scaf_strength / sc,
					y1 = C.p1.y * scaf_strength / sc,
					x2 = (C.p2.x + T[k  ]) * scaf_strength / sc,
					y2 = (C.p2.y + T[k+1]) * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * scaf_strength / sc,
					y1 = C.p2.y * scaf_strength / sc,
					x2 = (C.p1.x + T[k  ]) * scaf_strength / sc,
					y2 = (C.p1.y + T[k+1]) * scaf_strength / sc;

			double	v[3]	= {  x1,  y1, scaf_strength };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	WriteSolveRead( X, LHS, RHS, "A-WtT", nproc, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );

	DevFromTrans( T, X );
}

/* --------------------------------------------------------------- */
/* Fill_myc ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static vector<vector<int> >	myc;

void MAffine::Fill_myc( const vector<double> &X )
{
	myc.resize( vRgn.size() );

	int	nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		int	itr;

		if( (itr = vRgn[C.r1].itr) < 0 )
			continue;

		if( X[itr*NX] == 999.0 )
			continue;

		if( (itr = vRgn[C.r2].itr) < 0 )
			continue;

		if( X[itr*NX] == 999.0 )
			continue;

		myc[C.r1].push_back( i );
		myc[C.r2].push_back( i );
	}
}

/* --------------------------------------------------------------- */
/* AFromIDB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::AFromIDB( vector<double> &X, int nTr )
{
	int	nvars = nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	X.resize( nvars );

	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		const Til2Img*	t2i;
		IDBT2ICacheNGet1( t2i, idb, R.z, R.id );
		t2i->T.CopyOut( &X[R.itr * NX] );
	}

	Fill_myc( X );
}

/* --------------------------------------------------------------- */
/* AFromTbl ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::AFromTbl( vector<double> &X, int nTr )
{
	int	nvars = nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	X.resize( nvars );

	map<MZIDR,TAffine>	M;
	set<int>			Z;

	LoadTAffineTbl_AllZ( M, Z, "../montage/TAffineTable.txt" );

	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		map<MZIDR,TAffine>::iterator	it;

		it = M.find( MZIDR( R.z, R.id, R.rgn ) );

		if( it != M.end() )
			memcpy( &X[R.itr*6], it->second.t, 6*sizeof(double) );
		else
			R.itr = -1;
	}

	Fill_myc( X );
}

/* --------------------------------------------------------------- */
/* AFromScf ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::AFromScf( vector<double> &X, int nTr )
{
	int	nvars = nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	int	z0, nz;

	LoadAffTable( X, z0, nz, nTr );
	UntwistScaffold( X, z0, nz );

	Fill_myc( X );
}

/* --------------------------------------------------------------- */
/* OnePass ------------------------------------------------------- */
/* --------------------------------------------------------------- */

#if 0
void MAffine::OnePass(
	vector<double>	&Xout,
	vector<double>	&Xin,
	vector<double>	&S,
	int				nTr,
	double			w )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

sc = 1;
scaf_strength = 0.1;

	Xout.resize( nvars );

/* ---------------- */
/* Solve rgn by rgn */
/* ---------------- */

	int	nr = vRgn.size();

	vector<double>	x( 6 );
	int				i1[3] = { 0, 1, 2 },
					i2[3] = { 3, 4, 5 };

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		vector<double>	RHS( 6, 0.0 );
		vector<LHSCol>	LHS( 6 );
		int				nc = myc[i].size();

		for( int j = 0; j < nc; ++j ) {

			const Constraint&	C = vAllC[myc[i][j]];
			Point				A, B, As, Bs;

			if( C.r1 == i ) {
				As = A = C.p1;
				Bs = B = C.p2;
				L2GPoint( B, Xin, vRgn[C.r2].itr );
				L2GPoint( Bs, S, vRgn[C.r2].itr );
			}
			else {
				As = A = C.p2;
				Bs = B = C.p1;
				L2GPoint( B, Xin, vRgn[C.r1].itr );
				L2GPoint( Bs, S, vRgn[C.r1].itr );
			}

			double	v[3] = { A.x/sc, A.y/sc, 1.0 },
					vs[3] = { scaf_strength*A.x/sc, scaf_strength*A.y/sc, scaf_strength };

			AddConstraint( LHS, RHS, 3, i1, v, B.x/sc );
			AddConstraint( LHS, RHS, 3, i2, v, B.y/sc );

//			AddConstraint( LHS, RHS, 3, i1, vs, scaf_strength*Bs.x/sc );
//			AddConstraint( LHS, RHS, 3, i2, vs, scaf_strength*Bs.y/sc );
		}

		WriteSolveRead( x, LHS, RHS, "A-RLX", 1, false );
		memcpy( &Xout[R.itr*NX], &x[0], NX*sizeof(double) );
	}

	RescaleAll( Xout, sc );

// like w = 0.75

	for( int i = 0; i < nvars; ++i ) {

		Xout[i] = (1-w) * Xin[i] + w * Xout[i];
	}
}
#endif

#if 1

class Cperr {
public:
	double	e;
	int		j;
public:
	Cperr( double e, int j ) : e(e), j(j) {};
	bool operator < ( const Cperr &rhs ) const
		{return e < rhs.e;};
};

void MAffine::OnePass(
	vector<double>	&Xout,
	vector<double>	&Xin,
	vector<double>	&S,
	int				nTr,
	double			w )
{
	int	nvars = nTr * NX;

	Xout.resize( nvars );

/* ---------------- */
/* Solve rgn by rgn */
/* ---------------- */

	int	nr = vRgn.size();

	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		int	nc = myc[i].size();

		if( nc < 3 )
			continue;

		double	*RHS = &Xout[R.itr*NX];
		double	LHS[6*6];
		TAffine	Ta( &Xin[R.itr * NX] );

		memset( RHS, 0, 6   * sizeof(double) );
		memset( LHS, 0, 6*6 * sizeof(double) );

		for( int j = 0; j < nc; ++j ) {

			const Constraint&	C = vAllC[myc[i][j]];
			Point				A, B;

			// like w = 0.9 (same layer), 0.9 (down)

			if( C.r1 == i ) {
				TAffine Tb( &Xin[vRgn[C.r2].itr * NX] );
				Tb.Transform( B = C.p2 );
				Ta.Transform( A = C.p1 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p1;
			}
			else {
				TAffine Tb( &Xin[vRgn[C.r1].itr * NX] );
				Tb.Transform( B = C.p1 );
				Ta.Transform( A = C.p2 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p2;
			}

			double	v[3] = { A.x, A.y, 1.0 };

			AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
			AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
		}

		Solve_Quick( LHS, RHS, 6 );
//----------------------------------------------------------
// Exper to sort errors, kick out top n, then resolve
#if 0
		TAffine			Tx( &Xout[R.itr * NX] );
		vector<Cperr>	ve;

		// collect errors

		for( int j = 0; j < nc; ++j ) {

			const Constraint&	C = vAllC[myc[i][j]];
			Point				A, B;

			if( C.r1 == i ) {
				TAffine Tb( &Xin[vRgn[C.r2].itr * NX] );
				Tb.Transform( B = C.p2 );
				Tx.Transform( A = C.p1 );
			}
			else {
				TAffine Tb( &Xin[vRgn[C.r1].itr * NX] );
				Tb.Transform( B = C.p1 );
				Tx.Transform( A = C.p2 );
			}

			ve.push_back( Cperr( A.DistSqr( B ), j ) );
		}

		// sort by error
		sort( ve.begin(), ve.end() );

	{
		int	ne = nc - 1;	// omit top n

		memset( RHS, 0, 6   * sizeof(double) );
		memset( LHS, 0, 6*6 * sizeof(double) );

		for( int j = 0; j < ne; ++j ) {

			const Constraint&	C = vAllC[myc[i][ve[j].j]];
			Point				A, B;

			// like w = 0.9 (same layer), 0.9 (down)

			if( C.r1 == i ) {
				TAffine Tb( &Xin[vRgn[C.r2].itr * NX] );
				Tb.Transform( B = C.p2 );
				Ta.Transform( A = C.p1 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p1;
			}
			else {
				TAffine Tb( &Xin[vRgn[C.r1].itr * NX] );
				Tb.Transform( B = C.p1 );
				Ta.Transform( A = C.p2 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p2;
			}

			double	v[3] = { A.x, A.y, 1.0 };

			AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
			AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
		}

		Solve_Quick( LHS, RHS, 6 );
	}
#endif
//----------------------------------------------------------
	}
}
#endif

/* --------------------------------------------------------------- */
/* KeepSLOnly ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Cull downs from given myc list.
//
static void KeepSLOnly( vector<int> &v, int nc )
{
	int id = 0;

	for( int is = 0; is < nc; ++is ) {

		const Constraint&	C = vAllC[v[is]];

		if( C.r1 == C.r2 )
			v[id++] = v[is];
	}

	v.resize( id );
}

/* --------------------------------------------------------------- */
/* OnePassTH ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CThrdat {
public:
	int			r0,
				rlim;
	vector<int>	Rslo;
	vector<int>	Rkil;
};

static vector<CThrdat>	vthr;
static vector<double>	*Xout;
static vector<double>	*Xin;
static double			w;
static int				gpass;


static void* _OnePass_AFromA( void* ithr )
{
	CThrdat	&me = vthr[(long)ithr];

	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

	for( int i = me.r0; i < me.rlim; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		int	nc = myc[i].size();

		if( nc < 3 )
			continue;

		double	*RHS = &(*Xout)[R.itr * 6];
		double	LHS[6*6];
		TAffine	Ta( &(*Xin)[R.itr * 6] );
		TAffine	Tb;
		int		lastb = -1;	// cache Tb

		memset( RHS, 0, 6   * sizeof(double) );
		memset( LHS, 0, 6*6 * sizeof(double) );

		for( int j = 0; j < nc; ++j ) {

			const Constraint&	C = vAllC[myc[i][j]];
			Point				A, B;

			// Mixing old and new solutions is related to
			// "successive over relaxation" methods in other
			// iterative solution schemes. Experimentally,
			// I like w = 0.9 (same layer), 0.9 (down).

			if( C.r1 == i ) {

				int	bitr = vRgn[C.r2].itr;

				if( bitr < 0 )
					continue;

				if( C.r2 != lastb ) {
					Tb.CopyIn( &(*Xin)[bitr * 6] );
					lastb = C.r2;
				}
				Tb.Transform( B = C.p2 );
				Ta.Transform( A = C.p1 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p1;
			}
			else {

				int	bitr = vRgn[C.r1].itr;

				if( bitr < 0 )
					continue;

				if( C.r1 != lastb ) {
					Tb.CopyIn( &(*Xin)[bitr * 6] );
					lastb = C.r1;
				}
				Tb.Transform( B = C.p1 );
				Ta.Transform( A = C.p2 );
				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p2;
			}

			double	v[3] = { A.x, A.y, 1.0 };

			AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
			AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
		}

		Solve_Quick( LHS, RHS, 6 );
	}

	return NULL;
}


static void AFromA_SLOnly( double *RHS, int i, int ithr )
{
	const RGN&	R = vRgn[i];

	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

	int	nc = myc[i].size();

	double	LHS[6*6];
	TAffine	Ta( &(*Xin)[R.itr * 6] );
	TAffine	Tb;
	int		lastb	= -1,	// cache Tb
			nSLc	= 0;

	memset( RHS, 0, 6   * sizeof(double) );
	memset( LHS, 0, 6*6 * sizeof(double) );

	for( int j = 0; j < nc; ++j ) {

		const Constraint&	C = vAllC[myc[i][j]];
		Point				A, B;

		// Mixing old and new solutions is related to
		// "successive over relaxation" methods in other
		// iterative solution schemes. Experimentally,
		// I like w = 0.9 (same layer), 0.9 (down).

		if( C.r1 == i ) {

			if( vRgn[C.r2].z != R.z )
				continue;

			int	bitr = vRgn[C.r2].itr;

			if( bitr < 0 )
				continue;

			if( C.r2 != lastb ) {
				Tb.CopyIn( &(*Xin)[bitr * 6] );
				lastb = C.r2;
			}
			Tb.Transform( B = C.p2 );
			Ta.Transform( A = C.p1 );
			B.x = w * B.x + (1 - w) * A.x;
			B.y = w * B.y + (1 - w) * A.y;
			A = C.p1;
		}
		else {

			if( vRgn[C.r1].z != R.z )
				continue;

			int	bitr = vRgn[C.r1].itr;

			if( bitr < 0 )
				continue;

			if( C.r1 != lastb ) {
				Tb.CopyIn( &(*Xin)[bitr * 6] );
				lastb = C.r1;
			}
			Tb.Transform( B = C.p1 );
			Ta.Transform( A = C.p2 );
			B.x = w * B.x + (1 - w) * A.x;
			B.y = w * B.y + (1 - w) * A.y;
			A = C.p2;
		}

		++nSLc;

		double	v[3] = { A.x, A.y, 1.0 };

		AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
		AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
	}

	if( nSLc < 3 ||
		!Solve_Quick( LHS, RHS, 6 ) ||
		X_AS_AFF( RHS, 0 ).Squareness() > SQRTOL ) {

		vthr[ithr].Rkil.push_back( i );
	}
	else
		vthr[ithr].Rslo.push_back( i );
}


static void* _OnePass_AFromA_stk( void* ithr )
{
	CThrdat	&me = vthr[(int)(long)ithr];

	int	i1[3] = { 0, 1, 2 },
		i2[3] = { 3, 4, 5 };

	for( int i = me.r0; i < me.rlim; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.itr < 0 )
			continue;

		int	nc = myc[i].size();

		if( nc < 3 )
			continue;

		double		*RHS = &(*Xout)[R.itr * 6];
		double		LHS[6*6];
		TAffine*	Ta = &X_AS_AFF( *Xin, R.itr );
		TAffine*	Tb;
		int			lastb = -1;	// cache Tb

		memset( RHS, 0, 6   * sizeof(double) );
		memset( LHS, 0, 6*6 * sizeof(double) );

		for( int j = 0; j < nc; ++j ) {

			const Constraint&	C = vAllC[myc[i][j]];
			Point				A, B;

			// Mixing old and new solutions is related to
			// "successive over relaxation" methods in other
			// iterative solution schemes. Experimentally,
			// I like w = 0.9 (same layer), 0.9 (down).

			if( C.r1 == i ) {

				int	bitr = vRgn[C.r2].itr;

				if( bitr < 0 )
					continue;

				if( C.r2 != lastb ) {
					Tb = &X_AS_AFF( *Xin, bitr );
					lastb = C.r2;
				}
				Tb->Transform( B = C.p2 );
				Ta->Transform( A = C.p1 );

				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p1;
			}
			else {

				int	bitr = vRgn[C.r1].itr;

				if( bitr < 0 )
					continue;

				if( C.r1 != lastb ) {
					Tb = &X_AS_AFF( *Xin, bitr );
					lastb = C.r1;
				}
				Tb->Transform( B = C.p1 );
				Ta->Transform( A = C.p2 );

				B.x = w * B.x + (1 - w) * A.x;
				B.y = w * B.y + (1 - w) * A.y;
				A = C.p2;
			}

			double	v[3] = { A.x, A.y, 1.0 };

			AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
			AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
		}

		if( gpass < EDITDELAY ) {
			Solve_Quick( LHS, RHS, 6 );
			continue;
		}

		if( !Solve_Quick( LHS, RHS, 6 ) ||
			X_AS_AFF( RHS, 0 ).Squareness() > SQRTOL ) {

			AFromA_SLOnly( RHS, i, (int)(long)ithr );
		}
	}

	return NULL;
}


// Remove and remark bad RGNs and Constraints.
//
static void Remark()
{
	int	ncut = 0, nkil = 0;
	int	nt = vthr.size();

	// for each thread's lists...
	for( int it = 0; it < nt; ++it ) {

		vector<int>&	vslo = vthr[it].Rslo;
		vector<int>&	vkil = vthr[it].Rkil;
		int				nr = vslo.size();

		ncut += nr;

		// for the SLOnlys...
		for( int ir = 0; ir < nr; ++ir ) {

			int				islo = vslo[ir];
			vector<int>&	vmine = myc[islo];
			int				nc = vmine.size();

			// for each constraint...
			for( int ic = 0; ic < nc; ++ic ) {

				Constraint&	C = vAllC[vmine[ic]];

				if( vRgn[C.r1].z == vRgn[C.r2].z )
					continue;

				// mark constraint
				C.used = false;

				// remove it from other myc list
				int	iother = (C.r1 == islo ? C.r2 : C.r1);
				vector<int>&	vother = myc[iother];
				int				mc = vother.size();

				for( int jc = 0; jc < mc; ++jc ) {

					if( vother[jc] == vmine[ic] ) {

						vother.erase( vother.begin()+jc );

						// mark other for kill?
						if( vother.size() < 3 )
							vkil.push_back( iother );

						break;
					}
				}

				// and remove it from my list
				vmine.erase( vmine.begin()+ic );
				--ic;
				--nc;
			}

			// kill me?
			if( nc < 3 ) {
				vkil.push_back( islo );
				--ncut;
			}
		}

		// for the kills...
		for( int ir = 0; ir < vkil.size(); ++ir ) {

			int	ikil = vkil[ir];

			// already dead?
			if( vRgn[ikil].itr < 0 )
				continue;

			vector<int>&	vmine = myc[ikil];
			int				nc = vmine.size();

			// for each constraint...
			for( int ic = 0; ic < nc; ++ic ) {

				Constraint&	C = vAllC[vmine[ic]];

				// mark constraint
				C.used = false;

				// remove it from other myc list
				int	iother = (C.r1 == ikil ? C.r2 : C.r1);
				vector<int>&	vother = myc[iother];
				int				mc = vother.size();

				for( int jc = 0; jc < mc; ++jc ) {

					if( vother[jc] == vmine[ic] ) {

						vother.erase( vother.begin()+jc );

						// mark other for kill?
						if( vother.size() < 3 )
							vkil.push_back( iother );

						break;
					}
				}
			}

			// kill my list
			vmine.clear();

			// mark me dead
			vRgn[ikil].itr = -1;
			++nkil;
		}
	}

	if( ncut || nkil )
		printf( "Tiles [cut, killed] = [%d, %d].\n", ncut, nkil );
}


static void OnePassTH(
	EZThreadproc	passproc,
	vector<double>	&shXout,
	vector<double>	&shXin,
	int				nXout,
	int				nthr,
	double			shw )
{
	int	nr = vRgn.size(),	// regions total
		nb;					// regions per thread

	if( nthr > nr )
		nthr = nr;

	Xout	= &shXout;
	Xin		= &shXin;
	w		= shw;
	nb		= nr / nthr;

	vthr.resize( nthr );
	Xout->resize( nXout );

	vthr[0].r0		= 0;
	vthr[0].rlim	= nb;

	for( int i = 1; i < nthr; ++i ) {
		CThrdat	&C = vthr[i];
		C.r0	= vthr[i-1].rlim;
		C.rlim	= (i == nthr-1 ? nr : C.r0 + nb);
	}

	if( !EZThreads( passproc, nthr, 1, "passproc" ) )
		exit( 42 );

	Remark();
	vthr.clear();
}

/* --------------------------------------------------------------- */
/* SolveSystemStandard ------------------------------------------- */
/* --------------------------------------------------------------- */

// Build and solve system of linear equations.
//
// Note:
// All translational variables are scaled down by 'scale' so they
// are sized similarly to the sine/cosine variables. This is only
// to stabilize solver algorithm. We undo the scaling on exit.
//
void MAffine::SolveSystemStandard( vector<double> &X, int nTr )
{
	double	scale	= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

/* ------------------ */
/* Get rough solution */
/* ------------------ */

	SetPointPairs( LHS, RHS, scale );

	if( unite_layer < 0 )
		SetIdentityTForm( LHS, RHS, nTr / 2 );
	else
		SetUniteLayer( LHS, RHS, scale );

	SolveWithSquareness( X, LHS, RHS, nTr );

/* ----------------------------------------- */
/* Use solution to add torsional constraints */
/* ----------------------------------------- */

	SolveWithUnitMag( X, LHS, RHS, nTr );

/* --------------------------- */
/* Rescale translational terms */
/* --------------------------- */

	RescaleAll( X, scale );
}

/* --------------------------------------------------------------- */
/* Test ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::Test()
{
	map<MZIDR,TAffine>	M;
	set<int>			Z;

//	LoadTAffineTbl_AllZ( M, Z, "../cross_wkspc/HiRes_TF.txt" );
	LoadTAffineTbl_AllZ( M, Z, "../s_iter200_Scaf_A_2400/TAffineTable.txt" );

	map<MZIDR,TAffine>::iterator	it;
	FILE	*f = FileOpenOrDie( "Angles_2400.txt", "w" );

	int	ntot = M.size();
	int	b5 = 0, g5 = 0;

	for( it = M.begin(); it != M.end(); ++it ) {

		const MZIDR&	W = it->first;
		const TAffine&	T = it->second;
		double			c = T.Squareness();

		if( c < .174 )
			++b5;
		else {
			++g5;
			fprintf( f, "%d\t%d\t%f\n", W.z, W.id, c );
		}
	}

	fclose( f );

	printf( "tot=%d b5=%d g5=%d\n", ntot, b5, g5 );
}

static void XWrite( const vector<double> &X )
{
	FILE	*f = FileOpenOrDie( "X_A.dat", "wb" );
	fwrite( &X[0], sizeof(double), X.size(), f );
	fclose( f );
}

static void XRead( vector<double> &X )
{
	FILE	*f = FileOpenOrDie( "X_A.dat", "rb" );
	fread( &X[0], sizeof(double), X.size(), f );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SolveSystem( vector<double> &X, int nTr )
{
	clock_t	t0 = StartTiming();

#if 0

// This older method works fine for small montages, especially
// if small in the y-dimension, but there will be some banana
// curvature visible. The newer scaffold method remains more
// square at any size, so is preferred.

	SolveSystemStandard( X, nTr );

#else

	if( priorafftbl ) {

		fflush( stdout );

#if 0	// stack with SLU

//		AffineFromFile( X, nTr );
		AffineFromFile2( X, nTr );

#elif 0	// stack iterative

		AFromScf( X, nTr );
		vector<double>	S = X;
		clock_t			t2 = StartTiming();
		for( int i = 0; i < 150; ++i ) {	// 25
			vector<double> Xin = X;
			OnePass( X, Xin, S, nTr, 0.9 );
			printf( "Done pass %d\n", i + 1 ); fflush( stdout );
		}
		PrintMagnitude( X );
		StopTiming( stdout, "Iters", t2 );
		myc.clear();

#else	// threaded stack iterative

		AFromScf( X, nTr );
		clock_t			t2 = StartTiming();
		for( gpass = 0; gpass < 4000; ++gpass ) {	// 25
			vector<double> Xin = X;
			OnePassTH( _OnePass_AFromA_stk, X, Xin, nTr*6, 16, 0.9 );
//			printf( "Done pass %d\n", gpass + 1 ); fflush( stdout );
		}
		PrintMagnitude( X );
		StopTiming( stdout, "Iters", t2 );
		myc.clear();

#endif
	}
	else {

#if 0	// montage with SLU

		AffineFromTransWt( X, nTr );

#elif 0	// montage iterative

		AFromIDB( X, nTr );
		vector<double>	S = X;
		clock_t			t2 = StartTiming();
		for( int i = 0; i < 10; ++i ) {	// 250
			vector<double> Xin = X;
			OnePass( X, Xin, S, nTr, 0.9 );
		}
		PrintMagnitude( X );
		StopTiming( stdout, "Iters", t2 );
		myc.clear();

#else	// threaded montage iterative

		vector<double>	Xin;
		clock_t			t2 = StartTiming();

		for( gpass = 0; gpass < 250; ++gpass ) {	// 250

			if( !gpass )
				AFromIDB( Xin, nTr );
			else
				XRead( Xin );

			OnePassTH( _OnePass_AFromA_stk, X, Xin, nTr*6, 8, 0.9 );
			XWrite( X );
		}
		PrintMagnitude( X );
		StopTiming( stdout, "Iters", t2 );
		myc.clear();

#endif
	}

#endif

	StopTiming( stdout, "SolveA", t0 );
}

/* --------------------------------------------------------------- */
/* UpdateScaffold ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::UpdateScaffold( vector<double> &X, int nTr )
{
	int	z0, nz;

	LoadAffTable( X, z0, nz, nTr );
	UntwistScaffold( X, z0, nz );

	double			xbnd, ybnd;
	vector<double>	lrbt;

	Bounds( xbnd, ybnd, X, lrbt, 0, NULL );
	WriteTransforms( X, false, NULL );

	WriteTrakEM( xbnd, ybnd, X, 0, 0, 0, 0 );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTransforms(
	const vector<double>	&X,
	int						bstrings,
	FILE					*FOUT )
{
	printf( "---- Write transforms ----\n" );

	FILE	*ft  = FileOpenOrDie( "TAffineTable.txt", "w" ),
			*fx  = FileOpenOrDie( "magnitude_outliers.txt", "w" );
	double	smin = 100.0,
			smax = 0.0,
			smag = 0.0;
	int		nr   = vRgn.size(), nTr = 0;

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[(*zs)[i].i];

		if( I.itr < 0 )
			continue;

		int	j = I.itr * NX;

		++nTr;

		fprintf( ft, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
		I.z, I.id, I.rgn,
		X[j  ], X[j+1], X[j+2],
		X[j+3], X[j+4], X[j+5] );

		if( FOUT ) {

			if( !bstrings ) {

				fprintf( FOUT, "TAFFINE %d.%d:%d %f %f %f %f %f %f\n",
				I.z, I.id, I.rgn,
				X[j  ], X[j+1], X[j+2],
				X[j+3], X[j+4], X[j+5] );
			}
			else {
				fprintf( FOUT, "TRANSFORM '%s::%d' %f %f %f %f %f %f\n",
				I.GetName(), I.rgn,
				X[j  ], X[j+1], X[j+2],
				X[j+3], X[j+4], X[j+5] );
			}
		}

		double	mag = Magnitude( X, I.itr );

		smag += mag;
		smin  = fmin( smin, mag );
		smax  = fmax( smax, mag );

		if( mag < 0.9 ) {
			fprintf( fx, "Low mag %f @ %d.%d:%d\n",
				mag, I.z, I.id, I.rgn );
		}
		else if( mag > 1.1 ) {
			fprintf( fx, "Hi  mag %f @ %d.%d:%d\n",
				mag, I.z, I.id, I.rgn );
		}
	}

	fclose( fx );
	fclose( ft );

	printf(
	"Average magnitude=%f, min=%f, max=%f, max/min=%f\n\n",
	smag/nTr, smin, smax, smax/smin );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTrakEM(
	double					xmax,
	double					ymax,
	const vector<double>	&X,
	double					trim,
	int						xml_type,
	int						xml_min,
	int						xml_max )
{
	FILE	*f = FileOpenOrDie( "MultLayAff.xml", "w" );

	int	oid = 3;

	fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

	TrakEM2WriteDTD( f );

	fprintf( f, "<trakem2>\n" );

	fprintf( f,
	"\t<project\n"
	"\t\tid=\"0\"\n"
	"\t\ttitle=\"Project\"\n"
	"\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
	"\t\tn_mipmap_threads=\"8\"\n"
	"\t/>\n" );

	fprintf( f,
	"\t<t2_layer_set\n"
	"\t\toid=\"%d\"\n"
	"\t\ttransform=\"matrix(1,0,0,1,0,0)\"\n"
	"\t\ttitle=\"Top level\"\n"
	"\t\tlayer_width=\"%.2f\"\n"
	"\t\tlayer_height=\"%.2f\"\n"
	"\t>\n",
	oid++, xmax, ymax );

	int	prev	= -1;	// will be previously written layer
	int	offset	= int(2 * trim + 0.5);
	int	nr		= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[(*zs)[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		// changed layer
		if( (*zs)[i].z != prev ) {

			if( prev != -1 )
				fprintf( f, "\t\t</t2_layer>\n" );

			fprintf( f,
			"\t\t<t2_layer\n"
			"\t\t\toid=\"%d\"\n"
			"\t\t\tthickness=\"0\"\n"
			"\t\t\tz=\"%d\"\n"
			"\t\t>\n",
			oid++, (*zs)[i].z );

			prev = (*zs)[i].z;
		}

		const char	*path;
		char		title[128];
		DisplayStrings( title, path, I );

		// fix origin : undo trimming
		int		j = I.itr * NX;
		double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
		double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

		fprintf( f,
		"\t\t\t<t2_patch\n"
		"\t\t\t\toid=\"%d\"\n"
		"\t\t\t\twidth=\"%d\"\n"
		"\t\t\t\theight=\"%d\"\n"
		"\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
		"\t\t\t\ttitle=\"%s\"\n"
		"\t\t\t\ttype=\"%d\"\n"
		"\t\t\t\tfile_path=\"%s\"\n"
		"\t\t\t\to_width=\"%d\"\n"
		"\t\t\t\to_height=\"%d\"\n",
		oid++, gW - offset, gH - offset,
		X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
		title, xml_type, path, gW - offset, gH - offset );

		if( xml_min < xml_max ) {

			fprintf( f,
			"\t\t\t\tmin=\"%d\"\n"
			"\t\t\t\tmax=\"%d\"\n"
			"\t\t\t/>\n",
			xml_min, xml_max );
		}
		else
			fprintf( f, "\t\t\t/>\n" );
	}

	if( nr > 0 )
		fprintf( f, "\t\t</t2_layer>\n" );

	fprintf( f, "\t</t2_layer_set>\n" );
	fprintf( f, "</trakem2>\n" );
	fclose( f );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteJython(
	const vector<double>	&X,
	double					trim,
	int						Ntr )
{
	FILE	*f = FileOpenOrDie( "JythonTransforms.txt", "w" );

	fprintf( f, "transforms = {\n" );

	int	nr = vRgn.size();

	for( int i = 0, itrf = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[(*zs)[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		++itrf;

		const char	*path;
		DisplayStrings( NULL, path, I );

		// fix origin : undo trimming
		int		j = I.itr * NX;
		double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
		double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

		fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f]%s\n",
			path, X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
			(itrf == Ntr ? "" : ",") );
	}

	fprintf( f, "}\n" );
	fclose( f );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::G2LPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	TAffine	I, T( &X[itr * NX] );
	I.InverseOf( T );
	I.Transform( p );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::L2GPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	TAffine	T( &X[itr * NX] );
	T.Transform( p );
}


void MAffine::L2GPoint(
	vector<Point>			&p,
	const vector<double>	&X,
	int						itr )
{
	TAffine	T( &X[itr * NX] );
	T.Transform( p );
}


