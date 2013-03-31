

#include	"lsq_Trans.h"
#include	"lsq_Affine.h"
#include	"lsq_Hmgphy.h"

#include	"TrakEM2_UTL.h"
#include	"PipeFiles.h"
#include	"File.h"
#include	"Maths.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SetPointPairs(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc,
	double			same_strength )
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

		// H1(p1) - H2(p2) = 0

		double	va[6] = { x1,  y1,  fz, -x2, -y2, -fz};
		int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
		int		i2[6] = {j+3, j+4, j+5, k+3, k+4, k+5};

		AddConstraint( LHS, RHS, 6, i1, va, 0.0 );
		AddConstraint( LHS, RHS, 6, i2, va, 0.0 );

		double	vb[4] = { x1,  y1, -x2, -y2};
		int		i3[4] = {j+6, j+7, k+6, k+7};

		AddConstraint( LHS, RHS, 4, i3, vb, 0.0 );
	}
}

/* --------------------------------------------------------------- */
/* SetIdentityTForm ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Explicitly set some TForm to Identity.
// @@@ Does it matter which one we use?
//
void MHmgphy::SetIdentityTForm(
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
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;

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
// solution output file gArgs.tfm_file.
//
void MHmgphy::SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc,
	int				unite_layer,
	const char		*tfm_file )
{
/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

	map<MZIDR,THmgphy>	M;

	LoadTHmgphyTbl_ThisZ( M, unite_layer, tfm_file );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

	double	stiff	= 10.0;

	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	R = vRgn[i];

		if( R.z != unite_layer || R.itr < 0 )
			continue;

		map<MZIDR,THmgphy>::iterator	it;

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
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[6] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[7] );		j++;
	}
}

/* --------------------------------------------------------------- */
/* SolveWithSquareness ------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SolveWithSquareness(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTr,
	double			square_strength )
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

	printf( "Solve with [transform squareness].\n" );
	WriteSolveRead( X, LHS, RHS, false );
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
void MHmgphy::SolveWithUnitMag(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTR,
	double			scale_strength )
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

	printf( "Solve with [unit magnitude].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	printf( "\t\t\t\t" );
	PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* RescaleAll ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::RescaleAll(
	vector<double>	&X,
	double			sc )
{
	THmgphy	D, U;
	int		nr = vRgn.size();

	U.NUSetScl( sc );
	D.NUSetScl( 1/sc );

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		THmgphy	T( &X[itr] );

		T = U * (T * D);

		T.CopyOut( &X[itr] );
	}
}

/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::RotateAll(
	vector<double>	&X,
	double			degcw )
{
	THmgphy	T, R;
	int		nr	= vRgn.size();

	R.SetCWRot( degcw, Point(0,0) );

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		THmgphy	t( &X[itr] );

		T = R * t;
		T.CopyOut( &X[itr] );
	}
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::NewOriginAll(
	vector<double>	&X,
	double			xorg,
	double			yorg )
{
	THmgphy	M( 1,0,-xorg, 0,1,-yorg, 0,0 );
	int		nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		itr *= NX;

		THmgphy	T( &X[itr] );

		T = M * T;

		T.CopyOut( &X[itr] );
	}
}

/* --------------------------------------------------------------- */
/* HmgphyEquHmgphy2 ---------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::HmgphyEquHmgphy2(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Affines A

	vector<double>	A;
	HmgphyEquHmgphy( A, nTr, gW, gH, same_strength, square_strength );

// SetPointPairs: H(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = A(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Set identity

	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr, square_strength );
	//SolveWithUnitMag( X, LHS, RHS, nTr, square_strength );

	printf( "Solve with [2nd hmgphied points].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyEquHmgphy ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::HmgphyEquHmgphy(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Affines A

	vector<double>	A;
	HmgphyEquAffine( A, nTr, gW, gH, same_strength, square_strength );

// SetPointPairs: H(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = A(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Set identity

	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr, square_strength );
	//SolveWithUnitMag( X, LHS, RHS, nTr, square_strength );

	printf( "Solve with [hmgphied points].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyEquAffine ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::HmgphyEquAffine(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the Affines A

	MAffine			M;
	vector<double>	A;
	M.SolveSystem( A, nTr, gW, gH, same_strength, square_strength, 0, -1, NULL );

// SetPointPairs: H(pi) = A(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = A(p2)
		{
			int		j  = vRgn[C.r1].itr * NX;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p2;

			M.L2GPoint( g2, A, vRgn[C.r2].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = A(p1)
		{
			int		j  = vRgn[C.r2].itr * NX;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2,
					y2;
			Point	g2 = C.p1;

			M.L2GPoint( g2, A, vRgn[C.r1].itr );
			x2 = g2.x / sc;
			y2 = g2.y / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Set identity

	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr, square_strength );
	//SolveWithUnitMag( X, LHS, RHS, nTr, square_strength );

	printf( "Solve with [affined points].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyEquTransWt ---------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::HmgphyEquTransWt(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Standard starting point

	SetPointPairs( LHS, RHS, sc, same_strength );
	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;
	M.SolveSystem( T, nTr, 0, 0, 0, 0, 0, -1, NULL );

// Relatively weighted H(pi) = T(pj)

	double	fz	= 0.01;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2 = (C.p2.x + T[k  ]) / sc,
					y2 = (C.p2.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2 = (C.p1.x + T[k  ]) / sc,
					y2 = (C.p1.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr, square_strength );
	//SolveWithUnitMag( X, LHS, RHS, nTr, square_strength );

	printf( "Solve with [fixed translation].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyEquTrans ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MHmgphy::HmgphyEquTrans(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			square_strength )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;
	M.SolveSystem( T, nTr, 0, 0, 0, 0, 0, -1, NULL );

// SetPointPairs: H(pi) = T(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// H(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2 = (C.p2.x + T[k  ]) / sc,
					y2 = (C.p2.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}

		// H(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2 = (C.p1.x + T[k  ]) / sc,
					y2 = (C.p1.y + T[k+1]) / sc;

			double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
			int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
					i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

			AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

			v[3] = -x1*y2;
			v[4] = -y1*y2;

			AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
		}
	}

// Set identity

	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr, square_strength );
	//SolveWithUnitMag( X, LHS, RHS, nTr, square_strength );

	printf( "Solve with [fixed translation].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* SolveSystemStandard ------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SolveSystemStandard(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength,
	double			scale_strength,
	int				unite_layer,
	const char		*tfm_file )
{
	double	scale	= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Hmg: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

/* ------------------ */
/* Get rough solution */
/* ------------------ */

	SetPointPairs( LHS, RHS, scale, same_strength );

	if( unite_layer < 0 )
		SetIdentityTForm( LHS, RHS, nTr / 2 );
	else
		SetUniteLayer( LHS, RHS, scale, unite_layer, tfm_file );

	SolveWithSquareness( X, LHS, RHS, nTr, square_strength );

/* ----------------------------------------- */
/* Use solution to add torsional constraints */
/* ----------------------------------------- */

	//if( gArgs.make_layer_square )
	//	SolveWithMontageSqr( X, LHS, RHS );

	SolveWithUnitMag( X, LHS, RHS, nTr, scale_strength );

/* --------------------------- */
/* Rescale translational terms */
/* --------------------------- */

	RescaleAll( X, scale );
}

/* --------------------------------------------------------------- */
/* WriteSideRatios ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to guage trapezoidism by reporting the ratio of
// each image's left vertical side over its right side.
//
// A B
// D C
//
void MHmgphy::WriteSideRatios(
	const vector<zsort>		&zs,
	const vector<double>	&X )
{
	FILE	*f	= FileOpenOrDie( "SideRatios.txt", "w" );
	int		nr	= vRgn.size();
	MeanStd	M[4];

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[zs[i].i];

		if( I.itr < 0 )
			continue;

		int	j = I.itr * NX;

		THmgphy		T( &X[j] );
		Point		A( 0,    0 ), B( 2200,    0 ),
					D( 0, 2200 ), C( 2200, 2200 );
		double		d;
		const char	*c, *n = FileNamePtr( I.GetName() );
		int			cam = 0;

		T.Transform( A );
		T.Transform( B );
		T.Transform( C );
		T.Transform( D );

		d = D.Dist( A ) / C.Dist( B );

		B.x -= A.x;
		B.y -= A.y;
		C.x -= A.x;
		C.y -= A.y;
		D.x -= A.x;
		D.y -= A.y;
		A.x  = 0;
		A.y  = 0;

		if( c = strstr( n, "_cam" ) )
			cam = atoi( c + 4 );

		fprintf( f,
		"%d\t%g\t%g\t%g"
		"\t\t%g\t%g\t\t%g\t%g\t\t%g\t%g\t\t%g\t%g\n",
		cam, d, X[j+6], X[j+7],
		A.x, A.y, B.x, B.y, C.x, C.y, D.x, D.y );

		M[cam].Element( d );
	}

	for( int i = 0; i < 4; ++i ) {

		double	ave, std;

		M[i].Stats( ave, std );

		printf( "{cam,L/R,std}: {%d,%g,%g}\n", i, ave, std );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SolveSystem(
	vector<double>	&X,
	int				nTr,
	int				gW,
	int				gH,
	double			same_strength,
	double			square_strength,
	double			scale_strength,
	int				unite_layer,
	const char		*tfm_file )
{
	//SolveSystemStandard( X, nTr, gW, gH,
	//	same_strength, square_strength,
	//	scale_strength, unite_layer, tfm_file );

	//HmgphyEquTrans( X, nTr, gW, gH, square_strength );

	//HmgphyEquTransWt( X, nTr, gW, gH,
	//	same_strength, square_strength );

	//HmgphyEquAffine( X, nTr, gW, gH,
	//	same_strength, square_strength );

	HmgphyEquHmgphy( X, nTr, gW, gH,
		same_strength, square_strength );

	//HmgphyEquHmgphy2( X, nTr, gW, gH,
	//	same_strength, square_strength );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::WriteTransforms(
	const vector<zsort>		&zs,
	const vector<double>	&X,
	int						bstrings,
	FILE					*FOUT )
{
	printf( "---- Write transforms ----\n" );

	FILE	*f   = FileOpenOrDie( "THmgphyTable.txt", "w" );
	double	smin = 100.0,
			smax = 0.0,
			smag = 0.0;
	int		nr   = vRgn.size(), nTr = 0;

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[zs[i].i];

		if( I.itr < 0 )
			continue;

		int	j = I.itr * NX;

		++nTr;

		fprintf( f, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
		I.z, I.id, I.rgn,
		X[j  ], X[j+1], X[j+2],
		X[j+3], X[j+4], X[j+5],
		X[j+6], X[j+7] );

		if( !bstrings ) {

			fprintf( FOUT, "THMGPHY %d.%d:%d %f %f %f %f %f %f %f %f\n",
			I.z, I.id, I.rgn,
			X[j  ], X[j+1], X[j+2],
			X[j+3], X[j+4], X[j+5],
			X[j+6], X[j+7] );
		}
		else {
			fprintf( FOUT, "THMGPHY '%s::%d' %f %f %f %f %f %f %f %f\n",
			I.GetName(), I.rgn,
			X[j  ], X[j+1], X[j+2],
			X[j+3], X[j+4], X[j+5],
			X[j+6], X[j+7] );
		}

		double	mag = Magnitude( X, I.itr );

		smag += mag;
		smin  = fmin( smin, mag );
		smax  = fmax( smax, mag );
	}

	fclose( f );

	printf(
	"Average magnitude=%f, min=%f, max=%f, max/min=%f\n\n",
	smag/nTr, smin, smax, smax/smin );

	WriteSideRatios( zs, X );
}

/* --------------------------------------------------------------- */
/* TopLeft ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void TopLeft(
	double			&top,
	double			&left,
	const THmgphy	&T,
	int				gW,
	int				gH,
	double			trim )
{
	vector<Point>	cnr( 4 );

	cnr[0] = Point(      trim,      trim );
	cnr[1] = Point( gW-1-trim,      trim );
	cnr[2] = Point( gW-1-trim, gH-1-trim );
	cnr[3] = Point(      trim, gH-1-trim );

	T.Transform( cnr );

	top  = BIGD;
	left = BIGD;

	for( int k = 0; k < 4; ++k ) {

		top  = fmin( top,  cnr[k].y );
		left = fmin( left, cnr[k].x );
	}
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::WriteTrakEM(
	double					xmax,
	double					ymax,
	const vector<zsort>		&zs,
	const vector<double>	&X,
	int						gW,
	int						gH,
	double					trim,
	int						xml_type,
	int						xml_min,
	int						xml_max )
{
	FILE	*f = FileOpenOrDie( "MultLayHmgphy.xml", "w" );

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

		const RGN&	I = vRgn[zs[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		// changed layer
		if( zs[i].z != prev ) {

			if( prev != -1 )
				fprintf( f, "\t\t</t2_layer>\n" );

			fprintf( f,
			"\t\t<t2_layer\n"
			"\t\t\toid=\"%d\"\n"
			"\t\t\tthickness=\"0\"\n"
			"\t\t\tz=\"%d\"\n"
			"\t\t>\n",
			oid++, zs[i].z );

			prev = zs[i].z;
		}

		// trim trailing quotes and '::'
		// s = filename only
		char		buf[2048];
		strcpy( buf, I.GetName() );
		char		*p = strtok( buf, " ':\n" );
		const char	*s1 = FileNamePtr( p ),
					*s2	= FileDotPtr( s1 );

		// fix origin : undo trimming
		int		j = I.itr * NX;
		THmgphy	T( &X[j] );
		double	x_orig;
		double	y_orig;

		TopLeft( y_orig, x_orig, T, gW, gH, trim );

		fprintf( f,
		"\t\t\t<t2_patch\n"
		"\t\t\t\toid=\"%d\"\n"
		"\t\t\t\twidth=\"%d\"\n"
		"\t\t\t\theight=\"%d\"\n"
		"\t\t\t\ttransform=\"matrix(1,0,0,1,%f,%f)\"\n"
		"\t\t\t\ttitle=\"%.*s\"\n"
		"\t\t\t\ttype=\"%d\"\n"
		"\t\t\t\tfile_path=\"%s\"\n"
		"\t\t\t\to_width=\"%d\"\n"
		"\t\t\t\to_height=\"%d\"\n",
		oid++, gW - offset, gH - offset,
		x_orig, y_orig,
		s2 - s1, s1, xml_type, p, gW - offset, gH - offset );

		if( xml_min < xml_max ) {

			fprintf( f,
			"\t\t\t\tmin=\"%d\"\n"
			"\t\t\t\tmax=\"%d\"\n"
			"\t\t\t>\n",
			xml_min, xml_max );
		}
		else
			fprintf( f, "\t\t\t>\n" );

		fprintf( f,
		"\t\t\t<ict_transform"
		" class=\"mpicbg.trakem2.transform.HomographyModel2D\""
		" data=\"%g %g %g %g %g %g %g %g 1\"/>\n"
		"\t\t\t</t2_patch>\n",
		X[j  ], X[j+1], X[j+2], X[j+3],
		X[j+4], X[j+5], X[j+6], X[j+7] );
	}

	if( nr > 0 )
		fprintf( f, "\t\t</t2_layer>\n" );

	fprintf( f, "\t</t2_layer_set>\n" );
	fprintf( f, "</trakem2>\n" );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::WriteJython(
	const vector<zsort>		&zs,
	const vector<double>	&X,
	int						gW,
	int						gH,
	double					trim,
	int						Ntr )
{
	FILE	*f = FileOpenOrDie( "JythonTransforms.txt", "w" );

	fprintf( f, "transforms = {\n" );

	int	nr = vRgn.size();

	for( int i = 0, itrf = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[zs[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		++itrf;

		// trim trailing quotes and '::'
		char	buf[2048];
		strcpy( buf, I.GetName() );
		char	*p = strtok( buf, " ':\n" );

		// fix origin : undo trimming
		int		j = I.itr * NX;
		THmgphy	T( &X[j] );
		double	x_orig;
		double	y_orig;

		TopLeft( y_orig, x_orig, T, gW, gH, trim );

		fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f, %f, %f]%s\n",
			p,
			X[j+0], X[j+3], X[j+6],
			X[j+1], X[j+4], X[j+7],
			x_orig, y_orig,
			(itrf == Ntr ? "" : ",") );
	}

	fprintf( f, "}\n" );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MHmgphy::G2LPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	THmgphy	I, T( &X[itr * NX] );
	I.InverseOf( T );
	I.Transform( p );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MHmgphy::L2GPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	THmgphy	T( &X[itr * NX] );
	T.Transform( p );
}


void MHmgphy::L2GPoint(
	vector<Point>			&p,
	const vector<double>	&X,
	int						itr )
{
	THmgphy	T( &X[itr * NX] );
	T.Transform( p );
}


