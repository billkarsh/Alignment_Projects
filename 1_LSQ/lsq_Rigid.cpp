

#include	"lsq_Rigid.h"

#include	"PipeFiles.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::SetPointPairs(
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

		// T1(p1) - T2(p2) = 0
		//
		// {X0,X1,X2,X3} = {c,s,x,y}

		double	v[6]  = { x1, -y1,  fz, -x2,  y2, -fz};
		int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
		int		i2[6] = {j+1,   j, j+3, k+1,   k, k+3};

		AddConstraint( LHS, RHS, 6, i1, v, 0.0 );

		v[1] =  y1;
		v[4] = -y2;

		AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
	}
}

/* --------------------------------------------------------------- */
/* SetIdentityTForm ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Explicitly set some TForm to Identity.
// @@@ Does it matter which one we use?
//
void MRigid::SetIdentityTForm(
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
void MRigid::SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc,
	int				unite_layer,
	const char		*tfm_file )
{
/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

	map<MZIDR,TAffine>	M;

	LoadTAffineTbl_ThisZ( M, unite_layer, tfm_file );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

	double	stiff	= 10.0;

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
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[3] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[2] / sc );	j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[5] / sc );	j++;
	}
}

/* --------------------------------------------------------------- */
/* SolveWithSquareness ------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::SolveWithSquareness(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTr,
	double			square_strength )
{
/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

	printf( "Solve with [rigid only].\n" );
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
void MRigid::SolveWithUnitMag(
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
		double	s = X[j+1];
		double	m = sqrt( c*c + s*s );

		// c*x/m + s*y/m = 1

		double	V[2] = {c * stiff, s * stiff};
		int		I[2] = {j, j+1};

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

void MRigid::RescaleAll(
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
		X[itr+3] *= sc;
	}
}

/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::RotateAll(
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

		TAffine	t(
			X[itr  ], -X[itr+1], X[itr+2],
			X[itr+1],  X[itr  ], X[itr+3] );

		T = R * t;

		X[itr]		= T.t[0];
		X[itr+1]	= T.t[3];
		X[itr+2]	= T.t[2];
		X[itr+3]	= T.t[5];
	}
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::NewOriginAll(
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
		X[itr+3] -= yorg;
	}
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::WriteTransforms(
	const vector<zsort>		&zs,
	const vector<double>	&X,
	int						bstrings,
	FILE					*FOUT )
{
	printf( "---- Write transforms ----\n" );

	FILE	*f   = FileOpenOrDie( "TAffineTable.txt", "w" );
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

		fprintf( f, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
		I.z, I.id, I.rgn,
		X[j  ], -X[j+1], X[j+2],
		X[j+1],  X[j  ], X[j+3] );

		if( !bstrings ) {

			fprintf( FOUT, "TAFFINE %d.%d:%d %f %f %f %f %f %f\n",
			I.z, I.id, I.rgn,
			X[j  ], -X[j+1], X[j+2],
			X[j+1],  X[j  ], X[j+3] );
		}
		else {
			fprintf( FOUT, "TRANSFORM '%s::%d' %f %f %f %f %f %f\n",
			I.GetName(), I.rgn,
			X[j  ], -X[j+1], X[j+2],
			X[j+1],  X[j  ], X[j+3] );
		}

		double	mag = sqrt( X[j]*X[j] - X[j+1]*X[j+1] );

		smag += mag;
		smin  = fmin( smin, mag );
		smax  = fmax( smax, mag );
	}

	fclose( f );

	printf(
	"Average magnitude=%f, min=%f, max=%f, max/min=%f\n\n",
	smag/nTr, smin, smax, smax/smin );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MRigid::L2GPoint(
	Point			&p,
	vector<double>	&X,
	int				itr )
{
	int		j = itr * NX;
	TAffine	T( X[j], -X[j+1], X[j+2], X[j+1], X[j], X[j+3] );
	T.Transform( p );
}


void MRigid::L2GPoint(
	vector<Point>	&p,
	vector<double>	&X,
	int				itr )
{
	int		j = itr * NX;
	TAffine	T( X[j], -X[j+1], X[j+2], X[j+1], X[j], X[j+3] );
	T.Transform( p );
}


