

#include	"lsq_MTrans.h"
#include	"lsq_MAffine.h"

#include	"TrakEM2_UTL.h"
#include	"PipeFiles.h"
#include	"File.h"

#include	<math.h>


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
void MAffine::SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
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
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[1] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[2] / sc );	j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[3] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[4] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[5] / sc );	j++;
	}
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

	printf( "Solve [affines with transform squareness].\n" );
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

	printf( "Solve [affines with unit magnitude].\n" );
	WriteSolveRead( X, LHS, RHS, false );
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
/* DeviantAffines ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Experiment to see how much translation terms of each affine
// have moved from the trans-only starting values. We just list
// all tiles with dev > XXX, but do nothing with that for now.
//
void MAffine::DeviantAffines(
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

		if( (dx = sqrt( dx*dx + dy*dy  )) > 200 )
			printf( "Dev: %d/%d dr= %d\n", I.z, I.id, int(dx) );
	}
}

/* --------------------------------------------------------------- */
/* AffineFromTransWt --------------------------------------------- */
/* --------------------------------------------------------------- */

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
	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;

	M.SetModelParams( gW, gH,
		same_strength,
		square_strength,
		scale_strength,
		unite_layer, tfm_file );

	M.SolveSystem( T, nTr );

// Relatively weighted A(pi) = T(pj)

	double	fz	= 0.01;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2 = (C.p2.x + T[k  ]) * fz / sc,
					y2 = (C.p2.y + T[k+1]) * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2 = (C.p1.x + T[k  ]) * fz / sc,
					y2 = (C.p1.y + T[k+1]) * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	printf( "Solve [affines from translations].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );

	DeviantAffines( T, X );
}

/* --------------------------------------------------------------- */
/* AffineFromTrans ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::AffineFromTrans( vector<double> &X, int nTr )
{
	double	sc		= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Aff: %d unknowns; %d constraints.\n",
		nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

// Get the pure translations T

	MTrans			M;
	vector<double>	T;

	M.SetModelParams( gW, gH,
		same_strength,
		square_strength,
		scale_strength,
		unite_layer, tfm_file );

	M.SolveSystem( T, nTr );

// SetPointPairs: A(pi) = T(pj)

	double	fz	= 1.0;
	int		nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		// A(p1) = T(p2)
		{
			int		j  = vRgn[C.r1].itr * NX,
					k  = vRgn[C.r2].itr * 2;
			double	x1 = C.p1.x * fz / sc,
					y1 = C.p1.y * fz / sc,
					x2 = (C.p2.x + T[k  ]) * fz / sc,
					y2 = (C.p2.y + T[k+1]) * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}

		// A(p2) = T(p1)
		{
			int		j  = vRgn[C.r2].itr * NX,
					k  = vRgn[C.r1].itr * 2;
			double	x1 = C.p2.x * fz / sc,
					y1 = C.p2.y * fz / sc,
					x2 = (C.p1.x + T[k  ]) * fz / sc,
					y2 = (C.p1.y + T[k+1]) * fz / sc;

			double	v[3]	= {  x1,  y1,  fz };
			int		i1[3]	= {   j, j+1, j+2 },
					i2[3]	= { j+3, j+4, j+5 };

			AddConstraint( LHS, RHS, 3, i1, v, x2 );
			AddConstraint( LHS, RHS, 3, i2, v, y2 );
		}
	}

// Set identity

	SetIdentityTForm( LHS, RHS, nTr / 2 );

// Solve

	//SolveWithSquareness( X, LHS, RHS, nTr );
	//SolveWithUnitMag( X, LHS, RHS, nTr );

	printf( "Solve [affines from translations].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X );

	RescaleAll( X, sc );
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
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SolveSystem( vector<double> &X, int nTr )
{
	//SolveSystemStandard( X, nTr );

	//AffineFromTrans( X, nTr );

	AffineFromTransWt( X, nTr );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTransforms(
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
		X[j  ], X[j+1], X[j+2],
		X[j+3], X[j+4], X[j+5] );

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

		double	mag = Magnitude( X, I.itr );

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
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTrakEM(
	double					xmax,
	double					ymax,
	const vector<zsort>		&zs,
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
		double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
		double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

		fprintf( f,
		"\t\t\t<t2_patch\n"
		"\t\t\t\toid=\"%d\"\n"
		"\t\t\t\twidth=\"%d\"\n"
		"\t\t\t\theight=\"%d\"\n"
		"\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
		"\t\t\t\ttitle=\"%.*s\"\n"
		"\t\t\t\ttype=\"%d\"\n"
		"\t\t\t\tfile_path=\"%s\"\n"
		"\t\t\t\to_width=\"%d\"\n"
		"\t\t\t\to_height=\"%d\"\n",
		oid++, gW - offset, gH - offset,
		X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
		s2 - s1, s1, xml_type, p, gW - offset, gH - offset );

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
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteJython(
	const vector<zsort>		&zs,
	const vector<double>	&X,
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
		double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
		double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

		fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f]%s\n",
			p, X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
			(itrf == Ntr ? "" : ",") );
	}

	fprintf( f, "}\n" );
	fclose( f );
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


