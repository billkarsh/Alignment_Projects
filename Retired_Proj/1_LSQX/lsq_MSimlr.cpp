

#include	"lsq_MSimlr.h"

#include	"TrakEM2_UTL.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MSimlr::SetPointPairs(
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

		// R1(p1) - R2(p2) = 0
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
void MSimlr::SetIdentityTForm(
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
// solution output file gArgs.unt_file.
//
void MSimlr::SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
{
/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

	map<MZIDR,TAffine>	M;

	LoadTAffineTbl_RngZ( M, unite_layer, unite_layer, unt_file );

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
/* SolveFromPoints ---------------------------------------------- */
/* --------------------------------------------------------------- */

void MSimlr::SolveFromPoints(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTr )
{
/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

	WriteSolveRead( X, LHS, RHS, "S-FrmPts", nproc, false );
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
void MSimlr::SolveWithUnitMag(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nTR )
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

	WriteSolveRead( X, LHS, RHS, "S-Unimag", nproc, false );
	printf( "\t\t\t\t" );
	PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* RescaleAll ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void MSimlr::RescaleAll(
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

void MSimlr::RotateAll(
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

void MSimlr::NewOriginAll(
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
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Build and solve system of linear equations.
//
// Note:
// All translational variables are scaled down by 'scale' so they
// are sized similarly to the sine/cosine variables. This is only
// to stabilize solver algorithm. We undo the scaling on exit.
//
void MSimlr::SolveSystem( vector<double> &X, int nTr )
{
	double	scale	= 2 * max( gW, gH );
	int		nvars	= nTr * NX;

	printf( "Rgd: %d unknowns; %d constraints.\n",
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

	SolveFromPoints( X, LHS, RHS, nTr );

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
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MSimlr::WriteTransforms(
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
		X[j  ], -X[j+1], X[j+2],
		X[j+1],  X[j  ], X[j+3] );

		if( !bstrings ) {

			fprintf( FOUT, "TAFFINE %d.%d-%d %f %f %f %f %f %f\n",
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

		double	mag = Magnitude( X, I.itr );

		smag += mag;
		smin  = fmin( smin, mag );
		smax  = fmax( smax, mag );

		if( mag < 0.9 ) {
			fprintf( fx, "Low mag %f @ %d.%d-%d\n",
				mag, I.z, I.id, I.rgn );
		}
		else if( mag > 1.1 ) {
			fprintf( fx, "Hi  mag %f @ %d.%d-%d\n",
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

void MSimlr::WriteTrakEM(
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
		double	x_orig = X[j  ]*trim - X[j+1]*trim + X[j+2];
		double	y_orig = X[j+1]*trim + X[j  ]*trim + X[j+3];

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
		X[j], -X[j+1], X[j+1], X[j], x_orig, y_orig,
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

void MSimlr::WriteJython(
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
		double	x_orig = X[j  ]*trim - X[j+1]*trim + X[j+2];
		double	y_orig = X[j+1]*trim + X[j  ]*trim + X[j+3];

		fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f]%s\n",
			path, X[j], -X[j+1], X[j+1], X[j], x_orig, y_orig,
			(itrf == Ntr ? "" : ",") );
	}

	fprintf( f, "}\n" );
	fclose( f );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MSimlr::G2LPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	int		j = itr * NX;
	TAffine	I, T( X[j], -X[j+1], X[j+2], X[j+1], X[j], X[j+3] );
	I.InverseOf( T );
	I.Transform( p );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MSimlr::L2GPoint(
	Point					&p,
	const vector<double>	&X,
	int						itr )
{
	int		j = itr * NX;
	TAffine	T( X[j], -X[j+1], X[j+2], X[j+1], X[j], X[j+3] );
	T.Transform( p );
}


void MSimlr::L2GPoint(
	vector<Point>			&p,
	const vector<double>	&X,
	int						itr )
{
	int		j = itr * NX;
	TAffine	T( X[j], -X[j+1], X[j+2], X[j+1], X[j], X[j+3] );
	T.Transform( p );
}


