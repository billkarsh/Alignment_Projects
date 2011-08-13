

#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"ImageIO.h"	// @@@ temp testing

#include	<fftw3.h>

#include	<math.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	_debugcorr	0






/* --------------------------------------------------------------- */
/* IntegrateImage ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void IntegrateImage(
	vector<double>			&S,
	vector<double>			&S2,
	vector<int>				&nz,
	int						wS,
	int						hS,
	const vector<double>	&I,
	int						wI )
{
	double	t;
	int		x, y, i, j, k, m, N = wS * hS;

	if( S.size() ) {

		S.resize(  N );
		S2.resize( N );
		nz.resize( N );

		memset( &S[0],  0, N * sizeof(double) );
		memset( &S2[0], 0, N * sizeof(double) );
		memset( &nz[0], 0, N * sizeof(int) );
	}
	else {
		S.resize(  N, 0.0 );
		S2.resize( N, 0.0 );
		nz.resize( N, 0 );
	}

// first element
	S[0]	= (t = I[0]);
	S2[0]	= t * t;
	nz[0]	= (t != 0.0);

// first row
	for( x = 1; x < wS; ++x ) {

		S[x]	= S[x-1]  + (t = I[x]);
		S2[x]	= S2[x-1] + t * t;
		nz[x]	= nz[x-1] + (t != 0.0);
	}

// first column
	for( y = 1; y < hS; ++y ) {

		j = wS * y;
		k = wI * y;

		S[j]	= S[j-wS]  + (t = I[k]);
		S2[j]	= S2[j-wS] + t * t;
		nz[j]	= nz[j-wS] + (t != 0.0);
	}

// now all rows and columns
	for( y = 1; y < hS; ++y ) {

		i = wS * y;
		k = wI * y;

		for( x = 1; x < wS; ++x ) {

			j = i + x;
			m = k + x;

			S[j]	= S[j-1]  + S[j-wS]  - S[j-wS-1]  + (t = I[m]);
			S2[j]	= S2[j-1] + S2[j-wS] - S2[j-wS-1] + t * t;
			nz[j]	= nz[j-1] + nz[j-wS] - nz[j-wS-1] + (t != 0.0);
		}
	}
}


static void IntegrateImage(
	vector<int>				&nz,
	int						wS,
	int						hS,
	const vector<double>	&I,
	int						wI )
{
	double	t;
	int		x, y, i, j, k, m, N = wS * hS;

	if( nz.size() ) {

		nz.resize( N );
		memset( &nz[0], 0, N * sizeof(int) );
	}
	else
		nz.resize( N, 0 );

// first element
	nz[0] = (I[0] != 0.0);

// first row
	for( x = 1; x < wS; ++x )
		nz[x] = nz[x-1] + (I[x] != 0.0);

// first column
	for( y = 1; y < hS; ++y ) {

		j = wS * y;
		k = wI * y;

		nz[j] = nz[j-wS] + (I[k] != 0.0);
	}

// now all rows and columns
	for( y = 1; y < hS; ++y ) {

		i = wS * y;
		k = wI * y;

		for( x = 1; x < wS; ++x ) {

			j = i + x;
			m = k + x;

			nz[j] = nz[j-1] + nz[j-wS] - nz[j-wS-1] + (I[m] != 0.0);
		}
	}
}

/* --------------------------------------------------------------- */
/* IntegralTable ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Look up in a cumulative table of doubles.
//
// Inverts IntegrateImage bookkeeping.
//
static double IntegralTable(
	vector<double>	&t,
	int				wS,
	const IBox		&B )
{
	double	rslt = t[B.R + wS*B.T];

	if( B.L > 0 )
		rslt -= t[B.L-1 + wS*B.T];

	if( B.B > 0 )
		rslt -= t[B.R + wS*(B.B-1)];

	if( B.L > 0 && B.B > 0 )
		rslt += t[(B.L-1) + wS*(B.B-1)];

	return rslt;
}


// Look up in a cumulative table of ints
//
// Inverts IntegrateImage bookkeeping.
//
static int IntegralTable(
	vector<int>		&t,
	int				wS,
	const IBox		&B )
{
	int		rslt = t[B.R + wS*B.T];

	if( B.L > 0 )
		rslt -= t[B.L-1 + wS*B.T];

	if( B.B > 0 )
		rslt -= t[B.R + wS*(B.B-1)];

	if( B.L > 0 && B.B > 0 )
		rslt += t[(B.L-1) + wS*(B.B-1)];

	return rslt;
}

/* --------------------------------------------------------------- */
/* DebugLinearCorr ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Compute normalized cross correlation straight from
// the definition. Not efficient - only for debugging.
//
static double DebugLinearCorr(
	FILE			*flog,
	vector<double>	&I1,
	vector<double>	&I2,
	int				wI,
	const IBox		&B1,
	const IBox		&B2 )
{
	double	suma	= 0.0,
			sumb	= 0.0;
	int		nx		= B1.R-B1.L+1;
	int		ny		= B1.T-B1.B+1;
	int		x, y;

	for( x = 0; x < nx; ++x ) {

		for( y = 0; y < ny; ++y ) {

			suma += I1[B1.L+x + wI*(B1.B+y)];
			sumb += I2[B2.L+x + wI*(B2.B+y)];
		}
	}

	double	avga	= suma/nx/ny;
	double	avgb	= sumb/nx/ny;
	double	sumn	= 0.0,
			sumd1	= 0.0,
			sumd2	= 0.0;

	for( x = 0; x < nx; ++x ) {

		for( y = 0; y < ny; ++y ) {

			double	va = I1[B1.L+x + wI*(B1.B+y)] - avga;
			double	vb = I2[B2.L+x + wI*(B2.B+y)] - avgb;

			sumn	+= va * vb;
			sumd1	+= va * va;
			sumd2	+= vb * vb;
		}
	}

	double	prod	= sumd1 * sumd2;
	double	rslt	= (prod < 1.0e-9 ? 0.0 : sumn/sqrt(prod));

	fprintf( flog,
	"\nDebugLinearCorr: s %f %f, a %f %f, ss %f %f %f, r %f.\n",
	suma, sumb, avga, avgb, sumn, sumd1, sumd2, rslt );

	return rslt;
}

/* --------------------------------------------------------------- */
/* LookFFT ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Get FFT value at signed (x,y) coordinate.
//
static double LookFFT(
	const double	*fr,
	int				nx,
	int				ny,
	int				x,
	int				y )
{
	if( x < 0 ) x += nx;
	if( y < 0 ) y += ny;

	return fr[x + nx*y] / (nx * ny);
}

/* --------------------------------------------------------------- */
/* LkFFT --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Get FFT value at signed (x,y) coordinate (no normalization).
//
static double LkFFT(
	const double	*fr,
	int				nx,
	int				ny,
	int				x,
	int				y )
{
	if( x < 0 ) x += nx;
	if( y < 0 ) y += ny;

	return fr[x + nx*y];
}

/* --------------------------------------------------------------- */
/* FFTSize ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return correct storage size for 1D FFT-based correlation.
//
// For 2D and higher cases, the dimensions are independent,
// so call this function separately for each.
//
// Inputs:
//	- nsrcdat:	size of source data
//	- ntrgdat:	size of target data
//	- nneglag:	number of strictly negative lags sought
//	- nposlag:	number of positive lags sought (include zero)
//
// Meaningful values for nneglag are: [1,...,nsrcdat-1].
// Meaningful values for nposlag are: [1,...,ntrgdat].
//
// Returns:
//	Storage allocation size that should be used for each of
//	{source data, target data, correlation result}. The value
//	accounts for needed zero padding and is a power of two.
//
// Choosing storage size
// ---------------------
// The first important constraint is that the two input data
// arrays and the single result array have the same length, so
// we seek a single value (N) that satisfies all requirements.
//
// Obviously N must be large enough to hold either input data:
//
//	(1) N >= max( nsrcdat, ntrgdat ).
//
// Next consider how results (R) are packed into N elements:
//	R[0]		= lag 0
//	R[1]		= lag +1
//	...
//	R[N/2-1]	= lag +(N/2-1)	=> N/2 lags are >= zero.
//
//
//	R[N-1]		= lag -1
//	R[N-2]		= lag -2
//	...
//	R[N/2]		= lag -(N/2)	=> N/2 lags are < zero.
//
// We have no choice about the fact that one half is positive and
// the other negative, so if we want a given number of positive or
// negative values, we have to increase N so they fit in their
// designated half. Which means:
//
//	(2) N/2 >= max( nneglag, nposlag ).
//
// Therefore, we have:
//
//	(3) N = max( nsrcdat, ntrgdat, nneglag*2, nposlag*2 ).
//
// Finally, we increase the size to the nearest power of two
// for best performance (more below).
//
// Addressing result array
// -----------------------
// The valid ranges of the result array---
//
//	Negative lags:	R[N-nneglag] ... R[N-1], inclusive
//	Positive lags:	R[0] ... R[nposlag-1], inclusive.
//
// Why a power of two?
// -------------------
// Although the FFTW library works correctly for storage sizes
// that are just large enough to accommodate wrap-around, the
// performance is much better for powers of two. Here are some
// representative time measurements (arbs)---
//
//	power of 2 dims:	1
//	even dims:			1.75
//	odd dims:			2.20+
//
int FFTSize( int nsrcdat, int ntrgdat, int nneglag, int nposlag )
{
// minimum padded size

	int	n = 2 * max( nneglag, nposlag );

	n = max( n, max( nsrcdat, ntrgdat ) );

// nearest power of two

	return CeilPow2( n );
}

/* --------------------------------------------------------------- */
/* FFTSizeSQR ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Older version of FFTSize that returns value to be used
// for both axes, that is, for square FFT workspaces.
//
int FFTSizeSQR( int w1, int h1, int w2, int h2 )
{
	return CeilPow2( max( w1 + w2, h1 + h2 ) - 1 );
}

/* --------------------------------------------------------------- */
/* FFT_2D -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Forward FFT of 2D data (real to complex).
//
// Assumes input data in row-major order. That is,
// ordered like a C-array: in[Nslow][Nfast].
//
int FFT_2D(
	vector<CD>				&out,
	const vector<double>	&in,
	int						Nfast,
	int						Nslow )
{
	fftw_plan	p;
	int			M = Nslow * (Nfast/2 + 1);

	out.resize( M );

	p = fftw_plan_dft_r2c_2d( Nslow, Nfast, (double*)&in[0],
			(double (*)[2])&out[0], FFTW_ESTIMATE );

	fftw_execute( p );
	fftw_destroy_plan( p );

	return M;
}

/* --------------------------------------------------------------- */
/* IFT_2D -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Inverse FFT of 2D data (complex to real).
//
// Creates output data in row-major order. That is,
// ordered like a C-array: out[Nslow][Nfast].
//
void IFT_2D(
	vector<double>			&out,
	const vector<CD>		&in,
	int						Nfast,
	int						Nslow )
{
	fftw_plan	p;
	int			N = Nslow * Nfast;

	out.resize( N );

	p = fftw_plan_dft_c2r_2d( Nslow, Nfast, (double (*)[2])&in[0],
			&out[0], FFTW_ESTIMATE );

	fftw_execute( p );
	fftw_destroy_plan( p );
}

/* --------------------------------------------------------------- */
/* ParabPeakFFT -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Improve estimate of FFT peak position (xpk, ypk) using
// fit of parabola to three points through peak. (d) is how
// many pixels away from peak to get samples.
//
// If, for XZ-plane, Z(x) = a(x - b)^2 + c, and,
//
// z0 = Z(xpk - d),
// z1 = Z(xpk),
// z2 = Z(xpk + d), then,
//
//        d      z0 - z2
// xpk' = - * -------------- + xpk
//        2    z0 + z2 - 2z1
//
// Note: Although solution is analytic, the quality of the result
// depends upon the quality of the data and wild displacements do
// occur. Therefore we limit displacements to some small multiple
// of the sampling distance.
//
void ParabPeakFFT(
	double			&xpk,
	double			&ypk,
	int				d,
	const double	*I,
	int				nx,
	int				ny )
{
	int		ix = (int)xpk,
			iy = (int)ypk;
	double	z1 = LkFFT( I, nx, ny, ix, iy ),
			z0,
			z2;

	if( (z0 = LkFFT( I, nx, ny, ix - d, iy )) &&
		(z2 = LkFFT( I, nx, ny, ix + d, iy )) ) {

		z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

		if( fabs( z0 ) <= 8 * d )
			xpk += z0;
	}

	if( (z0 = LkFFT( I, nx, ny, ix, iy - d )) &&
		(z2 = LkFFT( I, nx, ny, ix, iy + d )) ) {

		z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

		if( fabs( z0 ) <= 8 * d )
			ypk += z0;
	}
}

/* --------------------------------------------------------------- */
/* PrintCorLandscape --------------------------------------------- */
/* --------------------------------------------------------------- */

void PrintCorLandscape(
	double			biggest,
	int				bigx,
	int				bigy,
	int				Ox,
	int				Oy,
	int				radius,
	int				lim,
	int				step,
	const double	*I,
	int				nx,
	int				ny,
	double			norm,
	FILE			*flog )
{
	fprintf( flog,
		"Landscape: Max %f at (%d, %d); Ox Oy radius = %d %d %d.\n",
		biggest, bigx, bigy, Ox, Oy, radius );

	for( int iy = -lim; iy <= lim; iy += step ) {

		int	ay = bigy + iy;

		if( ay < 0 )
			ay += ny;

		for( int ix = -lim; ix <= lim; ix += step ) {

			int	ax = bigx + ix;

			if( ax < 0 )
				ax += nx;

			fprintf( flog, "%12.6f ", I[ax + nx*ay]/norm );
		}

		fprintf( flog, "\n" );
	}
}

/* --------------------------------------------------------------- */
/* RVectors ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return linear correlation coefficient.
//
double RVectors(
	const vector<double>	&a,
	const vector<double>	&b )
{
	double	A = 0.0, B = 0.0, AA = 0.0, BB = 0.0, AB = 0.0;
	int		n = a.size();

	for( int i = 0; i < n; ++i ) {

		A	+= a[i];
		B	+= b[i];
		AA	+= a[i] * a[i];
		BB	+= b[i] * b[i];
		AB	+= a[i] * b[i];
	}

	return (n*AB - A*B) / sqrt( (n*AA - A*A) * (n*BB - B*B) );
}

/* --------------------------------------------------------------- */
/* CLinCorr ------------------------------------------------------ */
/* --------------------------------------------------------------- */

class CLinCorr {

private:
	vector<double>	I1,    I2,
					i1sum, i1sum2,
					i2sum, i2sum2;
	vector<int>		i1nz,  i2nz;
	FILE			*flog;
	int				w1,  h1,
					w2,  h2,
					Nx,  Nxy;
	IBox			OL1, OL2;
	double*			rslt;
	int				ir,
					dx,  dy,
					olw, olh,
					i1c, i2c;

public:
	void Initialize(
		FILE					*flog,
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny );

	int CheckSize(
		vector<double>	&rslt,
		int				irslt,
		EvalType		LegalRgn,
		void*			arglr,
		int				dx,
		int				dy );

	int CheckDensity(
		EvalType		LegalCnt,
		void*			arglc );

	double GetCorr();
	int    SizeIndex();
};


void CLinCorr::Initialize(
		FILE					*flog,
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny )
{
	this->I1	= I1;
	this->I2	= I2;
	this->flog	= flog;
	this->w1	= w1;
	this->h1	= h1;
	this->w2	= w2;
	this->h2	= h2;
	this->Nx	= Nx;
	Nxy			= Nx * Ny;

	IntegrateImage( i1sum, i1sum2, i1nz, w1, h1, I1, Nx );
	IntegrateImage( i2sum, i2sum2, i2nz, w2, h2, I2, Nx );
}


int CLinCorr::CheckSize(
	vector<double>	&rslt,
	int				irslt,
	EvalType		LegalRgn,
	void*			arglr,
	int				dx,
	int				dy )
{
	int		ok = true;

	this->rslt	= &rslt[0];
	ir			= irslt;
	this->dx	= dx;
	this->dy	= dy;

	BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

	olw = OL1.R - OL1.L + 1;
	olh = OL1.T - OL1.B + 1;

	if( LegalRgn && !LegalRgn( olw, olh, arglr ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


int CLinCorr::CheckDensity(
	EvalType		LegalCnt,
	void*			arglc )
{
	int		ok = true;

	i1c = IntegralTable( i1nz, w1, OL1 );
	i2c = IntegralTable( i2nz, w2, OL2 );

	if( LegalCnt && !LegalCnt( i1c, i2c, arglc ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


double CLinCorr::GetCorr()
{
	double	n = olw * olh;

	double	im1sum	= IntegralTable( i1sum,  w1, OL1 );
	double	im1sum2	= IntegralTable( i1sum2, w1, OL1 );
	double	im2sum	= IntegralTable( i2sum,  w2, OL2 );
	double	im2sum2	= IntegralTable( i2sum2, w2, OL2 );

	double	num	= n * rslt[ir] / Nxy - im1sum * im2sum;
	double	d1	= n * im1sum2 - im1sum * im1sum;
	double	d2	= n * im2sum2 - im2sum * im2sum;
	double	d	= d1 * d2;
	double	r	= (d < n * n * 1.0E-9 ? 0.0 : num / sqrt( d ));

	if( abs( r ) > 1.0001 ) {

		fprintf( flog,
		"NormCorr: Very odd - ir=%d rslt[i]=%f olap_area=%d.\n",
		ir, r, n );

		fprintf( flog,
		"NormCorr: shift %d %d, i1 (%d %d) to (%d %d),"
		" i2 (%d %d) to (%d %d).\n",
		dx, dy,
		OL1.L, OL1.B, OL1.R, OL1.T,
		OL2.L, OL2.B, OL2.R, OL2.T );

		fprintf( flog,
		"NormCorr: sums:      %f %f %f %f.\n"
		"NormCorr: num d1 d2: %f %f %f.\n",
		im1sum, im1sum2, im2sum, im2sum2, num, d1, d2 );

		DebugLinearCorr( flog, I1, I2, Nx, OL1, OL2 );
		exit(44);
	}

	return rslt[ir] = r;
}


int CLinCorr::SizeIndex()
{
	return (int)(10.0 * log( max( 1, min(i1c,i2c) ) ));
}

/* --------------------------------------------------------------- */
/* CCrossCorr ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CCrossCorr {

private:
	vector<int>		i1nz,  i2nz;
	int				w1,  h1,
					w2,  h2,
					Nxy;
	IBox			OL1, OL2;
	double*			rslt;
	int				ir,
					i1c, i2c;

public:
	void Initialize(
		FILE					*flog,
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny );

	int CheckSize(
		vector<double>	&rslt,
		int				irslt,
		EvalType		LegalRgn,
		void*			arglr,
		int				dx,
		int				dy );

	int CheckDensity(
		EvalType		LegalCnt,
		void*			arglc );

	double GetCorr();
	int    SizeIndex();
};


void CCrossCorr::Initialize(
		FILE					*flog,
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny )
{
	this->w1	= w1;
	this->h1	= h1;
	this->w2	= w2;
	this->h2	= h2;
	Nxy			= Nx * Ny;

	IntegrateImage( i1nz, w1, h1, I1, Nx );
	IntegrateImage( i2nz, w2, h2, I2, Nx );
}


int CCrossCorr::CheckSize(
	vector<double>	&rslt,
	int				irslt,
	EvalType		LegalRgn,
	void*			arglr,
	int				dx,
	int				dy )
{
	int		olw, olh, ok = true;

	this->rslt	= &rslt[0];
	ir			= irslt;

	BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

	olw = OL1.R - OL1.L + 1;
	olh = OL1.T - OL1.B + 1;

	if( LegalRgn && !LegalRgn( olw, olh, arglr ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


int CCrossCorr::CheckDensity(
	EvalType		LegalCnt,
	void*			arglc )
{
	int		ok = true;

	i1c = IntegralTable( i1nz, w1, OL1 );
	i2c = IntegralTable( i2nz, w2, OL2 );

	if( LegalCnt && !LegalCnt( i1c, i2c, arglc ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


double CCrossCorr::GetCorr()
{
	return rslt[ir] /= (double)Nxy * i1c;
}


int CCrossCorr::SizeIndex()
{
	return (int)(10.0 * log( max( 1, min(i1c,i2c) ) ));
}

/* --------------------------------------------------------------- */
/* CorrPatches --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return normalized cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// Search confined to disc: (origin, radius) = {Ox, Oy, radius}.
//
// 'fft2' is a cache of the patch2 FFT. On entry, if fft2 has
// the correct size it is used. Otherwise recomputed here.
//
double CorrPatches(
	FILE					*flog,
	int						verbose,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						Ox,
	int						Oy,
	int						radius,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 )
{
// Bounding boxes of point lists

	IBox	B1, B2;
	int		w1, h1, w2, h2;

	BBoxFromPoints( B1, ip1 );
	BBoxFromPoints( B2, ip2 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

	w2 = B2.R - B2.L + 1;
	h2 = B2.T - B2.B + 1;

	if( verbose ) {

		fprintf( flog,
		"NormCorr: region size is [%d %d] in x, [%d %d] in y.\n",
		B1.L, B1.R, B1.B, B1.T );

		fprintf( flog,
		"NormCorr: target size is [%d %d] in x, [%d %d] in y.\n",
		B2.L, B2.R, B2.B, B2.T );
	}

// Get array sizes (Nx,Ny) and FFT size M

	int	nnegx	= w1 - 1,
		nnegy	= h1 - 1,
		nposx	= w2,
		nposy	= h2,
		Nx		= FFTSize( w1, w2, nnegx, nposx ),
		Ny		= FFTSize( h1, h2, nnegy, nposy ),
		M		= Ny*(Nx/2+1);

	if( verbose )
		fprintf( flog, "NormCorr: Nx = %d, Ny = %d.\n", Nx, Ny );

// Create images from point lists.

	vector<double>	i1, i2;

	ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
	ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1;

	if( fft2.size() != M )
		FFT_2D( fft2, i2, Nx, Ny );

	FFT_2D( fft1, i1, Nx, Ny );

	for( int i = 0; i < M; ++i )
		fft1[i] = fft2[i] * conj( fft1[i] );

	IFT_2D( rslt, fft1, Nx, Ny );

// Create array indexed by 'size': int( 10 * log( overlap_size ) ).
// Each element contains the best rslt[i] at that size index.
// The idea is to look for correlation peaks that are nearly
// as good as the best, but that derive from larger overlaps.

	int lmax = (int)(10.0*log(double(Nx*Ny))) + 1;
	vector<double>	max_by_size( lmax, 0.0 );

// Prepare correlation calculator

	CLinCorr	ccalc;	// uses lin corr coeff
//	CCrossCorr	ccalc;	// uses 1/n * SUM(a*b)

	ccalc.Initialize( flog, i1, w1, h1, i2, w2, h2, Nx, Ny );

// Now for the conventional biggest

	double	biggest	= -1.0E30;
	int		bigx	= -1;
	int		bigy	= -1;
	int		rsqr	= radius * radius;

	for( int iy = 0; iy < Ny; ++iy ) {

		int	dy, y = iy;

		if( y >= Ny - nnegy )	// negative zone
			y -= Ny;
		else if( y >= nposy ) {	// dead zone
			memset( &rslt[Nx*iy], 0, Nx*sizeof(double) );
			continue;
		}

		dy = y - Oy;

		for( int ix = 0; ix < Nx; ++ix ) {

			int	dx, x = ix, i = ix+Nx*iy;

			if( x >= Nx - nnegx )	// negative zone
				x -= Nx;
			else if( x >= nposx )	// dead zone
				goto skip;

			dx = x - Ox;

			if( dx*dx + dy*dy > rsqr ) {
skip:
				rslt[i] = 0.0;
				continue;
			}

			// Linear correlation coeff for pixel {x, y}

			if( !ccalc.CheckSize( rslt, i, LegalRgn, arglr, x, y ) )
				continue;

			if( !ccalc.CheckDensity( LegalCnt, arglc ) )
				continue;

			double	r = ccalc.GetCorr();

			if( r > biggest ) {
				biggest	= r;
				bigx	= x;
				bigy	= y;
			}

			// Update max_by_size

			int im = ccalc.SizeIndex();

			if( r > max_by_size[im] )
				max_by_size[im] = r;
		}
	}

// Reports

	if( biggest < -2.0 ) {

		fprintf( flog, "NormCorr: No legal subregions at all...\n");

		// for sake of grep searches

		fprintf( flog,
		"NormCorr: Maximum correlation of %f at [%d,%d].\n",
		0.0, 0, 0 );

		return 0.0;
	}

	if( verbose ) {

		PrintCorLandscape( biggest, bigx, bigy, Ox, Oy, radius,
		3, 1, &rslt[0], Nx, Ny, 1.0, flog );

		fprintf( flog, "----v-center-v---\n" );
		PrintCorLandscape( biggest, 32, 32, Ox, Oy, radius,
		3, 1, &rslt[0], Nx, Ny, 1.0, flog );
	}

// For debugging, print the results of correlation by region size.
// If we find a larger region with a correlation almost as good,
// that's potentially a better match.

	double	bc		= -10.0;		// biggest correlation
	int		we		= 0;			// which entry had it?
	int		nmax	= max_by_size.size();

	for( int i = 0; i < nmax; ++i ) {

		if( max_by_size[i] > bc ) {

			bc = max_by_size[i];
			we = i;
		}
	}

// Now print the entries with bigger size and comparable correlation
// Mpy by 64 to get into pixel counts in the 2K working image.

	double PT = 0.8 * bc;	// Print threshold

	for( int i = we + 1; i < nmax; ++i ) {

		if( max_by_size[i] >= PT ) {

			int i1 = (int)ceil(  64.0 * exp( i/10.0 ) );
			int i2 = (int)floor( 64.0 * exp( (i+1)/10.0 ) );

			fprintf( flog,
			"NormCorr: Possible bigger area match:"
			" %8d - %8d : %8.2f.\n", i1, i2, max_by_size[i] );
		}
	}

	fprintf( flog,
	"NormCorr: Maximum correlation of %f at [%d,%d].\n",
	biggest, bigx, bigy );

// Local peak test

#if 0
	double limit = 0.9 * biggest;

	for( int x = -w1 + 1; x < w2 - 1; ++x ) {
		for( int y = -h1 + 1; y < h2 - 1; ++y ) {
			double a = LkFFT( rslt, Nx, Ny, x, y );
			if( a > limit &&
				a > LkFFT( rslt, Nx, Ny, x-1, y ) &&
				a > LkFFT( rslt, Nx, Ny, x+1, y ) &&
				a > LkFFT( rslt, Nx, Ny, x, y-1 ) &&
				a > LkFFT( rslt, Nx, Ny, x, y+1 ) ) {

				fprintf( flog,
				"NormCorr: Local max at %d %d, value %f.\n",
				x, y, a );
			}
		}
	}
#endif

// Interpolate peak

	dx = bigx;
	dy = bigy;
	ParabPeakFFT( dx, dy, 1, &rslt[0], Nx, Ny );

	dx += B2.L - B1.L;
	dy += B2.B - B1.B;

	return biggest;
}

/* --------------------------------------------------------------- */
/* CorrPatchToImage ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Return the normalized cross-correlation and (dx, dy), which
// must be added to points in patch1 to match image2.
//
// Search confined to disc: (origin, radius) = {Ox, Oy, radius}.
//
// IMPORTANT:
// Although a normalized final result is returned, the intermediate
// values are not normalized, so comparison over widely differing
// insection areas is unreliable. CorrPatches is preferred.
//
double CorrPatchToImage(
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<double>	&i2,
	int						Ox,
	int						Oy,
	int						radius,
	bool					bFilter )
{
	int		Nsqr	= i2.size(),
			np1		= ip1.size(),
			N		= (int)sqrt( Nsqr ),
			i, M;
	double	norm	= (double)np1 * Nsqr;

// Bounding box of ip1

	IBox	B1;
	int		w1, h1;

	BBoxFromPoints( B1, ip1 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

// Effective {w2,h2} of caller-embedded i2 data

	int		w2 = 0, h2 = 0;

	for( i = 0; i < Nsqr; ++i ) {

		if( i2[i] != 0.0 ) {

			int	y = i / N;
			int	x = i - N * y;

			if( y > h2 )
				h2 = y;

			if( x > w2 )
				w2 = x;
		}
	}

	++w2;
	++h2;

// Image1 from point list

	vector<double>	i1( Nsqr, 0.0 );

	for( i = 0; i < np1; ++i ) {

		int	x = (int)ip1[i].x;
		int	y = (int)ip1[i].y;

		i1[x + N*y] = iv1[i];
	}

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1, fft2;

	M =	FFT_2D( fft2, i2, N, N );
		FFT_2D( fft1, i1, N, N );

	if( bFilter ) {

		// an experiment to filter out high freq noise

		int	yrow = N/2 + 1;

		for( i = 0; i < M; ++i ) {

			double	rad2, filt;
			int		x = i / yrow;
			int		y = i - yrow*x;

			// note that y never exceeds 2048

			if( x >= N/2 )
				x -= N;

			// keep first 30 harmonics or so...

			rad2 = ((double)x*x + (double)y*y) / 10000.0;

			if( rad2 > 10.0 )
				filt = 0.0;
			else
				filt = exp( -rad2 );

			fft1[i] = fft2[i] * conj( fft1[i] ) * filt;
		}
	}
	else {

		for( i = 0; i < M; ++i )
			fft1[i] = fft2[i] * conj( fft1[i] );
	}

	IFT_2D( rslt, fft1, N, N );

// Max lag

	double	biggest	= -1.0E30;
	int		bigx	= -1;
	int		bigy	= -1;
	int		rsqr	= radius * radius;

	for( i = 0; i < Nsqr; ++i ) {

		int	y = i / N;
		int	x = i - N * y;
		int	dx, dy;

		if( y > N - h1 )	// negative zone
			y -= N;
		else if( y >= h2 )	// dead zone
			continue;

		if( x > N - w1 )	// negative zone
			x -= N;
		else if( x >= w2 )	// dead zone
			continue;

		dx = x - Ox;
		dy = y - Oy;

		if( dx*dx + dy*dy > rsqr )
			continue;

		if( rslt[i] > biggest ) {

			biggest	= rslt[i];
			bigx	= x;
			bigy	= y;
		}
	}

	biggest /= norm;

	if( biggest > 1.0001 ) {

		printf( "FindCor: Very odd - norm-rslt x y %f %d %d.\n",
			biggest, bigx, bigy );
		exit(44);
	}

// Reports

	PrintCorLandscape( biggest, bigx, bigy, Ox, Oy, radius,
		3, 1, &rslt[0], N, N, norm, stdout );

	printf( "FindCor: Maximum correlation %f at (%d, %d).\n",
		biggest, bigx, bigy );

// Interpolate peak

	dx = bigx;
	dy = bigy;
	ParabPeakFFT( dx, dy, 1, &rslt[0], N, N );

	printf( "FindCor: Interpolated max at (%.3f, %.3f).\n", dx, dy );

	return biggest;
}

/* --------------------------------------------------------------- */
/* GradDescStep -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Let region A be represented in barycentric coords using
// {ac = control_points, am = multipliers, av = values}.
//
// Let region B be represented as a raster {bimg, w, h}.
//
// The objective is to optimize the alignment of A to B by tweaking
// the control points to MAXIMIZE cross correlation measure R.
//
// A gradient descent scheme is implemented by three functions:
//
// (1) ImproveControlPts() is the driver that adjusts step size and
// monitors convergence.
//
// (2) Cross corr. is evaluated using R = CorrVectors( av, bnew ).
//
// (3) The updated set of B-values (bnew), and a gradient vector
// dR/dc are both generated by GradDescStep().
//
// Notes
// -----
// The caller must first transform the {ac} from A- to B-coords
// before calling driver ImproveControlPts().
//
// As noted in the comments for ImproveControlPts(), the {am}
// for the ith point must include multipliers for each control
// point {ac}, where all except three are expected to be zero.
//
// At first glance one might guess that moving the control points
// ac produces an updated set of av. However, one should think of
// the av and the am as invariants in the optimization. Rather, as
// the A control points move, the values av move with the deformed
// A-region, but those values now line up with new (interpolated)
// B-points.
//
// R = SUMi(ai * bi). Remembering that it is really the bi that
// depend upon the control points, we get that the dependence of
// R on the x-coordinate of the k-th control point is:
//
// dR/dc(k,x) = SUMi(dR/dbi * dbi/dx * dx/dc(k))
//
// = SUMi(ai * dbi/dx * am(i,k)).
//
void GradDescStep(
	vector<double>					&bnew,
	vector<Point>					&dRdc,
	const vector<Point>				&ac,
	const vector<vector<double> >	&am,
	const vector<double>			&av,
	const vector<double>			&bimg,
	int								w,
	int								h )
{
	int		nm = am.size();
	int		nc = ac.size();

// Initialize

	bnew.resize( nm );
	dRdc.resize( nc );

	Point	zero( 0.0, 0.0 );

	for( int i = 0; i < nc; ++i )
		dRdc[i] = zero;

// Sum over all points in A

	for( int i = 0; i < nm; ++i ) {

		double	x = 0.0, y = 0.0;

		for( int j = 0; j < nc; ++j ) {
			x += am[i][j]*ac[j].x;
			y += am[i][j]*ac[j].y;
		}

		// the 4 points surrounding (x,y)

		if( x <  0.0	||
			x >= w - 1	||
			y <  0.0	||
			y >= h - 1 ) {

			bnew[i] = 0.0;	// anything outside bimg is 0.0
			continue;
		}

		int	xl = (int)x;
		int	xr = xl + 1;
		int	yl = (int)y;
		int	yu = yl + 1;

		// interpolate

		double alpha	= x - xl;
		double beta		= y - yl;
		double ll		= bimg[w*yl + xl];
		double ul		= bimg[w*yu + xl];
		double ur		= bimg[w*yu + xr];
		double lr		= bimg[w*yl + xr];
		double t		= lr - ll;
		double u		= ul - ll;
		double v		= ll - lr - ul + ur;

		bnew[i] = ll + alpha * t + beta * u + alpha*beta * v;

		// dbi/dx

		double dbdx = t +  beta * v;
		double dbdy = u + alpha * v;

		// dR/dc

		dbdx *= av[i];
		dbdy *= av[i];

		for( int j = 0; j < nc; ++j ) {
			dRdc[j].x += am[i][j]*dbdx;
			dRdc[j].y += am[i][j]*dbdy;
		}
	}

	NormalizeNonZeros( bnew );
}

/* --------------------------------------------------------------- */
/* CorrVectors --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return cross correlation between two vectors,
// assuming each has mean 0 and std 1.
//
// Since we are using this for correlation,
// just check the non-zero pixels.
//
// Returns the number of non-zero pixels, if requested.
//
double CorrVectors(
	FILE					*flog,
	const vector<double>	&a,
	const vector<double>	&b,
	int						*nnz )
{
	double	sum2	= 0.0;
	int		nz		= 0;	// number of 0 entries
	int		N		= a.size();

	if( !flog )
		flog = stdout;

	if( N != b.size() ) {
		fprintf( flog,
		"CorrVectors: Sizes differ! %d %d.\n", N, b.size() );
		exit( 42 );
	}

	for( int i = 0; i < N; ++i ) {

		double	prod = a[i]*b[i];

		sum2 += prod;

		if( fabs( prod ) < 1.0E-8 )
			++nz;
	}

//	fprintf( flog, "CorrVectors: %d of %d were small.\n", nz, N );

	N -= nz;

	if( nnz )
		*nnz = N;

	return sum2 / N;
}

/* --------------------------------------------------------------- */
/* ImproveControlPts --------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintControlPoints( FILE *flog, const vector<Point> &c )
{
	int		nc = c.size();

	fprintf( flog, "\nControl points:" );

	for( int i = 0; i < nc; ++i )
		fprintf( flog, "(%.3f %.3f) ", c[i].x, c[i].y );

	fprintf( flog, "\n" );
}


static void PrintPixels(
	FILE	*flog,
	const vector<vector<double> >	&am,
	const vector<double>			&av,
	const vector<double>			&bnew )
{
	int	nm = am.size();

	for( int i = 0; i < nm; ++i ) {

		int	nc = am[i].size();

		fprintf( flog, "---i=%d\n", i );

		for( int j = 0; j < nc; ++j )
			fprintf( flog, "%.4f ", am[i][j] );

		fprintf( flog, "\n av=%f bnew=%f\n", av[i], bnew[i] );
	}
}


// Driver function to improve correlation. Uses gradient descent
// to tweak locations of control points (See GradDescStep).
//
// Return best correlation obtained.
//
// ac			- A-region control points (in B-coord system)
// am			- control point multipliers (see IMPORTANT note)
// av			- A-values; << MUST BE NORMALIZED >>
// bimg			- B-raster mapped to
// w, h			- B-raster dims
// flog			- log file
// describe		- string describing caller context
// iniThresh	- required initial threshold
// finThresh	- if negative, flag to disable deformation...
//				- if positive, information in printed messages
//
// IMPORTANT:
// The usual expectation is that there would be exactly three
// multipliers per point (assuming the triangle is known). But
// in this code we carry as many multipliers as control points
// and set them all zero except the relevant three. This is done
// so that each point is expressed as a function of all control
// points, and we can thereby calculate changes in correlation
// as a function of changes in control points (mesh distortion).
//
double ImproveControlPts(
	vector<Point>					&ac,
	const vector<vector<double> >	&am,
	const vector<double>			&av,
	const vector<double>			&bimg,
	int								w,
	int								h,
	FILE							*flog,
	const char						*describe,
	double							iniThresh,
	double							finThresh )
{
	vector<double>	bnew;
	vector<Point>	dRdc;
	double			corr, corr_last;
	int				nc		= ac.size(),
					inarow	= 0;

// Initial state

	GradDescStep( bnew, dRdc, ac, am, av, bimg, w, h );
	corr_last = corr = CorrVectors( flog, av, bnew );

	fprintf( flog,
	"STAT: ImproveCpt: Initial %s correlation %f (%d pixels).\n",
	describe, corr, av.size() );

// Plausibility check

	if( corr < iniThresh ) {

		fprintf( flog,
		"STAT: ImproveCpt: Correlation %f less than %f at start.\n",
		corr, iniThresh );

		PrintControlPoints( flog, ac );
		//PrintPixels( flog, am, av, bnew );

		return 0.0;
	}

// Skip optimizing if finThresh < 0

	if( finThresh < 0 ) {

		fprintf( flog,
		"STAT: ImproveCpt: Skipping optimizer; final corr %f.\n",
		corr );

		return corr;
	}

// Try to tweak the control points for a good match

	for( double step = 10.0; step > 0.05; ) {

//		PrintControlPoints( flog, ac );

		fprintf( flog, "corr=%f\tstep=%f\n", corr, step );

		// compute gradient length factor S(step size)

		double	S = 0.0;

		for( int i = 0; i < nc; ++i )
			S += dRdc[i].RSqr();

		if( !S ) {
			fprintf( flog, "*** ALL DERIVS ZERO.\n" );
			break;
		}

		S = step / sqrt( S );

		// update control points and correlation

		vector<Point>	ac_new = ac;
		vector<Point>	dRdc_new;
		double			c_new;

		for( int i = 0; i < nc; ++i ) {
			ac_new[i].x += S * dRdc[i].x;
			ac_new[i].y += S * dRdc[i].y;
		}

		// new spot better?

		GradDescStep( bnew, dRdc_new, ac_new, am, av, bimg, w, h );
		c_new = CorrVectors( flog, av, bnew );

		if( c_new > corr ) {

			ac			= ac_new;
			dRdc		= dRdc_new;
			corr_last	= corr;
			corr		= c_new;
		}
		else
			step /= 2.0;

		// converged?

		if( step <= 2.0 && corr - corr_last < 0.0001 ) {

			if( ++inarow >= 3 )
				break;
		}
		else
			inarow = 0;
	}

	fprintf( flog,
	"STAT: ImproveCpt: Final %s correlation %f, (threshold %f).\n",
	describe, corr, finThresh );

	return corr;
}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

/* --------------------------------------------------------------- */
/* CRQ ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CRQ {

private:
	vector<double>	i1sum, i1sum2,
					i2sum, i2sum2;
	vector<int>		i1nz,  i2nz;
	int				w1,  h1,
					w2,  h2,
					Nxy;
	IBox			OL1, OL2;
	int				olw, olh,
					i1c, i2c;

public:
	void Initialize(
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny );

	int CheckSize(
		EvalType		LegalRgn,
		void*			arglr,
		int				dx,
		int				dy );

	int CheckDensity(
		EvalType		LegalCnt,
		void*			arglc );

	double CalcR( double rslt );
};


void CRQ::Initialize(
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny )
{
	this->w1	= w1;
	this->h1	= h1;
	this->w2	= w2;
	this->h2	= h2;
	Nxy			= Nx * Ny;

	IntegrateImage( i1sum, i1sum2, i1nz, w1, h1, I1, Nx );
	IntegrateImage( i2sum, i2sum2, i2nz, w2, h2, I2, Nx );
}


int CRQ::CheckSize(
	EvalType		LegalRgn,
	void*			arglr,
	int				dx,
	int				dy )
{
	int		ok = true;

	BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

	olw = OL1.R - OL1.L + 1;
	olh = OL1.T - OL1.B + 1;

	if( LegalRgn && !LegalRgn( olw, olh, arglr ) )
		ok = false;

	return ok;
}


int CRQ::CheckDensity(
	EvalType		LegalCnt,
	void*			arglc )
{
	int		ok = true;

	i1c = IntegralTable( i1nz, w1, OL1 );
	i2c = IntegralTable( i2nz, w2, OL2 );

	if( LegalCnt && !LegalCnt( i1c, i2c, arglc ) )
		ok = false;

	return ok;
}


double CRQ::CalcR( double rslt )
{
	double	n = olw * olh;

	double	im1sum	= IntegralTable( i1sum,  w1, OL1 );
	double	im1sum2	= IntegralTable( i1sum2, w1, OL1 );
	double	im2sum	= IntegralTable( i2sum,  w2, OL2 );
	double	im2sum2	= IntegralTable( i2sum2, w2, OL2 );

	double	num	= n * rslt / Nxy - im1sum * im2sum;
	double	d1	= n * im1sum2 - im1sum * im1sum;
	double	d2	= n * im2sum2 - im2sum * im2sum;
	double	d	= d1 * d2;
	double	r	= (d < n * n * 1.0E-9 ? 0.0 : num / sqrt( d ));

	return r;
}

/* --------------------------------------------------------------- */
/* CXCQ ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CXCQ {

private:
	vector<int>		i1nz,  i2nz;
	int				w1,  h1,
					w2,  h2,
					Nxy;
	IBox			OL1, OL2;
	int				i1c, i2c;

public:
	void Initialize(
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny );

	int CheckSize(
		EvalType		LegalRgn,
		void*			arglr,
		int				dx,
		int				dy );

	int CheckDensity(
		EvalType		LegalCnt,
		void*			arglc );

	double GetNorm();
};


void CXCQ::Initialize(
		const vector<double>	I1,
		int						w1,
		int						h1,
		const vector<double>	I2,
		int						w2,
		int						h2,
		int						Nx,
		int						Ny )
{
	this->w1	= w1;
	this->h1	= h1;
	this->w2	= w2;
	this->h2	= h2;
	Nxy			= Nx * Ny;

	IntegrateImage( i1nz, w1, h1, I1, Nx );
	IntegrateImage( i2nz, w2, h2, I2, Nx );
}


int CXCQ::CheckSize(
	EvalType		LegalRgn,
	void*			arglr,
	int				dx,
	int				dy )
{
	int		olw, olh, ok = true;

	BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

	olw = OL1.R - OL1.L + 1;
	olh = OL1.T - OL1.B + 1;

	if( LegalRgn && !LegalRgn( olw, olh, arglr ) )
		ok = false;

	return ok;
}


int CXCQ::CheckDensity(
	EvalType		LegalCnt,
	void*			arglc )
{
	int		ok = true;

	i1c = IntegralTable( i1nz, w1, OL1 );
	i2c = IntegralTable( i2nz, w2, OL2 );

	if( LegalCnt && !LegalCnt( i1c, i2c, arglc ) )
		ok = false;

	return ok;
}


double CXCQ::GetNorm()
{
	return (double)Nxy * i1c;
}

/* --------------------------------------------------------------- */
/* OffsetResults ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void OffsetResults( vector<double> &R )
{
	double	ave	= 0.0;
	int		na	= 0,
			nr	= R.size();

	for( int i = 0; i < nr; ++i ) {

		if( R[i] ) {
			ave += R[i];
			++na;
		}
	}

	if( na ) {

		ave /= na;

		for( int i = 0; i < nr; ++i ) {

			if( R[i] )
				R[i] -= ave;
		}
	}
}

/* --------------------------------------------------------------- */
/* ParabPeak ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Improve estimate of peak position (xpk, ypk) using fit of
// parabola to three points through peak. (d) is how many pixels
// away from peak to get samples.
//
// If, for XZ-plane, Z(x) = a(x - b)^2 + c, and,
//
// z0 = Z(xpk - d),
// z1 = Z(xpk),
// z2 = Z(xpk + d), then,
//
//        d      z0 - z2
// xpk' = - * -------------- + xpk
//        2    z0 + z2 - 2z1
//
// Uses standard C ordering.
//
// Note: Although solution is analytic, the quality of the result
// depends upon the quality of the data and wild displacements do
// occur. Therefore we limit displacements to some small multiple
// of the sampling distance.
//
static void ParabPeak(
	double			&xpk,
	double			&ypk,
	int				d,
	const double	*I,
	int				N )
{
	int		ix = (int)xpk,
			iy = (int)ypk;
	double	z1 = I[ix + N*iy],
			z0,
			z2;

	if( (z0 = I[ix-d + N*iy]) &&
		(z2 = I[ix+d + N*iy]) ) {

		z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

		if( fabs( z0 ) <= 8 * d )
			xpk += z0;
	}

	if( (z0 = I[ix + N*(iy-d)]) &&
		(z2 = I[ix + N*(iy+d)]) ) {

		z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

		if( fabs( z0 ) <= 8 * d )
			ypk += z0;
	}
}

/* --------------------------------------------------------------- */
/* CorrPatchesRQ ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return normalized cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// X-interval of offsets (lags) searched is [nnegx, nposx).
// Y-interval of offsets (lags) searched is [nnegy, nposy).
//
// For example, following the FFTSize() discussion, if the
// X-widths of the images are w1, w2; to search the largest
// sensible magnitudes for negative and positive offsets, set
// [w1-1, w2)...A negative magnitude >= w1 just pushes patch1
// all the way off patch2. For the positive case, the largest
// sensible value is similarly w2-1. Hovever, these parameters
// are used to allocate storage to hold that many offsets and
// zero is counted among the positives, so you must add one to
// get w2. Tighter search ranges are allowed. The minumum is
// [nnegx, nposx) = [1, 1).
//
// 'fft2' is a cache of the patch2 FFT. On entry, if fft2 has
// the correct size it is used. Otherwise recomputed here.
//
// Find peak using Q...correlation values are Pearson's R.
//
double CorrPatchesRQ(
	FILE					*flog,
	int						verbose,
	double					&Q,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 )
{
// Bounding boxes of point lists

	IBox	B1, B2;
	int		w1, h1, w2, h2;

	BBoxFromPoints( B1, ip1 );
	BBoxFromPoints( B2, ip2 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

	w2 = B2.R - B2.L + 1;
	h2 = B2.T - B2.B + 1;

	if( verbose ) {

		fprintf( flog,
		"RQ: region size is [%d %d] in x, [%d %d] in y.\n",
		B1.L, B1.R, B1.B, B1.T );

		fprintf( flog,
		"RQ: target size is [%d %d] in x, [%d %d] in y.\n",
		B2.L, B2.R, B2.B, B2.T );
	}

// Get array sizes (Nx,Ny) and FFT size M

	int	Nx	= FFTSize( w1, w2, nnegx, nposx ),
		Ny	= FFTSize( h1, h2, nnegy, nposy ),
		M	= Ny*(Nx/2+1);

	if( verbose )
		fprintf( flog, "RQ: Nx = %d, Ny = %d.\n", Nx, Ny );

// Create images from point lists.

	vector<double>	i1, i2;

	ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
	ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1;

	if( fft2.size() != M )
		FFT_2D( fft2, i2, Nx, Ny );

	FFT_2D( fft1, i1, Nx, Ny );

	for( int i = 0; i < M; ++i )
		fft1[i] = fft2[i] * conj( fft1[i] );

	IFT_2D( rslt, fft1, Nx, Ny );

// Prepare correlation calculator

	CRQ	ccalc;	// uses Pearson's r

	ccalc.Initialize( i1, w1, h1, i2, w2, h2, Nx, Ny );

// Reorganize valid entries of rslt image so that (dx,dy)=(0,0)
// is at the image center and (cx,cy)+(dx,dy) indexes all pixels
// of any sign.

	int	cx	= nnegx,
		cy	= nnegy,
		wR	= nnegx + nposx,
		hR	= nnegy + nposy,
		nR	= wR * hR;

	vector<double>	R( nR, 0.0 );

	for( int y = -nnegy; y < nposy; ++y ) {

		int	iy = Nx * (y >= 0 ? y : Ny + y);

		for( int x = -nnegx; x < nposx; ++x ) {

			if( !ccalc.CheckSize( LegalRgn, arglr, x, y ) )
				continue;

			if( !ccalc.CheckDensity( LegalCnt, arglc ) )
				continue;

			int	ix = (x >= 0 ? x : Nx + x);

			R[cx+x + wR*(cy+y)] = ccalc.CalcR( rslt[ix+iy] );
		}
	}

	OffsetResults( R );

//-----------------------------
#if	_debugcorr == 1
	CorrThmToTif8( "thmA.tif", i1, Nx, w1, h1 );
	CorrThmToTif8( "thmB.tif", i2, Nx, w2, h2 );
	vector<uint8> tif( nR );
	double mx = 0;
	for( int i = 0; i < nR; ++i ) if( R[i] > mx ) mx = R[i];
	for( int i = 0; i < nR; ++i ) {

		if( R[i] > 0.0 )
			tif[i] = int(250 * R[i] / mx);
		else
			tif[i] = 0;
	}
	Raster8ToTif8( "corr.tif", &tif[0], wR, hR );
	printf( "Center = (%d %d)\n", cx, cy );
#endif
//-----------------------------

// Scan R image for the best peak such that:
// - peak not adjacent to border.
// - peak not adjacent to any zero (4-way).
// - peak has maximal Q = [4*pk - SUM(4 neighbors)].
// The latter metric is average slope-like.

	int		qx	= -1,
			qy	= -1;

	Q = -1E30;

	for( int y = 1; y < hR - 1; ++y ) {

		for( int x = 1; x < wR - 1; ++x ) {

			double	q, t;
			int		iq = x + wR*y;

			if( (q = R[iq]) <= 0.0 )
				continue;

			q *= 4;

			if( !(t = R[iq-1]) )
				continue;

			q -= t;

			if( !(t = R[iq+1]) )
				continue;

			q -= t;

			if( !(t = R[iq-wR]) )
				continue;

			q -= t;

			if( !(t = R[iq+wR]) )
				continue;

			q -= t;

			if( q > Q ) {
				Q	= q;
				qx	= x;
				qy	= y;
			}
		}
	}

// Now for the conventional biggest

	double	bigR	= R[qx + wR*qy];
	int		bigx	= qx - cx;
	int		bigy	= qy - cy;

// Reports

	if( qx == -1 ) {

		fprintf( flog, "RQ: No legal subregions at all...\n");

		// for sake of grep searches

		fprintf( flog,
		"RQ: Max corr %f, max Q %f at [%d,%d].\n",
		0.0, 0.0, 0, 0 );

		Q	= 0.0;
		dx	= 0.0;
		dy	= 0.0;

		return 0.0;
	}

	if( verbose ) {

		fprintf( flog,
		"RQ: Max corr %f, max Q %f at [%d,%d] (%d,%d).\n",
		bigR, Q, bigx, bigy, qx, qy );
	}

// Interpolate peak

	dx = qx;
	dy = qy;
	ParabPeak( dx, dy, 1, &R[0], wR );

	dx += B2.L - B1.L - cx;
	dy += B2.B - B1.B - cy;

	return bigR;
}

/* --------------------------------------------------------------- */
/* CorrPatchesMaxQ ----------------------------------------------- */
/* --------------------------------------------------------------- */

#if	_debugcorr == 1
// Descending sort of these data
class SortQ{
public:
	double	q,  r;
	int		qx, qy;
};


static bool SortQ_q_dec( const SortQ &A, const SortQ &B )
{
	return A.q > B.q;
}


static bool SortQ_r_dec( const SortQ &A, const SortQ &B )
{
	return A.r > B.r;
}
#endif


// Return normalized cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// X-interval of offsets (lags) searched is [nnegx, nposx).
// Y-interval of offsets (lags) searched is [nnegy, nposy).
//
// For example, following the FFTSize() discussion, if the
// X-widths of the images are w1, w2; to search the largest
// sensible magnitudes for negative and positive offsets, set
// [w1-1, w2)...A negative magnitude >= w1 just pushes patch1
// all the way off patch2. For the positive case, the largest
// sensible value is similarly w2-1. Hovever, these parameters
// are used to allocate storage to hold that many offsets and
// zero is counted among the positives, so you must add one to
// get w2. Tighter search ranges are allowed. The minumum is
// [nnegx, nposx) = [1, 1).
//
// 'fft2' is a cache of the patch2 FFT. On entry, if fft2 has
// the correct size it is used. Otherwise recomputed here.
//
// Find peak using Q...correlation values are normed cross corrs.
//
double CorrPatchesMaxQ(
	FILE					*flog,
	int						verbose,
	double					&Q,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 )
{
// Bounding boxes of point lists

	IBox	B1, B2;
	int		w1, h1, w2, h2;

	BBoxFromPoints( B1, ip1 );
	BBoxFromPoints( B2, ip2 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

	w2 = B2.R - B2.L + 1;
	h2 = B2.T - B2.B + 1;

	if( verbose ) {

		fprintf( flog,
		"MaxQ: region size is [%d %d] in x, [%d %d] in y.\n",
		B1.L, B1.R, B1.B, B1.T );

		fprintf( flog,
		"MaxQ: target size is [%d %d] in x, [%d %d] in y.\n",
		B2.L, B2.R, B2.B, B2.T );
	}

// Get array sizes (Nx,Ny) and FFT size M

	int	Nx	= FFTSize( w1, w2, nnegx, nposx ),
		Ny	= FFTSize( h1, h2, nnegy, nposy ),
		M	= Ny*(Nx/2+1);

	if( verbose )
		fprintf( flog, "MaxQ: Nx = %d, Ny = %d.\n", Nx, Ny );

// Create images from point lists.

	vector<double>	i1, i2;

	ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
	ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1;

	if( fft2.size() != M )
		FFT_2D( fft2, i2, Nx, Ny );

	FFT_2D( fft1, i1, Nx, Ny );

	for( int i = 0; i < M; ++i )
		fft1[i] = fft2[i] * conj( fft1[i] );

	IFT_2D( rslt, fft1, Nx, Ny );

// Prepare correlation calculator

	CXCQ	ccalc;	// uses 1/n * SUM(a*b)

	ccalc.Initialize( i1, w1, h1, i2, w2, h2, Nx, Ny );

// Reorganize valid entries of rslt image so that (dx,dy)=(0,0)
// is at the image center and (cx,cy)+(dx,dy) indexes all pixels
// of any sign.

	int	cx	= nnegx,
		cy	= nnegy,
		wR	= nnegx + nposx,
		hR	= nnegy + nposy,
		nR	= wR * hR;

	vector<double>	R( nR, 0.0 );

	for( int y = -nnegy; y < nposy; ++y ) {

		int	iy = Nx * (y >= 0 ? y : Ny + y);

		for( int x = -nnegx; x < nposx; ++x ) {

			if( !ccalc.CheckSize( LegalRgn, arglr, x, y ) )
				continue;

			if( !ccalc.CheckDensity( LegalCnt, arglc ) )
				continue;

			int	ix = (x >= 0 ? x : Nx + x);

			R[cx+x + wR*(cy+y)] = rslt[ix+iy] / ccalc.GetNorm();
		}
	}

	OffsetResults( R );

//-----------------------------
#if	_debugcorr == 1
	CorrThmToTif8( "thmA.tif", i1, Nx, w1, h1 );
	CorrThmToTif8( "thmB.tif", i2, Nx, w2, h2 );
	vector<uint8> tif( nR );
	double mx = 0;
	for( int i = 0; i < nR; ++i ) if( R[i] > mx ) mx = R[i];
	for( int i = 0; i < nR; ++i ) {

		if( R[i] > 0.0 )
			tif[i] = int(250 * R[i] / mx);
		else
			tif[i] = 0;
	}
	Raster8ToTif8( "corr.tif", &tif[0], wR, hR );
	printf( "Center = (%d %d)\n", cx, cy );
#endif
//-----------------------------

// Scan R image for the best peak such that:
// - peak not adjacent to border.
// - peak not adjacent to any zero (4-way).
// - peak has maximal Q = [4*pk - SUM(4 neighbors)].
// The latter metric is average slope-like.

	int		qx	= -1,
			qy	= -1;

	Q = -1E30;

//-----------------------------
#if	_debugcorr == 1
	vector<SortQ> SQ;
#endif
//-----------------------------

	for( int y = 1; y < hR - 1; ++y ) {

		for( int x = 1; x < wR - 1; ++x ) {

			double	q, t;
			int		iq = x + wR*y;

			if( (q = R[iq]) <= 0.0 )
				continue;

			q *= 4;

			if( !(t = R[iq-1]) )
				continue;

			q -= t;

			if( !(t = R[iq+1]) )
				continue;

			q -= t;

			if( !(t = R[iq-wR]) )
				continue;

			q -= t;

			if( !(t = R[iq+wR]) )
				continue;

			q -= t;

//-----------------------------
#if	_debugcorr == 1
// look at top scoring points
			SortQ	sq;
			sq.q	= q;
			sq.r	= R[x + wR*y];
			sq.qx	= x;
			sq.qy	= y;
			SQ.push_back( sq );
#endif
//-----------------------------

			if( q > Q ) {
				Q	= q;
				qx	= x;
				qy	= y;
			}
		}
	}

//-----------------------------
#if	_debugcorr == 1
	int	nSQ = SQ.size();
	if( nSQ > 20 )
		nSQ = 20;

	sort( SQ.begin(), SQ.end(), SortQ_q_dec );
	printf( "top %d by Q:\n", nSQ );
	for( int i = 0; i < nSQ; ++i ) {
		printf( "Q %.3f R %.3f << %d, %d >>\n",
		SQ[i].q, SQ[i].r, SQ[i].qx, SQ[i].qy );
	}

	sort( SQ.begin(), SQ.end(), SortQ_r_dec );
	printf( "top %d by R:\n", nSQ );
	for( int i = 0; i < nSQ; ++i ) {
		printf( "Q %.3f R %.3f << %d, %d >>\n",
		SQ[i].q, SQ[i].r, SQ[i].qx, SQ[i].qy );
	}
#endif
//-----------------------------

// Now for the conventional biggest

	double	bigR	= R[qx + wR*qy];
	int		bigx	= qx - cx;
	int		bigy	= qy - cy;

// Reports

	if( qx == -1 ) {

		fprintf( flog, "MaxQ: No legal subregions at all...\n");

		// for sake of grep searches

		fprintf( flog,
		"MaxQ: Max corr %f, max Q %f at [%d,%d].\n",
		0.0, 0.0, 0, 0 );

		Q	= 0.0;
		dx	= 0.0;
		dy	= 0.0;

		return 0.0;
	}

	if( verbose ) {

		fprintf( flog,
		"MaxQ: Max corr %f, max Q %f at [%d,%d] (%d,%d).\n",
		bigR, Q, bigx, bigy, qx, qy );
	}

// Interpolate peak

	dx = qx;
	dy = qy;
	ParabPeak( dx, dy, 1, &R[0], wR );

	dx += B2.L - B1.L - cx;
	dy += B2.B - B1.B - cy;

	return bigR;
}

/* --------------------------------------------------------------- */
/* CorrPatchesMaxR ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Return normalized cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// X-interval of offsets (lags) searched is [nnegx, nposx).
// Y-interval of offsets (lags) searched is [nnegy, nposy).
//
// For example, following the FFTSize() discussion, if the
// X-widths of the images are w1, w2; to search the largest
// sensible magnitudes for negative and positive offsets, set
// [w1-1, w2)...A negative magnitude >= w1 just pushes patch1
// all the way off patch2. For the positive case, the largest
// sensible value is similarly w2-1. Hovever, these parameters
// are used to allocate storage to hold that many offsets and
// zero is counted among the positives, so you must add one to
// get w2. Tighter search ranges are allowed. The minumum is
// [nnegx, nposx) = [1, 1).
//
// 'fft2' is a cache of the patch2 FFT. On entry, if fft2 has
// the correct size it is used. Otherwise recomputed here.
//
// Find peak using R...correlation values are normed cross corrs.
//
double CorrPatchesMaxR(
	FILE					*flog,
	int						verbose,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 )
{
// Bounding boxes of point lists

	IBox	B1, B2;
	int		w1, h1, w2, h2;

	BBoxFromPoints( B1, ip1 );
	BBoxFromPoints( B2, ip2 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

	w2 = B2.R - B2.L + 1;
	h2 = B2.T - B2.B + 1;

	if( verbose ) {

		fprintf( flog,
		"MaxR: region size is [%d %d] in x, [%d %d] in y.\n",
		B1.L, B1.R, B1.B, B1.T );

		fprintf( flog,
		"MaxR: target size is [%d %d] in x, [%d %d] in y.\n",
		B2.L, B2.R, B2.B, B2.T );
	}

// Get array sizes (Nx,Ny) and FFT size M

	int	Nx	= FFTSize( w1, w2, nnegx, nposx ),
		Ny	= FFTSize( h1, h2, nnegy, nposy ),
		M	= Ny*(Nx/2+1);

	if( verbose )
		fprintf( flog, "MaxR: Nx = %d, Ny = %d.\n", Nx, Ny );

// Create images from point lists.

	vector<double>	i1, i2;

	ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
	ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1;

	if( fft2.size() != M )
		FFT_2D( fft2, i2, Nx, Ny );

	FFT_2D( fft1, i1, Nx, Ny );

	for( int i = 0; i < M; ++i )
		fft1[i] = fft2[i] * conj( fft1[i] );

	IFT_2D( rslt, fft1, Nx, Ny );

// Prepare correlation calculator

	CXCQ	ccalc;	// uses 1/n * SUM(a*b)

	ccalc.Initialize( i1, w1, h1, i2, w2, h2, Nx, Ny );

// Reorganize valid entries of rslt image so that (dx,dy)=(0,0)
// is at the image center and (cx,cy)+(dx,dy) indexes all pixels
// of any sign.

	int	cx	= nnegx,
		cy	= nnegy,
		wR	= nnegx + nposx,
		hR	= nnegy + nposy,
		nR	= wR * hR;

	vector<double>	R( nR, 0.0 );

	for( int y = -nnegy; y < nposy; ++y ) {

		int	iy = Nx * (y >= 0 ? y : Ny + y);

		for( int x = -nnegx; x < nposx; ++x ) {

			if( !ccalc.CheckSize( LegalRgn, arglr, x, y ) )
				continue;

			if( !ccalc.CheckDensity( LegalCnt, arglc ) )
				continue;

			int	ix = (x >= 0 ? x : Nx + x);

			R[cx+x + wR*(cy+y)] = rslt[ix+iy] / ccalc.GetNorm();
		}
	}

	OffsetResults( R );

//-----------------------------
#if	_debugcorr == 1
	CorrThmToTif8( "thmA.tif", i1, Nx, w1, h1 );
	CorrThmToTif8( "thmB.tif", i2, Nx, w2, h2 );
	vector<uint8> tif( nR );
	double mx = 0;
	for( int i = 0; i < nR; ++i ) if( R[i] > mx ) mx = R[i];
	for( int i = 0; i < nR; ++i ) {

		if( R[i] > 0.0 )
			tif[i] = int(250 * R[i] / mx);
		else
			tif[i] = 0;
	}
	Raster8ToTif8( "corr.tif", &tif[0], wR, hR );
	printf( "Center = (%d %d)\n", cx, cy );
#endif
//-----------------------------

// Scan R image for the best peak such that:
// - peak not adjacent to border.
// - peak not adjacent to any zero (4-way).
// - peak indeed higher than neighbors (4-way).

	double	bigR	= 0.0;
	int		rx		= -1;
	int		ry		= -1;

	for( int y = 1; y < hR - 1; ++y ) {

		for( int x = 1; x < wR - 1; ++x ) {

			double	r, t;
			int		ir = x + wR*y;

			if( (r = R[ir]) < bigR )
				continue;

			if( !(t = R[ir-1]) || t >= r )
				continue;

			if( !(t = R[ir+1]) || t >= r )
				continue;

			if( !(t = R[ir-wR]) || t >= r )
				continue;

			if( !(t = R[ir+wR]) || t >= r )
				continue;

			bigR	= r;
			rx		= x;
			ry		= y;
		}
	}

// Reports

	if( rx == -1 ) {

		fprintf( flog, "MaxR: No legal subregions at all...\n");

		// for sake of grep searches

		fprintf( flog,
		"MaxR: Max corr %f at [%d,%d].\n",
		0.0, 0, 0 );

		dx	= 0.0;
		dy	= 0.0;

		return 0.0;
	}

	int	bigx = rx - cx;
	int	bigy = ry - cy;

	if( verbose ) {

		fprintf( flog,
		"MaxR: Max corr %f at [%d,%d] (%d,%d).\n",
		bigR, bigx, bigy, rx, ry );
	}

// Interpolate peak

	dx = rx;
	dy = ry;
	ParabPeak( dx, dy, 1, &R[0], wR );

	dx += B2.L - B1.L - cx;
	dy += B2.B - B1.B - cy;

	return bigR;
}


