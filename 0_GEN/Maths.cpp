

#include	"Maths.h"

#include	<math.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	BINVAL( i )	(datamin + (i) * bwid)






/* --------------------------------------------------------------- */
/* Numerics ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Given n >= 0; return p such that: p=2^k && p >= n.
//
int CeilPow2( int n )
{
	int	p = 1;

	while( p < n )
		p *= 2;

	return p;
}

/* --------------------------------------------------------------- */
/* MeanStd::Stats ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MeanStd::Stats( double &avg, double &std )
{
	avg	= sum;
	std	= 0.0;

	if( n > 1 ) {
		avg /= n;
#ifdef TINYSTAT
		std  = sqrt( (sum2 - n*avg*avg) / n );
#else
		std  = sqrt( (sum2 - n*avg*avg) / (n - 1.0) );
#endif
	}
}

/* --------------------------------------------------------------- */
/* Stats --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Finds statistics on a vector of doubles
//
void Stats( const vector<double> &v, double &avg, double &std )
{
	MeanStd	m;
	int		n = v.size();

	for( int i = 0; i < n; ++i )
		m.Element( v[i] );

	m.Stats( avg, std );
}

// Finds statistics on an array of N uint8
//
void Stats( const uint8 *v, int N, double &avg, double &std )
{
	MeanStd	m;

	for( int i = 0; i < N; ++i )
		m.Element( v[i] );

	m.Stats( avg, std );
}

/* --------------------------------------------------------------- */
/* StatsRasterNonZeros ------------------------------------------- */
/* --------------------------------------------------------------- */

// Finds statistics on uint8 raster
//
void StatsRasterNonZeros(
	const uint8*	raster,
	int				npix,
	double			&avg,
	double			&std )
{
	MeanStd	m;

	for( int i = 0; i < npix; ++i ) {

		if( raster[i] )
			m.Element( raster[i] );
	}

	m.Stats( avg, std );

	printf( "RasterStats: %d nonzero, %f percent.\n",
		m.HowMany(), m.HowMany() * 100.0 / npix );

	printf( "RasterStats: For nonzeros: mean = %f, std dev = %f,\n",
		avg, std );
}

/* --------------------------------------------------------------- */
/* Normalize ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Normalizes vector to mean 0 and standard deviation 1
//
void Normalize( vector<double> &v )
{
	double	avg, std;
	int		n = v.size();

	Stats( v, avg, std );

	//printf( "Normalize: Average=%f, std dev=%f.\n", avg, std );

	for( int i = 0; i < n; ++i )
		v[i] = (v[i] - avg) / std;
}


// Normalizes doubles array to mean 0 and standard deviation 1
//
void Normalize( double *a, int n )
{
	MeanStd	m;
	double	avg, std;
	int		i;

	for( i = 0; i < n; ++i )
		m.Element( a[i] );

	m.Stats( avg, std );

	//printf( "Normalize: Average=%f, std dev=%f.\n", avg, std );

	for( i = 0; i < n; ++i )
		a[i] = (a[i] - avg) / std;
}

/* --------------------------------------------------------------- */
/* NormalizeNonZeros --------------------------------------------- */
/* --------------------------------------------------------------- */

// Normalizes vector to mean 0 and standard deviation 1,
// but only the non-zero elements.
//
void NormalizeNonZeros( vector<double> &v )
{
	MeanStd	m;
	double	avg, std;
	int		n = v.size();

	for( int i = 0; i < n; ++i ) {

		if( v[i] )
			m.Element( v[i] );
	}

	m.Stats( avg, std );

	for( int i = 0; i < n; ++i ) {

		if( v[i] )
			v[i] = (v[i] - avg) / std;
	}
}

/* --------------------------------------------------------------- */
/* CoNormalize --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Normalizes two vectors to common scale: mn,std = 0,1
//
void CoNormalize( vector<double> &v1, vector<double> &v2 )
{
	MeanStd	m;
	double	avg, std;
	int		i, n1 = v1.size(), n2 = v2.size();

	for( i = 0; i < n1; ++i )
		m.Element( v1[i] );

	for( i = 0; i < n2; ++i )
		m.Element( v2[i] );

	m.Stats( avg, std );

	//printf( "CoNormalize: Average=%f, std dev=%f.\n", avg, std );

	for( i = 0; i < n1; ++i )
		v1[i] = (v1[i] - avg) / std;

	for( i = 0; i < n2; ++i )
		v2[i] = (v2[i] - avg) / std;
}

/* --------------------------------------------------------------- */
/* CoUpperTail --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Like CoNormalize, but all values below mean + nsd sigma
// are zeroed.
//
// Recommend calling CoNormalize again afterward.
//
void CoUpperTail(
	vector<double>	&v1,
	vector<double>	&v2,
	double			nsigma )
{
	MeanStd	m;
	double	avg, std;
	int		i, n1 = v1.size(), n2 = v2.size();

	for( i = 0; i < n1; ++i )
		m.Element( v1[i] );

	for( i = 0; i < n2; ++i )
		m.Element( v2[i] );

	m.Stats( avg, std );
	nsigma = avg + std * nsigma;

	for( i = 0; i < n1; ++i ) {

		if( v1[i] >= nsigma )
			v1[i] = (v1[i] - avg) / std;
		else
			v1[i] = 0.0;
	}

	for( i = 0; i < n2; ++i ) {

		if( v2[i] >= nsigma )
			v2[i] = (v2[i] - avg) / std;
		else
			v2[i] = 0.0;
	}
}

/* --------------------------------------------------------------- */
/* CoExcludeMiddle ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Like CoNormalize, but all values within nsd sigma of mean
// are zeroed.
//
// Recommend calling CoNormalize again afterward.
//
void CoExcludeMiddle(
	vector<double>	&v1,
	vector<double>	&v2,
	double			nsigma )
{
	MeanStd	m;
	double	avg, std;
	int		i, n1 = v1.size(), n2 = v2.size();

	for( i = 0; i < n1; ++i )
		m.Element( v1[i] );

	for( i = 0; i < n2; ++i )
		m.Element( v2[i] );

	m.Stats( avg, std );
	nsigma *= std;

	for( i = 0; i < n1; ++i ) {

		if( fabs( v1[i] - avg ) >= nsigma )
			v1[i] = (v1[i] - avg) / std;
		else
			v1[i] = 0.0;
	}

	for( i = 0; i < n2; ++i ) {

		if( fabs( v2[i] - avg ) >= nsigma )
			v2[i] = (v2[i] - avg) / std;
		else
			v2[i] = 0.0;
	}
}

/* --------------------------------------------------------------- */
/* FirstNonzero -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Returns -1 if all zero, else index of first non-zero value.
//
int FirstNonzero( const double *v, int n )
{
	int	idx = -1;

	for( int i = 0; i < n; ++i ) {

		if( v[i] ) {
			idx	= i;
			break;
		}
	}

	return idx;
}

/* --------------------------------------------------------------- */
/* IndexOfMaxVal ------------------------------------------------- */
/* --------------------------------------------------------------- */

int IndexOfMaxVal( const double *v, int n )
{
	double	max = v[0];
	int		idx = 0;

	for( int i = 1; i < n; ++i ) {

		if( v[i] > max ) {
			max	= v[i];
			idx	= i;
		}
	}

	return idx;
}

/* --------------------------------------------------------------- */
/* MedianVal ----------------------------------------------------- */
/* --------------------------------------------------------------- */

double MedianVal( vector<double> &V )
{
	int	n = V.size();

	if( !n )
		return 0.0;

	sort( V.begin(), V.end() );

	if( n & 1 )
		return V[n/2];
	else
		return (V[n/2 - 1] + V[n/2]) / 2.0;
}

/* --------------------------------------------------------------- */
/* LineFit ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fit data pairs (x,y) to line on interval [i0,ilim).
//
void LineFit(
	double			*icpt,
	double			*slope,
	double			*lincor,
	const double	*x,
	const double	*y,
	int				i0,
	int				ilim )
{
	double	r, varX, varY,
			n		= 0.0,
			sumX	= 0.0,
			sumY	= 0.0,
			sumXX	= 0.0,
			sumYY	= 0.0,
			sumXY	= 0.0;

/* --------------------------- */
/* Collect sums over [i0,ilim) */
/* --------------------------- */

	for( int i = i0; i < ilim; ++i ) {

		double	X = x[i],
				Y = y[i];

		sumX	+= X;
		sumY	+= Y;
		sumXX	+= X * X;
		sumYY	+= Y * Y;
		sumXY	+= X * Y;
		++n;
	}

/* --------------- */
/* Compute results */
/* --------------- */

	r		= n * sumXY - sumX * sumY;
	varX	= n * sumXX - sumX * sumX;
	varY	= n * sumYY - sumY * sumY;

	if( icpt )
		*icpt = (sumXX * sumY - sumX * sumXY) / varX;

	if( slope )
		*slope = r / varX;

	if( lincor )
		*lincor = r / sqrt( varX * varY );
}

/* --------------------------------------------------------------- */
/* Histogram ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void Histogram(
	double			&uflo,
	double			&oflo,
	double			*bins,
	int				nbins,
	double			datamin,
	double			datamax,
	const uint16	*data,
	int				ndata,
	bool			breset )
{
	double	dtob = nbins / (datamax - datamin);

	if( breset ) {
		uflo	= 0.0;
		oflo	= 0.0;
		memset( bins, 0, nbins * sizeof(double) );
	}

	for( int i = 0; i < ndata; ++i ) {

		if( data[i] < datamin )
			++uflo;
		else if( data[i] >= datamax )
			++oflo;
		else
			++bins[int((data[i] - datamin)*dtob)];
	}
}

/* --------------------------------------------------------------- */
/* PercentileBin ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return bin index I such that integrated
// count over [0..I] is >= frac * total.
//
int PercentileBin(
	const double	*bins,
	int				nbins,
	double			frac )
{
	double	cnts = 0.0, sum = 0.0;

	for( int i = 0; i < nbins; ++i )
		cnts += bins[i];

	cnts *= frac;

	for( int i = 0; i < nbins; ++i ) {

		sum += bins[i];

		if( sum >= cnts )
			return i;
	}

	return nbins - 1;
}

/* --------------------------------------------------------------- */
/* IsoDataThresh ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return threshold value that lies midway between the average
// of all binned data below and the average of all binned data
// above that threshold.
//
// Returns -1 if histogram empty.
//
double IsoDataThresh(
	const double	*bins,
	int				nbins,
	double			datamin,
	double			datamax )
{
	double	bwid = (datamax - datamin) / nbins;
	double	nlo, dlo, nhi, dhi, T, err;
	int		ilo, ihi;

	ilo = FirstNonzero( bins, nbins );
	ihi = -1;

// find lowest non-zero bin

	if( ilo < 0 )
		return -1.0;

// find highest non-zero bin

	for( int i = nbins - 1; i >= ilo; --i ) {

		if( bins[i] ) {
			ihi = i;
			break;
		}
	}

	if( ihi == ilo )
		return BINVAL( ilo );

// start with all data on high side

	dlo	= 0.0;
	nlo	= 0.0;
	dhi	= 0.0;
	nhi	= 0.0;
	err = 1e30;

	for( int i = ilo; i <= ihi; ++i ) {
		dhi += bins[i];
		nhi += bins[i] * BINVAL( i );
	}

// for each thresh slide one bin to left

	for( int i = ilo; i < ihi; ++i ) {

		double	q, t;

		q    = bins[i];
		dlo += q;
		dhi -= q;

		q   *= BINVAL( i );
		nlo += q;
		nhi -= q;

		t = 0.5*(nlo/dlo + nhi/dhi);
		q = fabs( BINVAL( i + 1 ) - t );

		if( q < err ) {
			T	= t;
			err	= q;
		}
	}

	return T;
}

/* --------------------------------------------------------------- */
/* MinSepThresh -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return threshold value minimizing the distance between the
// averages of all binned data below and above that threshold.
//
// Returns -1 if histogram empty.
//
double MinSepThresh(
	const double	*bins,
	int				nbins,
	double			datamin,
	double			datamax )
{
	double	bwid = (datamax - datamin) / nbins;
	double	nlo, dlo, nhi, dhi, T, sep;
	int		ilo, ihi;

	ilo = FirstNonzero( bins, nbins );
	ihi = -1;

// find lowest non-zero bin

	if( ilo < 0 )
		return -1.0;

// find highest non-zero bin

	for( int i = nbins - 1; i >= ilo; --i ) {

		if( bins[i] ) {
			ihi = i;
			break;
		}
	}

	if( ihi == ilo )
		return BINVAL( ilo );

// start with all data on high side

	dlo	= 0.0;
	nlo	= 0.0;
	dhi	= 0.0;
	nhi	= 0.0;
	sep = 1e30;

	for( int i = ilo; i <= ihi; ++i ) {
		dhi += bins[i];
		nhi += bins[i] * BINVAL( i );
	}

// for each thresh slide one bin to left

	for( int i = ilo; i < ihi; ++i ) {

		double	q, t;

		q    = bins[i];
		dlo += q;
		dhi -= q;

		q   *= BINVAL( i );
		nlo += q;
		nhi -= q;

		q = fabs( nlo/dlo - nhi/dhi );

		if( q < sep ) {
			T	= nlo/dlo + nhi/dhi;
			sep	= q;
		}
	}

	return T * 0.5;
}

/* --------------------------------------------------------------- */
/* OtsuThresh ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return threshold value minimizing the summed variances about
// the means of data below and above that threshold.
//
// Returns -1 if histogram empty.
//
double OtsuThresh(
	const double	*bins,
	int				nbins,
	double			datamin,
	double			datamax )
{
	double	bwid = (datamax - datamin) / nbins;
	double	tlo, slo, nlo, thi, shi, nhi, T, var;
	int		ilo, ihi;

	ilo = FirstNonzero( bins, nbins );
	ihi = -1;

// find lowest non-zero bin

	if( ilo < 0 )
		return -1.0;

// find highest non-zero bin

	for( int i = nbins - 1; i >= ilo; --i ) {

		if( bins[i] ) {
			ihi = i;
			break;
		}
	}

	if( ihi == ilo )
		return BINVAL( ilo );

// start with all data on high side

	nlo	= 0.0;
	slo	= 0.0;
	nhi	= 0.0;
	shi	= 0.0;
	var	= 0;

	for( int i = ilo; i <= ihi; ++i ) {
		nhi	+= bins[i];
		shi	+= bins[i] * BINVAL( i );
	}

// for each thresh slide one bin to left

	for( int i = ilo; i < ihi; ++i ) {

		double	q;

		q	 = bins[i];
		nlo	+= q;
		nhi	-= q;

		q	*= BINVAL( i );
		slo	+= q;
		shi	-= q;

		tlo	 = slo / nlo;
		thi	 = shi / nhi;

		q	 = tlo - thi;
		q	*= nlo * nhi * q;

		if( q > var ) {
			T	= tlo + thi;
			var	= q;
		}
	}

	return T * 0.5;
}

/* --------------------------------------------------------------- */
/* PrintVectorAsMat ---------------------------------------------- */
/* --------------------------------------------------------------- */

void PrintVectorAsMat( FILE* f, const vector<double> &v, int w )
{
	int		n = v.size();

	for( int i = 0; i < n; ++i ) {

		if( !(i % w) )
			fprintf( f, "\n" );

		fprintf( f, "%5.1f", v[i] );
	}

	fprintf( f, "\n" );
}

/* --------------------------------------------------------------- */
/* Print3x3Matrix ------------------------------------------------ */
/* --------------------------------------------------------------- */

void Print3x3Matrix( FILE* f, const double a[3][3] )
{
	for( int i = 0; i < 3; ++i ) {

		for( int j = 0; j < 3; ++j )
			fprintf( f, "%12.6f ", a[i][j] );

		fprintf( f, "\n" );
	}
}

/* --------------------------------------------------------------- */
/* Invert3x3Matrix ----------------------------------------------- */
/* --------------------------------------------------------------- */

double Invert3x3Matrix( double i[3][3], const double a[3][3] )
{
	double	det = a[0][0]*a[1][1]*a[2][2]
				+ a[0][1]*a[1][2]*a[2][0]
				+ a[1][0]*a[2][1]*a[0][2]
				- a[0][2]*a[1][1]*a[2][0]
				- a[0][1]*a[1][0]*a[2][2]
				- a[1][2]*a[2][1]*a[0][0];

	i[0][0] =  (a[1][1]*a[2][2] - a[1][2]*a[2][1]) / det;
	i[1][0] = -(a[1][0]*a[2][2] - a[1][2]*a[2][0]) / det;
	i[2][0] =  (a[1][0]*a[2][1] - a[1][1]*a[2][0]) / det;
	i[0][1] = -(a[0][1]*a[2][2] - a[0][2]*a[2][1]) / det;
	i[1][1] =  (a[0][0]*a[2][2] - a[0][2]*a[2][0]) / det;
	i[2][1] = -(a[0][0]*a[2][1] - a[0][1]*a[2][0]) / det;
	i[0][2] =  (a[0][1]*a[1][2] - a[0][2]*a[1][1]) / det;
	i[1][2] = -(a[0][0]*a[1][2] - a[0][2]*a[1][0]) / det;
	i[2][2] =  (a[0][0]*a[1][1] - a[0][1]*a[1][0]) / det;

	return det;
}

/* --------------------------------------------------------------- */
/* Invert4x4Matrix ----------------------------------------------- */
/* --------------------------------------------------------------- */

double Invert4x4Matrix( double i[4][4], const double a[4][4] )
{
	double	det = a[0][3]*a[1][2]*a[2][1]*a[3][0]
				- a[0][2]*a[1][3]*a[2][1]*a[3][0]
				- a[0][3]*a[1][1]*a[2][2]*a[3][0]
				+ a[0][1]*a[1][3]*a[2][2]*a[3][0]
				+ a[0][2]*a[1][1]*a[2][3]*a[3][0]
				- a[0][1]*a[1][2]*a[2][3]*a[3][0]
				- a[0][3]*a[1][2]*a[2][0]*a[3][1]
				+ a[0][2]*a[1][3]*a[2][0]*a[3][1]
				+ a[0][3]*a[1][0]*a[2][2]*a[3][1]
				- a[0][0]*a[1][3]*a[2][2]*a[3][1]
				- a[0][2]*a[1][0]*a[2][3]*a[3][1]
				+ a[0][0]*a[1][2]*a[2][3]*a[3][1]
				+ a[0][3]*a[1][1]*a[2][0]*a[3][2]
				- a[0][1]*a[1][3]*a[2][0]*a[3][2]
				- a[0][3]*a[1][0]*a[2][1]*a[3][2]
				+ a[0][0]*a[1][3]*a[2][1]*a[3][2]
				+ a[0][1]*a[1][0]*a[2][3]*a[3][2]
				- a[0][0]*a[1][1]*a[2][3]*a[3][2]
				- a[0][2]*a[1][1]*a[2][0]*a[3][3]
				+ a[0][1]*a[1][2]*a[2][0]*a[3][3]
				+ a[0][2]*a[1][0]*a[2][1]*a[3][3]
				- a[0][0]*a[1][2]*a[2][1]*a[3][3]
				- a[0][1]*a[1][0]*a[2][2]*a[3][3]
				+ a[0][0]*a[1][1]*a[2][2]*a[3][3];

	i[0][0] = (-a[1][3]*a[2][2]*a[3][1] + a[1][2]*a[2][3]*a[3][1] + a[1][3]*a[2][1]*a[3][2] - a[1][1]*a[2][3]*a[3][2] - a[1][2]*a[2][1]*a[3][3] + a[1][1]*a[2][2]*a[3][3]) / det;
	i[0][1] = ( a[0][3]*a[2][2]*a[3][1] - a[0][2]*a[2][3]*a[3][1] - a[0][3]*a[2][1]*a[3][2] + a[0][1]*a[2][3]*a[3][2] + a[0][2]*a[2][1]*a[3][3] - a[0][1]*a[2][2]*a[3][3]) / det;
	i[0][2] = (-a[0][3]*a[1][2]*a[3][1] + a[0][2]*a[1][3]*a[3][1] + a[0][3]*a[1][1]*a[3][2] - a[0][1]*a[1][3]*a[3][2] - a[0][2]*a[1][1]*a[3][3] + a[0][1]*a[1][2]*a[3][3]) / det;
	i[0][3] = ( a[0][3]*a[1][2]*a[2][1] - a[0][2]*a[1][3]*a[2][1] - a[0][3]*a[1][1]*a[2][2] + a[0][1]*a[1][3]*a[2][2] + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3]) / det;

	i[1][0] = ( a[1][3]*a[2][2]*a[3][0] - a[1][2]*a[2][3]*a[3][0] - a[1][3]*a[2][0]*a[3][2] + a[1][0]*a[2][3]*a[3][2] + a[1][2]*a[2][0]*a[3][3] - a[1][0]*a[2][2]*a[3][3]) / det;
	i[1][1] = (-a[0][3]*a[2][2]*a[3][0] + a[0][2]*a[2][3]*a[3][0] + a[0][3]*a[2][0]*a[3][2] - a[0][0]*a[2][3]*a[3][2] - a[0][2]*a[2][0]*a[3][3] + a[0][0]*a[2][2]*a[3][3]) / det;
	i[1][2] = ( a[0][3]*a[1][2]*a[3][0] - a[0][2]*a[1][3]*a[3][0] - a[0][3]*a[1][0]*a[3][2] + a[0][0]*a[1][3]*a[3][2] + a[0][2]*a[1][0]*a[3][3] - a[0][0]*a[1][2]*a[3][3]) / det;
	i[1][3] = (-a[0][3]*a[1][2]*a[2][0] + a[0][2]*a[1][3]*a[2][0] + a[0][3]*a[1][0]*a[2][2] - a[0][0]*a[1][3]*a[2][2] - a[0][2]*a[1][0]*a[2][3] + a[0][0]*a[1][2]*a[2][3]) / det;

	i[2][0] = (-a[1][3]*a[2][1]*a[3][0] + a[1][1]*a[2][3]*a[3][0] + a[1][3]*a[2][0]*a[3][1] - a[1][0]*a[2][3]*a[3][1] - a[1][1]*a[2][0]*a[3][3] + a[1][0]*a[2][1]*a[3][3]) / det;
	i[2][1] = ( a[0][3]*a[2][1]*a[3][0] - a[0][1]*a[2][3]*a[3][0] - a[0][3]*a[2][0]*a[3][1] + a[0][0]*a[2][3]*a[3][1] + a[0][1]*a[2][0]*a[3][3] - a[0][0]*a[2][1]*a[3][3]) / det;
	i[2][2] = (-a[0][3]*a[1][1]*a[3][0] + a[0][1]*a[1][3]*a[3][0] + a[0][3]*a[1][0]*a[3][1] - a[0][0]*a[1][3]*a[3][1] - a[0][1]*a[1][0]*a[3][3] + a[0][0]*a[1][1]*a[3][3]) / det;
	i[2][3] = ( a[0][3]*a[1][1]*a[2][0] - a[0][1]*a[1][3]*a[2][0] - a[0][3]*a[1][0]*a[2][1] + a[0][0]*a[1][3]*a[2][1] + a[0][1]*a[1][0]*a[2][3] - a[0][0]*a[1][1]*a[2][3]) / det;

	i[3][0] = ( a[1][2]*a[2][1]*a[3][0] - a[1][1]*a[2][2]*a[3][0] - a[1][2]*a[2][0]*a[3][1] + a[1][0]*a[2][2]*a[3][1] + a[1][1]*a[2][0]*a[3][2] - a[1][0]*a[2][1]*a[3][2]) / det;
	i[3][1] = (-a[0][2]*a[2][1]*a[3][0] + a[0][1]*a[2][2]*a[3][0] + a[0][2]*a[2][0]*a[3][1] - a[0][0]*a[2][2]*a[3][1] - a[0][1]*a[2][0]*a[3][2] + a[0][0]*a[2][1]*a[3][2]) / det;
	i[3][2] = ( a[0][2]*a[1][1]*a[3][0] - a[0][1]*a[1][2]*a[3][0] - a[0][2]*a[1][0]*a[3][1] + a[0][0]*a[1][2]*a[3][1] + a[0][1]*a[1][0]*a[3][2] - a[0][0]*a[1][1]*a[3][2]) / det;
	i[3][3] = (-a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]) / det;

	return det;
}

/* --------------------------------------------------------------- */
/* LegPolyCreate ------------------------------------------------- */
/* --------------------------------------------------------------- */

void LegPolyCreate(
	vector<vector<double> >	&L,
	int						maxOrder,
	int						nSamples )
{
	double	n = nSamples, n_1 = n - 1;
	int		io, is;

	L.resize( maxOrder + 1 );

	for( io = 0; io <= maxOrder; ++io )
		L[io].resize( nSamples );

	for( is = 0; is < nSamples; ++is ) {

		double	x = (2.0 * is - n_1) / n;	// x in range (-1,1)

		switch( maxOrder ) {

			default:
				L[6][is] = 0.0625 * (231.0*x*x*x*x*x*x - 315.0*x*x*x*x + 105.0*x*x - 5.0);
			case 5:
				L[5][is] = 0.125 * (63.0*x*x*x*x*x - 70.0*x*x*x + 15.0*x);
			case 4:
				L[4][is] = 0.125 * (35.0*x*x*x*x - 30.0*x*x + 3.0);
			case 3:
				L[3][is] = 2.5*x*x*x - 1.5*x;
			case 2:
				L[2][is] = 1.5*x*x - 0.5;
			case 1:
				L[1][is] = x;
			case 0:
				L[0][is] = 1.0;
		}
	}
}

/* --------------------------------------------------------------- */
/* LegPolyFlatten ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Project all src pixels onto low order Legendre polynomials
// and remove those components. Place results in vals.
//
void LegPolyFlatten(
	vector<double>		&vals,
	const uint8*		src,
	int					w,
	int					h,
	int					maxOrder )
{
// Copy indicated points to vals

	int		npts = w * h;

	vals.resize( npts );

	for( int i = 0; i < npts; ++i )
		vals[i] = src[i];

// Correct vals

	if( maxOrder > 0 ) {

		// Create Legendre basis polynomials

		vector<vector<double> > xpolys;
		vector<vector<double> > ypolys;

		LegPolyCreate( xpolys, maxOrder, w );
		LegPolyCreate( ypolys, maxOrder, h );

		// Remove projections onto polys up to maxOrder

		for( int xo = 0; xo <= maxOrder; ++xo ) {

			const double*	XP = &xpolys[xo][0];

			for( int yo = 0; yo <= maxOrder; ++yo ) {

				// the (0,0) term is just normalization
				// which we do later anyway

				if( xo + yo == 0 )
					continue;

				const double*	YP = &ypolys[yo][0];

				// get projection coefficient

				double	coef = 0.0, norm = 0.0;

				for( int i = 0; i < npts; ++i ) {

					int		y = i / w,
							x = i - w * y;
					double	F = XP[x] * YP[y];

					coef	+= F * vals[i];
					norm	+= F * F;
				}

				coef /= norm;

				//printf( "LegFlatten: xo=%d yo=%d cof=%f.\n",
				//	xo, yo, coef );

				// now remove that component

				for( int i = 0; i < npts; ++i ) {

					int		y = i / w,
							x = i - w * y;

					vals[i] -= coef * XP[x] * YP[y];
				}
			}
		}
	}

// Final normalization

	Normalize( vals );
}

/* --------------------------------------------------------------- */
/* LegPolyFlatten ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Project src pixels (those in pts list) onto low order Legendre
// polynomials and remove those components. Place results in vals.
//
void LegPolyFlatten(
	vector<double>		&vals,
	const vector<Point>	&pts,
	const uint8*		src,
	int					w,
	int					h,
	int					maxOrder )
{
// Copy indicated points to vals

	int		npts = pts.size();

	vals.resize( npts );

	for( int i = 0; i < npts; ++i ) {

		int	x = (int)pts[i].x,
			y = (int)pts[i].y;

		vals[i] = src[x + w * y];
	}

// Correct vals

	if( maxOrder > 0 ) {

		// Create Legendre basis polynomials

		vector<vector<double> > xpolys;
		vector<vector<double> > ypolys;

		LegPolyCreate( xpolys, maxOrder, w );
		LegPolyCreate( ypolys, maxOrder, h );

		// Remove projections onto polys up to maxOrder

		for( int xo = 0; xo <= maxOrder; ++xo ) {

			const double*	XP = &xpolys[xo][0];

			for( int yo = 0; yo <= maxOrder; ++yo ) {

				// the (0,0) term is just normalization
				// which we do later anyway

				if( xo + yo == 0 )
					continue;

				const double*	YP = &ypolys[yo][0];

				// get projection coefficient

				double	coef = 0.0, norm = 0.0;

				for( int i = 0; i < npts; ++i ) {

					double F =
					XP[(int)pts[i].x] * YP[(int)pts[i].y];

					coef	+= F * vals[i];
					norm	+= F * F;
				}

				coef /= norm;

				//printf( "LegFlatten: xo=%d yo=%d cof=%f.\n",
				//	xo, yo, coef );

				// now remove that component

				for( int i = 0; i < npts; ++i ) {

					vals[i] -= coef *
					XP[(int)pts[i].x] * YP[(int)pts[i].y];
				}
			}
		}
	}

// Final normalization

	Normalize( vals );
}

/* --------------------------------------------------------------- */
/* CopyRaster ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CopyRaster(
	double*			D,
	int				wD,
	const double*	S,
	int				wS,
	int				wCopy,
	int				hCopy )
{
	int		N = wCopy * sizeof(double);

	for( int i = 0; i < hCopy; ++i, D += wD, S += wS )
		memcpy( D, S, N );
}

/* --------------------------------------------------------------- */
/* DecimateVector ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Bin point and value lists by linear dimension lbin.
//
// IMPORTANT: On entry and exit, corner of point list's
// bounding box is (0,0).
//
void DecimateVector(
	vector<Point>	&p,
	vector<double>	&v,
	int				w,
	int				h,
	int				lbin )
{
// Sum vector data into value and count rasters
// New dims set to lowest multiple of lbin.

	w = (w + lbin - 1) / lbin;
	h = (h + lbin - 1) / lbin;

	int				nr = w * h,
					np = p.size();
	vector<double>	vr( nr, 0.0 );
	vector<uint8>	cr( nr, 0 );

	for( int ip = 0; ip < np; ++ip ) {

		int	x	= int(p[ip].x) / lbin,
			y	= int(p[ip].y) / lbin,
			ir	= x + w*y;

		vr[ir] += v[ip];
		++cr[ir];
	}

// Gather and average raster points having non-zero counts

	int		k = 0;

	for( int ir = 0; ir < nr; ++ir ) {

		int		cnt;

		if( cnt = cr[ir] ) {

			int	y = ir / w,
				x = ir - w * y;

			v[k] = vr[ir] / cnt;
			p[k] = Point( x, y );
			++k;
		}
	}

	p.resize( k );
	v.resize( k );
}

/* --------------------------------------------------------------- */
/* InterpolatePixel ---------------------------------------------- */
/* --------------------------------------------------------------- */

double InterpolatePixel(
	double					x,
	double					y,
	const uint8*			image,
	int						w )
{
	int		ix		= (int)x;
	int		iy		= (int)y;
	int		i		= ix + w * iy;
	double	alpha	= x - ix;
	double	beta	= y - iy;

	return	(1.0-alpha)*(1.0-beta)*image[i] +
		        (alpha)*(1.0-beta)*image[i + 1] +
			(1.0-alpha)*(    beta)*image[i + w] +
			    (alpha)*(    beta)*image[i + 1 + w];
}


double InterpolatePixel(
	double					x,
	double					y,
	const vector<double>	&image,
	int						w )
{
	int		ix		= (int)x;
	int		iy		= (int)y;
	int		i		= ix + w * iy;
	double	alpha	= x - ix;
	double	beta	= y - iy;

	return	(1.0-alpha)*(1.0-beta)*image[i] +
		        (alpha)*(1.0-beta)*image[i + 1] +
			(1.0-alpha)*(    beta)*image[i + w] +
			    (alpha)*(    beta)*image[i + 1 + w];
}

/* --------------------------------------------------------------- */
/* DistributePixel ----------------------------------------------- */
/* --------------------------------------------------------------- */

void DistributePixel(
	double			x,
	double			y,
	double			val,
	float*			image,
	int				w )
{
	int		ix		= (int)x;
	int		iy		= (int)y;
	int		i		= ix + w * iy;
	double	alpha	= x - ix;
	double	beta	= y - iy;

	image[i]			+= (1.0-alpha)*(1.0-beta)*val;
	image[i + 1]		+=     (alpha)*(1.0-beta)*val;
	image[i + w]		+= (1.0-alpha)*(    beta)*val;
	image[i + 1 + w]	+=     (alpha)*(    beta)*val;
}

/* --------------------------------------------------------------- */
/* BiCubicInterp ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Cubic interpolation:
// data[0] .. data[1] .. x .. data[2] .. data[3],
// where x is a fractional value.
//
template<class T>
static double _CubicInterp( const T *data, double x )
{
	double	a = data[0];
	double	b = data[1];
	double	c = data[2];
	double	d = data[3];
	double	x2 = x*x;
	double	x3 = x2*x;

	return	x3 * (-0.5*a + 1.5*b - 1.5*c + 0.5*d) +
			x2 * (a - 2.5*b + 2.0*c - 0.5*d) +
			 x * (-0.5*a + 0.5*c) +
			 b;
}


double BiCubicInterp( const double* image, int w, Point p )
{
	int		x	= int(p.x);
	int		y	= int(p.y);
	double	dx	= p.x - x;
	int		n	= (x-1) + w*(y-1);
	double	yi[4];

	yi[0] = _CubicInterp( image + n, dx );	// interpolate in X
	n += w;									// move up a line
	yi[1] = _CubicInterp( image + n, dx );	// interpolate in X
	n += w;									// etc.
	yi[2] = _CubicInterp( image + n, dx );
	n += w;
	yi[3] = _CubicInterp( image + n, dx );

	return _CubicInterp( yi, p.y - y );		// Now interpolate in y
}


double BiCubicInterp( const uint8* image, int w, Point p )
{
	int		x	= int(p.x);
	int		y	= int(p.y);
	double	dx	= p.x - x;
	int		n	= (x-1) + w*(y-1);
	double	yi[4];

	yi[0] = _CubicInterp( image + n, dx );	// interpolate in X
	n += w;									// move up a line
	yi[1] = _CubicInterp( image + n, dx );	// interpolate in X
	n += w;									// etc.
	yi[2] = _CubicInterp( image + n, dx );
	n += w;
	yi[3] = _CubicInterp( image + n, dx );

	return _CubicInterp( yi, p.y - y );		// Now interpolate in y
}


// Testing function
//
static void TestInterp()
{
	double	d[4] = {2,4,2,3};

	for( double x=-1; x<=2.01; x += 0.5 )
		printf( "%f -> %f\n", x, _CubicInterp( d, x ) );
}

/* --------------------------------------------------------------- */
/* ImageGradients ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Find the gradients across a set of points in a picture.
//
// Uses standard fit to line with free {slope, intercept}.
//
void ImageGradients(
	double				&slopex,
	double				&slopey,
	const vector<uint8>	&image,
	int					w,
	const vector<Point>	&pts,
	FILE				*flog )
{
	double	sum_p	= 0.0;
	double	sum_xp	= 0.0;
	double	sum_x	= 0.0;
	double	sum_x2	= 0.0;
	double	sum_yp	= 0.0;
	double	sum_y	= 0.0;
	double	sum_y2	= 0.0;
	int		N		= pts.size();

	for( int i = 0; i < N; ++i ) {

		int	ix		= int(pts[i].x);
		int	iy		= int(pts[i].y);
		int	pixel	= image[ix + w*iy];

		sum_p	+= pixel;
		sum_xp	+= ix*pixel;
		sum_x	+= ix;
		sum_x2	+= ix*ix;
		sum_yp	+= iy*pixel;
		sum_y	+= iy;
		sum_y2	+= iy*iy;
	}

	double	avgx = sum_x / N;
	slopex = (sum_xp - avgx * sum_p) / (sum_x2 - avgx * sum_x);

	double	avgy = sum_y / N;
	slopey = (sum_yp - avgy * sum_p) / (sum_y2 - avgy * sum_y);

	fprintf( flog,
	"Gradients: Slope bias per 1000 pixels: x %f, y %f.\n",
	slopex*1000.0, slopey*1000.0 );
}

/* --------------------------------------------------------------- */
/* IsLowContrast ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Test whether a 128x128 block is low contrast.
// Note that image (I) has already been pre-processed
// so has mean 0 and std 1.
//
bool IsLowContrast( vector<double> &I, double std )
{
	vector<double>	result;		// power as a function of scale
	vector<double>	smaller;
	vector<double>	tmp = I;
	int				i, x, y, whlf, w, lowcon = false;

	result.reserve( 8 );

	for( w = 128; w > 1; w = whlf ) {

		// sum squares at current size
		double mean, lstd;
		Stats( tmp, mean, lstd );
		result.push_back( lstd );

		// make a 1/2-size copy
		whlf = w / 2;
		smaller.resize( whlf * whlf );

		for( y = 0; y < whlf; ++y ) {

			for( x = 0; x < whlf; ++x ) {

				int		n = 2 * (x + w*y);

				smaller[x + whlf * y] =
				(tmp[n] + tmp[n+1] + tmp[n+w] + tmp[n+w+1]) / 4.0;
			}
		}

		tmp	= smaller;
	}

// print powers
	printf( "IsLowContrast: Amplitudes(%.2f):", std );

	for( i = 0; i < result.size(); ++i )
		printf( "%.3f ", result[i]*std );

	printf( "\n" );

// if the contrast is low in absolute terms ||
// after reducing 128x128 to 4x4 ||
// some other test??

	if( result[0] < 0.30 || result[5]*std < 2.0 )
		lowcon = true;

	return lowcon;
}


