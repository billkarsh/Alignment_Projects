

#include	"ImageIO.h"
#include	"Metrics.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"TAffine.h"

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static int FFTFileIdx = 0;






/* --------------------------------------------------------------- */
/* MeanSqrDiff --------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool MeanSqrDiff(
    const vector<double>	&av,
    const vector<double>	&bv,
    const char				*msg,
    FILE*					flog )
{
    double	sum	= 0.0;
    int		na	= av.size();

    if( !na ) {
        fprintf( flog, "Metrics: No data?\n" );
        return false;
    }

    for( int i = 0; i < na; ++i ) {

        double	d = av[i] - bv[i];

        sum += d*d;
    }

    fprintf( flog,
    "Metrics: %s: Mean square difference %f\n", msg, sum/na );

    return true;
}

/* --------------------------------------------------------------- */
/* SmallestFootprint --------------------------------------------- */
/* --------------------------------------------------------------- */

// Rotate set of points until they occupy the smallest area.
// (Allows more efficient image-based measurements).
//
// If zero is the best angle, then a pointer to the input vector
// is returned. If any other angle is better, then (newpts) gets
// the transformed points and a pointer to newpts is returned.
//
static const vector<Point>* SmallestFootprint(
    vector<Point>			&newpts,
    const vector<Point>		&pts,
    const char				*msg,
    FILE*					flog )
{
    DBox	B;
    int		deg = TightestBBox( B, pts );

    fprintf( flog,
    "Metrics: %s: Smallest footprint deg=%d, area=%g\n",
    msg, deg, (B.R - B.L) * (B.T - B.B) );

    if( deg == 0 )
        return &pts;
    else {

        TAffine	T;

        T.NUSetRot( deg*PI/180 );
        T.Apply_R_Part( newpts = pts );

        return &newpts;
    }
}

/* --------------------------------------------------------------- */
/* MakeMetricImagesFFT ------------------------------------------- */
/* --------------------------------------------------------------- */

// For the input (pts) region and value lists {av, bv):
// - Return images {i1, i2, diff} (and dims Nx, Ny).
// - Optionally write image files.
//
static void MakeMetricImagesFFT(
    vector<double>			&i1,
    vector<double>			&i2,
    vector<double>			&diff,
    int						&Nx,
    int						&Ny,
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    bool					write_images,
    const char				*msg,
    FILE*					flog )
{
/* ---------- */
/* Initialize */
/* ---------- */

    i1.clear();
    i2.clear();
    diff.clear();
    Nx	= 0;
    Ny	= 0;

/* -------------------- */
/* Set image dimensions */
/* -------------------- */

    IBox	B;

    BBoxFromPoints( B, pts );

    Nx = CeilPow2( B.R - B.L + 1 );
    Ny = CeilPow2( B.T - B.B + 1 );

    int	N2 = Nx * Ny;

    fprintf( flog,
    "MetricImages: Range x %d %d, y %d %d, use Nx=%d, Ny=%d\n",
    B.L, B.R, B.B, B.T, Nx, Ny );

/* ----------- */
/* Fill images */
/* ----------- */

    i1.resize( N2, 0.0 );
    i2.resize( N2, 0.0 );
    diff.resize( N2 );

    int	np = pts.size();

    for( int i = 0; i < np; ++i ) {

        double	x = pts[i].x - B.L,
                y = pts[i].y - B.B;

        DistributePixel( x, y, av[i], i1, Nx, Ny );
        DistributePixel( x, y, bv[i], i2, Nx, Ny );
    }

    for( int i = 0; i < N2; ++i )
        diff[i] = i1[i] - i2[i];

/* ------------ */
/* Write images */
/* ------------ */

    if( write_images ) {

        char	fname[32];

        sprintf( fname, "fft%d-a.tif", FFTFileIdx );
        VectorDblToTif8( fname, i1, Nx, Ny, flog );

        sprintf( fname, "fft%d-b.tif", FFTFileIdx );
        VectorDblToTif8( fname, i2, Nx, Ny, flog );

        sprintf( fname, "fft%d-d.tif", FFTFileIdx );
        VectorDblToTif8( fname, diff, Nx, Ny, flog );

        ++FFTFileIdx;

        double	e1 = 0.0, e2 = 0.0, ed = 0.0;

        for( int i = 0; i < N2; ++i ) {

            e1 += i1[i]*i1[i];
            e2 += i2[i]*i2[i];
            ed += diff[i]*diff[i];
        }

        fprintf( flog,
        "MetricImages: Energies (1,2,dif) %f %f %f\n",
        e1, e2, ed );
    }
}

/* --------------------------------------------------------------- */
/* MakeMetricImagesEMM ------------------------------------------- */
/* --------------------------------------------------------------- */

// For the input (pts) region and value lists {av, bv):
// - Return image diff (and dims Nx, Ny).
// - Optionally write image file.
//
static void MakeMetricImagesEMM(
    vector<double>			&diff,
    int						&Nx,
    int						&Ny,
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    bool					write_images,
    const char				*msg,
    FILE*					flog )
{
/* ---------- */
/* Initialize */
/* ---------- */

    diff.clear();
    Nx	= 0;
    Ny	= 0;

/* -------------------- */
/* Set image dimensions */
/* -------------------- */

    IBox	B;

    BBoxFromPoints( B, pts );

    Nx = B.R - B.L + 1;
    Ny = B.T - B.B + 1;

// make even

    Nx += (Nx & 1);
    Ny += (Ny & 1);

    int	N2 = Nx * Ny;

    fprintf( flog,
    "MetricImages: Range x %d %d, y %d %d, use Nx=%d, Ny=%d\n",
    B.L, B.R, B.B, B.T, Nx, Ny );

/* ----------- */
/* Fill images */
/* ----------- */

    diff.resize( N2, 0.0 );

    int	np = pts.size();

    for( int i = 0; i < np; ++i ) {

        double	x = pts[i].x - B.L,
                y = pts[i].y - B.B;

        DistributePixel( x, y, av[i] - bv[i], diff, Nx, Ny );
    }

/* ------------ */
/* Write images */
/* ------------ */

    if( write_images ) {

        char	fname[32];

        sprintf( fname, "fft%d-d.tif", FFTFileIdx );
        VectorDblToTif8( fname, diff, Nx, Ny, flog );

        ++FFTFileIdx;

        double	ed = 0.0;

        for( int i = 0; i < N2; ++i )
            ed += diff[i]*diff[i];

        fprintf( flog, "MetricImages: Energies (dif) %f\n", ed );
    }
}

/* --------------------------------------------------------------- */
/* UnnormHaarTf -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Un-normalized Haar transform
//
static void UnnormHaarTf( double *data, int n )
{
    vector<double>	tmp( n );

    while( n > 1 ) {

        n /= 2;

        for( int i = 0; i < n; ++i ) {

            int	k = 2 * i;

            tmp[i]		= data[k] + data[k+1];
            tmp[i+n]	= data[k] - data[k+1];
        }

        memcpy( &data[0], &tmp[0], n * 2 * sizeof(double) );
    }
}

/* --------------------------------------------------------------- */
/* UnnormHaarTf2D ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void UnnormHaarTf2D( vector<double> &data, int w, int h )
{
// transform every row

    for( int y = 0; y < h; ++y )
        UnnormHaarTf( &data[w*y], w );

//	PrintVectorAsMat( stdout, data, w );

// now every column

    vector<double>	col( h );

    for( int x = 0; x < w; ++x ) {

        for( int y = 0; y < h; ++y )
            col[y] = data[x + w*y];

        UnnormHaarTf( &col[0], h );

        for( int y = 0; y < h; ++y )
            data[x + w*y] = col[y];
    }
}

/* --------------------------------------------------------------- */
/* ApproxEarthMoversMetric --------------------------------------- */
/* --------------------------------------------------------------- */

static double ApproxEarthMoversMetric(
    vector<double>	&dd,
    int				wi,
    int				hi,
    FILE*			flog )
{
// find powers of 2 that are big enough

    int	w = CeilPow2( wi ),
        h = CeilPow2( hi );

    vector<double>	d( w*h, 0.0 );

    CopyRaster( &d[0], w, &dd[0], wi, wi, hi );

    UnnormHaarTf2D( d, w, h );

//	PrintVectorAsMat( stdout, d, 8 );

    vector<int>	cdist( w, 1 );
    vector<int>	rdist( h, 1 );
    int			mask;

    mask = w >> 1;  // only works for w = 2^n

    for( int x = 1; x < w; ++x ) {

        for( int tmp = x; !(tmp & mask); tmp <<= 1 )
            cdist[x] <<= 1;
    }

    mask = h >> 1;  // only works for h = 2^n

    for( int y = 1; y < h; ++y ) {

        for( int tmp = y; !(tmp & mask); tmp <<= 1 )
            rdist[y] <<= 1;
    }

// use of 1 as lower bound is right;
// do not care about sums, only differences

    double	approx = 0.0;  // approximate Earth Mover's metric

    for( int x = 1; x < w; ++x ) {

        for( int y = 1; y < h; ++y ) {

            int		dx		= cdist[x];
            int		dy		= rdist[y];
            double	dist	= min( dx, dy );

            approx += dist * abs( d[x + w*y] );
        }
    }

// Normalize:
// - divide by 2 deltas per move unit.
// - divide by N pixels (per pixel cost).
// - divide by Poisson noise on area N.

    int		N = wi * hi;

    approx /= 2.0 * N * sqrt( N );

    fprintf( flog,
    "Approximate EM metric %f for %d points.\n",
    approx, N );

    return approx;
}

/* --------------------------------------------------------------- */
/* FFT_r2c_lookup ------------------------------------------------ */
/* --------------------------------------------------------------- */

static CD FFT_r2c_lookup(
    vector<CD>	&c,
    int			Nx,
    int			Ny,
    int			x,
    int			y )
{
    int M = Nx/2+1;

    if( x < 0 ) x += Nx;
    if( y < 0 ) y += Ny;

    if( x > Nx/2 )
        return conj( c[(Nx-x) + M*y] );

    return c[x + M*y];
}

/* --------------------------------------------------------------- */
/* EarthMoversMetric --------------------------------------------- */
/* --------------------------------------------------------------- */

double EarthMoversMetric(
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    bool					write_images,
    const char				*msg,
    FILE*					flog )
{
/* ----------------------- */
/* Report difference power */
/* ----------------------- */

    if( !MeanSqrDiff( av, bv, msg, flog ) )
        return 0.0;

/* ----------- */
/* Make images */
/* ----------- */

    vector<double>	diff;
    int				Nx, Ny;

    {
        const vector<Point>	*pbest;
        vector<Point>		altpts;

        pbest = SmallestFootprint( altpts, pts, msg, flog );

        MakeMetricImagesEMM( diff, Nx, Ny,
            *pbest, av, bv, write_images, msg, flog );
    }

/* ------- */
/* Measure */
/* ------- */

    return ApproxEarthMoversMetric( diff, Nx, Ny, flog );
}

/* --------------------------------------------------------------- */
/* FourierMatch -------------------------------------------------- */
/* --------------------------------------------------------------- */

//static void AddJunkToDiff( vector<double> diff, int Nx, int Ny )
//{
//	double	av, sd;
//	Stats( diff, av, sd );
//
//	printf( "***** junk av, std = %f, %f\n", av, sd );
//
//	const int size = 30;
//	const int nobj = 1000;
//
//	for( int i = 0; i < nobj; ++i ) {
//
//		int	x0 = int( (Nx-size-1)*rand()/RAND_MAX );
//		int	y0 = int( (Ny-size-1)*rand()/RAND_MAX );
//
//		for( int y = y0; y < y0+size; ++y ) {
//
//			for( int x = x0; x < x0+size; ++x ) {
//
//				diff[x+Nx*y] = av + 3*sd;
//			}
//		}
//	}
//}


// See if two images match in the Fourier domain.
//
// Vectors av[] and bv[] should be normalized and same size.
//
// wvlen (in pixels) sets a minimum feature size (roughly).
//
double FourierMatch(
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    int						wvlen,
    bool					write_images,
    const char				*msg,
    FILE*					flog )
{
/* ----------------------- */
/* Report difference power */
/* ----------------------- */

    if( !MeanSqrDiff( av, bv, msg, flog ) )
        return 0.0;

/* ----------- */
/* Make images */
/* ----------- */

    vector<double>	i1, i2, diff;
    int				Nx, Ny;

    {
        const vector<Point>	*pbest;
        vector<Point>		altpts;

        pbest = SmallestFootprint( altpts, pts, msg, flog );

        MakeMetricImagesFFT( i1, i2, diff, Nx, Ny,
            *pbest, av, bv, write_images, msg, flog );
    }

/* ----------- */
/* Image power */
/* ----------- */

// Compare lowest few coefficients (longer than wvlen).
//
// Remember FFT x-elements are arranged like this:
//
//	elem	freq		wvlen
//	----	----		-----
//	0		0			DC const
//	1		1/Nx		Nx
//	2		2/Nx		Nx/2
//	3		3/Nx		Nx/3
// ...		...			...
//	Nx/2	(Nx/2)/Nx	2
//
// Since wvlen = Nx/elem; then elem = Nx/wvlen.

    vector<CD>	i1fft, i2fft, dfft;
    double		total = 0.0, dot = 0.0;
    int			xlim = Nx/wvlen, ylim = Ny/wvlen;

    FFT_2D( i1fft,  i1, Nx, Ny, false );
    FFT_2D( i2fft,  i2, Nx, Ny, false );
    FFT_2D( dfft, diff, Nx, Ny, false );

    for( int x = -xlim; x <= xlim; ++x ) {

        for( int y = -ylim; y <= ylim; ++y ) {

            CD	v1 = FFT_r2c_lookup( i1fft, Nx, Ny, x, y ),
                v2 = FFT_r2c_lookup( i2fft, Nx, Ny, x, y );

            total	+= sqrt( norm( v1 ) * norm( v2 ) );
            dot		+= v1.real()*v2.real() + v1.imag()*v2.imag();
        }
    }

    dot /= total;

    fprintf( flog, "FFT: norm-dot %f, energy %f\n", dot, total );

/* --------------------------------*/
/* Difference power: make spectrum */
/* --------------------------------*/

// Form spectrum with powers ordered from small to large wavelength.
// The objective will be to report what fraction of total power is
// covered by small-sized features in the diff image. The idea is
// that in a seriously bad mismatch, diff will get large features
// and the spectrum would be shifted out to long wavelength.

    int				M		= Nx/2 + 1;
    int				maxf	= int(sqrt( Nx*Nx + Ny*Ny )/2.0) + 1;
    vector<double>	pspectrum( maxf, 0.0 );

    for( int x = 1; x < M; ++x ) {

        double	wavex = double(Nx)/x;

        for( int y = 1; y < Ny; ++y ) {

            double	wavey = double(Ny)/(y > Ny/2 ? Ny-y : y);
            double	wave  = sqrt( wavex*wavex + wavey*wavey );

            int	iwave = int(wave);

            if( iwave > maxf - 1 )
                iwave = maxf - 1;

            pspectrum[iwave] += norm( dfft[x + M*y] );
        }
    }

/* --------------------------------------------------- */
/* Difference power: where cum power exceeds threshold */
/* --------------------------------------------------- */

    const double thresh = 0.50;	// fraction of total power

    double	tot = 0.0, cum = 0.0;
    int		i;

// tot = total power

    for( i = 0; i < maxf; ++i )
        tot += pspectrum[i];

// repeat summing, but only up to tot*thresh

    for( i = 0; cum < tot * thresh && i < maxf; ++i ) {

        cum += pspectrum[i];

        // report progress periodically
        if( i && !(i % 10) ) {
            fprintf( flog,
            "FFT: %s: cum frac to %d = %f\n", msg, i, cum/tot );
        }
    }

    fprintf( flog,
    "FFT: %s: Cum power exceeds %2d%% of %f"
    " at index %d/%d (frac %f).\n",
    msg, int(thresh*100.0), tot, i, maxf, cum/tot );

    return dot;
}

/* --------------------------------------------------------------- */
/* PercentYellow ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Try another metric. By eye, we say a match between a red and
// a green image is good if it has a lot of yellow points. Here we
// report a percentage of yellow pixels (that is, more yellow than
// red or green to some tolerance). Several tolerances are tried
// for now since best tolerance value is not well characterized.
// This calculation requires normalized data...
//
// IMPORTANT: Caller must call Normalize() on the data before
// passing them here.
//
double PercentYellow(
    const vector<double>	&a,
    const vector<double>	&b,
    FILE*					flog )
{
    int yellow, red, green, N = a.size();

    for( double tol = 1.05; tol < 1.255; tol += 0.05 ) {

        yellow	= 0;
        red		= 0;
        green	= 0;

        for( int i = 0; i < N; ++i ) {

            int	pa = 127 + int(a[i]*35.0);
            int	pb = 127 + int(b[i]*35.0);

            if( pa > tol * pb )
                ++red;
            else if( pa * tol < pb )
                ++green;
            else
                ++yellow;
        }

        fprintf( flog,
        "%%Yellow: Tol %6.2f  red %6.1f  yellow %6.1f  green %6.1f\n",
        tol, 100.0*red/N, 100.0*yellow/N, 100.0*green/N );
    }

    return (double)yellow / N;
}


