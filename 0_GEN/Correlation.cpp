

#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"ImageIO.h"
#include	"Debug.h"

#include	<math.h>
#include	<pthread.h>
#include	<stdlib.h>
#include	<string.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// mutex_fft guards {fftw interface, cached fft2 array}
static pthread_mutex_t	mutex_fft = PTHREAD_MUTEX_INITIALIZER;
static int _dbg_simgidx = 0;






/* --------------------------------------------------------------- */
/* FFT ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

#ifdef ALN_USE_MKL
#include	"Correlation_fft_mkl.cpp"
#else
#include	"Correlation_fft_fftw.cpp"
#endif

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

    S.assign(  N, 0.0 );
    S2.assign( N, 0.0 );
    nz.assign( N, 0 );

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
    int		x, y, i, j, k, m;

    nz.assign( wS * hS, 0 );

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
    "\nDebugLinearCorr: s %f %f, a %f %f, ss %f %f %f, r %f\n",
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
//	- n1:	size of source 'a' data
//	- n2:	size of target 'b' data
//
// Returns:
//	Storage allocation size that should be used for each array:
//	{'a' data, 'b' data, correlation results}. Return value
//	accounts for needed zero padding and is a power of two.
//
// Result data ordering
// --------------------
// Array sizes are N = n1 + n2 - 1.
//
//	R[0]		= lag 0
//	R[1]		= lag +1
//	...
//	R[n2-1]		= lag +(n2-1)	=> (n2) lags are >= zero.
//
//
//	R[N-1]		= lag -1
//	R[N-2]		= lag -2
//	...
//	R[N-(n1-1)]	= lag -(n1-1)	=> (n1-1) lags are < zero.
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
int FFTSize( int n1, int n2 )
{
    return CeilPow2( n1 + n2 - 1 );
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
/* Convolve ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Convolve src with filter (kernel) (response function) K.
// Optionally normalize values in K. If kfft is non-empty and
// correct size it is assumed valid.
//
// dst can be same as src if desired.
//
void Convolve(
    vector<double>			&dst,
    const vector<double>	&src,
    int						ws,
    int						hs,
    const double			*K,
    int						wk,
    int						hk,
    bool					kIsSymmetric,
    bool					preNormK,
    vector<CD>				&kfft,
    FILE					*flog )
{
    int	Ns	= ws * hs,
        Nk	= wk * hk,
        Nx	= FFTSize( ws, wk / (1 + kIsSymmetric) ),
        Ny	= FFTSize( hs, hk / (1 + kIsSymmetric) ),
        Nxy	= Nx * Ny,
        M	= Ny*(Nx/2+1);

// Prepare K fft

    if( kfft.size() != M ) {

        // workspace image
        vector<double>	KK( Nxy, 0.0 );
        double			nrm = 1.0;

        // normalization
        if( preNormK ) {

            nrm = 0.0;

            for( int i = 0; i < Nk; ++i )
                nrm += K[i];

            if( !nrm )
                nrm = 1.0;
        }

        // load workspace
        int	w2 = wk / 2,
            h2 = hk / 2;

        for( int i = -h2; i <= h2; ++i ) {

            int	y = (i >= 0 ? i : Ny + i);

            for( int j = -w2; j <= w2; ++j ) {

                int	x = (j >= 0 ? j : Nx + j);

                KK[x + Nx*y] = K[j+w2 + wk*(i+h2)] / nrm;
            }
        }

        // FFT
        FFT_2D( kfft, KK, Nx, Ny, false, flog );
    }

// Prepare src fft

    vector<double>	SS( Nxy, 0.0 );
    vector<CD>		sfft;

    CopyRaster( &SS[0], Nx, &src[0], ws, ws, hs );

    FFT_2D( sfft, SS, Nx, Ny, false, flog );

// Convolve

    for( int i = 0; i < M; ++i )
        sfft[i] *= kfft[i];

    IFT_2D( SS, sfft, Nx, Ny, flog );

// Copy back

    dst.resize( Ns );

    CopyRaster( &dst[0], ws, &SS[0], Nx, ws, hs );

// Normalize

    for( int i = 0; i < Ns; ++i )
        dst[i] /= Nxy;
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

    if( ix - d >= 0 &&
        ix + d < nx &&
        (z0 = LkFFT( I, nx, ny, ix - d, iy )) &&
        (z2 = LkFFT( I, nx, ny, ix + d, iy )) ) {

        z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

        if( fabs( z0 ) <= 8 * d )
            xpk += z0;
    }

    if( iy - d >= 0 &&
        iy + d < ny &&
        (z0 = LkFFT( I, nx, ny, ix, iy - d )) &&
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
        "Landscape: Max %f at (%d, %d); Ox Oy radius = %d %d %d\n",
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
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
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
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
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
        "NormCorr: Very odd - ir=%d rslt[i]=%f olap_area=%ld\n",
        ir, r, (long)n );

        fprintf( flog,
        "NormCorr: shift %d %d, i1 (%d %d) to (%d %d),"
        " i2 (%d %d) to (%d %d).\n",
        dx, dy,
        OL1.L, OL1.B, OL1.R, OL1.T,
        OL2.L, OL2.B, OL2.R, OL2.T );

        fprintf( flog,
        "NormCorr: sums:      %f %f %f %f\n"
        "NormCorr: num d1 d2: %f %f %f\n",
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
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
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
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
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

// Return cross-correlation and additive displacement
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

    int	Nx	= FFTSize( w1, w2 ),
        Ny	= FFTSize( h1, h2 ),
        M	= Ny*(Nx/2+1);

    if( verbose )
        fprintf( flog, "NormCorr: Nx = %d, Ny = %d\n", Nx, Ny );

// Create images from point lists.

    vector<double>	i1, i2;

    ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
    ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

    vector<double>	rslt;
    vector<CD>		fft1;

#ifdef ALN_USE_MKL

    vector<CD>	_fft2;

    FFT_2D( _fft2, i2, Nx, Ny, false, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    for( int i = 0; i < M; ++i )
        fft1[i] = _fft2[i] * conj( fft1[i] );

#else

    FFT_2D( fft2, i2, Nx, Ny, true, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    pthread_mutex_lock( &mutex_fft );

    for( int i = 0; i < M; ++i )
        fft1[i] = fft2[i] * conj( fft1[i] );

    pthread_mutex_unlock( &mutex_fft );

#endif

    IFT_2D( rslt, fft1, Nx, Ny, flog );

// Create array indexed by 'size': int( 10 * log( overlap_size ) ).
// Each element contains the best rslt[i] at that size index.
// The idea is to look for correlation peaks that are nearly
// as good as the best, but that derive from larger overlaps.

    int lmax = (int)(10.0*log(double(Nx*Ny))) + 1;
    vector<double>	max_by_size( lmax, 0.0 );

// Prepare correlation calculator

    CLinCorr	ccalc;	// uses Pearson's r
//	CCrossCorr	ccalc;	// uses 1/n * SUM(a*b)

    ccalc.Initialize( flog, i1, w1, h1, i2, w2, h2, Nx, Ny );

// Now for the conventional biggest

    double	biggest	= -1.0E30;
    int		bigx	= -1,
            bigy	= -1,
            nnegx	= w1 - 1,
            nnegy	= h1 - 1,
            nposx	= w2,
            nposy	= h2,
            rsqr	= radius * radius;

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

        if( verbose ) {

            fprintf( flog,
            "NormCorr: No legal subregions at all...\n" );

            fprintf( flog,
            "NormCorr: Maximum correlation of %f at [%d,%d].\n",
            0.0, 0, 0 );
        }

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
            " %8d - %8d : %8.2f\n", i1, i2, max_by_size[i] );
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
                "NormCorr: Local max at %d %d, value %f\n",
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

// Return the standard cross-correlation and (dx, dy), which
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

    M =	FFT_2D( fft2, i2, N, N, false );
        FFT_2D( fft1, i1, N, N, false );

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

        printf( "FindCor: Very odd - norm-rslt x y %f %d %d\n",
            biggest, bigx, bigy );
        exit( 44 );
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
    dRdc.assign( nc, Point( 0.0, 0.0 ) );

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
        "CorrVectors: Sizes differ! %d %ld\n", N, b.size() );
        exit( 42 );
    }

    for( int i = 0; i < N; ++i ) {

        double	prod = a[i]*b[i];

        sum2 += prod;

        if( fabs( prod ) < 1.0E-8 )
            ++nz;
    }

    if( !isfinite( sum2 ) ) {
        fprintf( flog,
        "CorrVectors: Likely all zero pixels, corr sum = %g"
        " ... returning corr = zero.\n", sum2 );

        nz		= N;
        sum2	= 0.0;
    }

//	fprintf( flog, "CorrVectors: %d of %d were small.\n", nz, N );

    N -= nz;

    if( nnz )
        *nnz = N;

    return (N ? sum2 / N : 0.0);
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
    "STAT: ImproveCpt: Initial %s correlation %f (%ld pixels).\n",
    describe, corr, av.size() );

// Plausibility check

    if( corr < iniThresh ) {

        fprintf( flog,
        "FAIL: ImproveCpt: Correlation %f less than %f at start.\n",
        corr, iniThresh );

        PrintControlPoints( flog, ac );
        //PrintPixels( flog, am, av, bnew );

        return 0.0;
    }

// Skip optimizing if finThresh < 0

    if( finThresh < 0 ) {

        fprintf( flog,
        "STAT: ImproveCpt: Skipping optimizer; final corr %f\n",
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
/* RCalc --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class RCalc {

private:
    vector<double>	i1sum, i1sum2,
                    i2sum, i2sum2;
    vector<int>		i1nz,  i2nz;
    int				w1,  h1,
                    w2,  h2,
                    Nx,  Ny,
                    Nxy;
    IBox			OL1, OL2;
    int				olw, olh;

public:
    void Initialize(
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
        int						w2,
        int						h2,
        int						Nx,
        int						Ny );

    bool Valid(
        EvalType	LegalRgn,
        void*		arglr,
        EvalType	LegalCnt,
        void*		arglc,
        int			dx,
        int			dy );

    double CalcR( double rslt );

    inline double CalcS( double rslt )
        {return rslt / Nxy;};
};


void RCalc::Initialize(
        const vector<double>	&I1,
        int						w1,
        int						h1,
        const vector<double>	&I2,
        int						w2,
        int						h2,
        int						Nx,
        int						Ny )
{
    this->w1	= w1;
    this->h1	= h1;
    this->w2	= w2;
    this->h2	= h2;
    this->Nx	= Nx;
    this->Ny	= Ny;
    Nxy			= Nx * Ny;

    IntegrateImage( i1sum, i1sum2, i1nz, w1, h1, I1, Nx );
    IntegrateImage( i2sum, i2sum2, i2nz, w2, h2, I2, Nx );
}


bool RCalc::Valid(
        EvalType	LegalRgn,
        void*		arglr,
        EvalType	LegalCnt,
        void*		arglc,
        int			dx,
        int			dy )
{
    int		ok = true;

    BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

    olw = OL1.R - OL1.L + 1;
    olh = OL1.T - OL1.B + 1;

    if( LegalRgn && !LegalRgn( olw, olh, arglr ) )
        ok = false;

// Large enough density of non-zero values?

    if( LegalCnt ) {

        int	i1c = IntegralTable( i1nz, w1, OL1 );
        int	i2c = IntegralTable( i2nz, w2, OL2 );

        if( !LegalCnt( i1c, i2c, arglc ) )
            ok = false;
    }

    return ok;
}


double RCalc::CalcR( double rslt )
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

    return (r > -1.0 && r < 1.0 ? r : 0.0);
}

/* --------------------------------------------------------------- */
/* CCorImg ------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CCorImg {

private:
    FILE	*flog;
    int		verbose;
    IBox	B1, B2;
    int		w1, h1, w2, h2,
            nnegx, nposx,
            nnegy, nposy,
            cx, cy, wR, hR, nR;

public:
    bool SetDims(
        FILE					*iflog,
        int						iverbose,
        const vector<Point>		&ip1,
        const vector<Point>		&ip2 );

    void MaskA(
        vector<uint8>			&A,
        int						Ox,
        int						Oy,
        int						Rx,
        int						Ry );

    void MakeRandA(
        vector<double>			&R,
        vector<uint8>			&A,
        const vector<Point>		&ip1,
        const vector<double>	&iv1,
        const vector<Point>		&ip2,
        const vector<double>	&iv2,
        EvalType				LegalRgn,
        void*					arglr,
        EvalType				LegalCnt,
        void*					arglc,
        int						Ox,
        int						Oy,
        int						Rx,
        int						Ry,
        vector<CD>				&fft2 );

    void MakeSandRandA(
        vector<double>			&S,
        vector<double>			&R,
        vector<uint8>			&A,
        const vector<Point>		&ip1,
        const vector<double>	&iv1,
        const vector<Point>		&ip2,
        const vector<double>	&iv2,
        EvalType				LegalRgn,
        void*					arglr,
        EvalType				LegalCnt,
        void*					arglc,
        int						Ox,
        int						Oy,
        int						Rx,
        int						Ry,
        vector<CD>				&fft2 );

    void MakeF(
        vector<double>			&F,
        vector<uint8>			&A,
        const vector<double>	&R );

    bool OrderF(
        vector<int>				&forder,
        const vector<double>	&F,
        const vector<uint8>		&A,
        const vector<double>	&R,
        double					mincor );

    bool FPeak(
        int						&rx,
        int						&ry,
        const vector<int>		&forder,
        const vector<double>	&F,
        const vector<double>	&R,
        double					nbmaxht );

    void SimpleMax(
        int						&rx,
        int						&ry,
        const vector<int>		&forder );

    double ReturnR(
        double					&dx,
        double					&dy,
        int						rx,
        int						ry,
        const vector<double>	&Rreport,
        const vector<double>	&Rtweak );
};

/* --------------------------------------------------------------- */
/* CCorImg::SetDims ---------------------------------------------- */
/* --------------------------------------------------------------- */

bool CCorImg::SetDims(
    FILE					*iflog,
    int						iverbose,
    const vector<Point>		&ip1,
    const vector<Point>		&ip2 )
{
    flog	= iflog;
    verbose	= iverbose;

    BBoxFromPoints( B1, ip1 );
    BBoxFromPoints( B2, ip2 );

    w1 = B1.R - B1.L + 1;
    h1 = B1.T - B1.B + 1;

    w2 = B2.R - B2.L + 1;
    h2 = B2.T - B2.B + 1;

    nnegx = w1 - 1;
    nposx = w2;
    nnegy = h1 - 1;
    nposy = h2;

    cx = nnegx,
    cy = nnegy,
    wR = nnegx + nposx,
    hR = nnegy + nposy,
    nR = wR * hR;

    if( verbose ) {

        fprintf( flog,
        "Corr: Region size is [%d %d] in x, [%d %d] in y.\n",
        B1.L, B1.R, B1.B, B1.T );

        fprintf( flog,
        "Corr: Target size is [%d %d] in x, [%d %d] in y.\n",
        B2.L, B2.R, B2.B, B2.T );
    }

// Need enough room to apply typical kernel

    const int mindelta = 11/2 + 1;

    if( nnegx < mindelta || nposx < mindelta ||
        nnegy < mindelta || nposy < mindelta ) {

        if( verbose )
            fprintf( flog, "Corr: Too small to evaluate.\n" );

        return false;
    }

    return true;
}

/* --------------------------------------------------------------- */
/* CCorImg::MaskA ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Oval described as x^2 + (a^2)*(y^2) = Rx^2, where a = Rx/Ry.
//
void CCorImg::MaskA(
    vector<uint8>			&A,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry )
{
    if( Rx <= 0 || Ry <= 0 )
        return;

    Ox += B1.L - B2.L + cx;
    Oy += B1.B - B2.B + cy;

    int	xmin = Ox - Rx,
        xmax = Ox + Rx,
        ymin = Oy - Ry,
        ymax = Oy + Ry;

    double	rx2 = Rx * Rx,
            asq = rx2 / ((double)Ry * Ry);

    for( int i = 0; i < nR; ++i ) {

        int	y = i / wR;

        if( y < ymin || y > ymax ) {
            A[i] = 0;
            continue;
        }

        int	x = i - wR * y;

        if( x < xmin || x > xmax ) {
            A[i] = 0;
            continue;
        }

        x -= Ox;
        y -= Oy;

        if( x*x + asq*y*y > rx2 )
            A[i] = 0;
    }
}

/* --------------------------------------------------------------- */
/* CCorImg::MakeRandA -------------------------------------------- */
/* --------------------------------------------------------------- */

void CCorImg::MakeRandA(
    vector<double>			&R,
    vector<uint8>			&A,
    const vector<Point>		&ip1,
    const vector<double>	&iv1,
    const vector<Point>		&ip2,
    const vector<double>	&iv2,
    EvalType				LegalRgn,
    void*					arglr,
    EvalType				LegalCnt,
    void*					arglc,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry,
    vector<CD>				&fft2 )
{
// Get array sizes (Nx,Ny) and FFT size M

    int	Nx, Ny, M;

    Nx = FFTSize( w1, w2 ),
    Ny = FFTSize( h1, h2 ),
    M  = Ny*(Nx/2+1);

    if( verbose )
        fprintf( flog, "Corr: Nx = %d, Ny = %d\n", Nx, Ny );

// Create images from point lists

    vector<double>	i1, i2;

    ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
    ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

    vector<double>	rslt;
    vector<CD>		fft1;

#ifdef ALN_USE_MKL

    vector<CD>	_fft2;

    FFT_2D( _fft2, i2, Nx, Ny, false, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    for( int i = 0; i < M; ++i )
        fft1[i] = _fft2[i] * conj( fft1[i] );

#else

    FFT_2D( fft2, i2, Nx, Ny, true, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    pthread_mutex_lock( &mutex_fft );

    for( int i = 0; i < M; ++i )
        fft1[i] = fft2[i] * conj( fft1[i] );

    pthread_mutex_unlock( &mutex_fft );

#endif

    IFT_2D( rslt, fft1, Nx, Ny, flog );

// Prepare correlation calculator

    RCalc	calc;

    calc.Initialize( i1, w1, h1, i2, w2, h2, Nx, Ny );

// Reorganize valid entries of rslt image so that (dx,dy)=(0,0)
// is at the image center and (cx,cy)+(dx,dy) indexes all pixels
// of any sign.

    R.resize( nR );
    A.resize( nR );

    double	vmin = 1e7, vmax = -1e7;

    for( int y = -nnegy; y < nposy; ++y ) {

        int	iy = Nx * (y >= 0 ? y : Ny + y);

        for( int x = -nnegx; x < nposx; ++x ) {

            int	ix = (x >= 0 ? x : Nx + x);
            int	ir = cx+x + wR*(cy+y);

            A[ir] = calc.Valid( LegalRgn, arglr,
                        LegalCnt, arglc, x, y );

            R[ir] = calc.CalcR( rslt[ix+iy] );

            if( R[ir] < vmin )
                vmin = R[ir];

            if( R[ir] > vmax )
                vmax = R[ir];
        }
    }

    if( dbgCor )
        fprintf( flog, "Corr: Center = (%d %d).\n", cx, cy );

    if( verbose ) {
        fprintf( flog,
        "Corr: R image range [%11.6f %11.6f].\n", vmin, vmax );
    }

    MaskA( A, Ox, Oy, Rx, Ry );

    if( dbgCor ) {

        char	simg[32];

        sprintf( simg, "thmA_%d.tif", _dbg_simgidx );
        CorrThmToTif8( simg, i1, Nx, w1, h1, flog );

        sprintf( simg, "thmB_%d.tif", _dbg_simgidx );
        CorrThmToTif8( simg, i2, Nx, w2, h2, flog );

        sprintf( simg, "corr_A_%d.tif", _dbg_simgidx );
        Raster8ToTif8( simg, &A[0], wR, hR, flog );

        sprintf( simg, "corr_R_%d.tif", _dbg_simgidx );
        RasterDblToTifFlt( simg, &R[0], wR, hR, flog );
    }
}

/* --------------------------------------------------------------- */
/* CCorImg::MakeSandRandA ---------------------------------------- */
/* --------------------------------------------------------------- */

void CCorImg::MakeSandRandA(
    vector<double>			&S,
    vector<double>			&R,
    vector<uint8>			&A,
    const vector<Point>		&ip1,
    const vector<double>	&iv1,
    const vector<Point>		&ip2,
    const vector<double>	&iv2,
    EvalType				LegalRgn,
    void*					arglr,
    EvalType				LegalCnt,
    void*					arglc,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry,
    vector<CD>				&fft2 )
{
// Get array sizes (Nx,Ny) and FFT size M

    int	Nx, Ny, M;

    Nx = FFTSize( w1, w2 ),
    Ny = FFTSize( h1, h2 ),
    M  = Ny*(Nx/2+1);

    if( verbose )
        fprintf( flog, "Corr: Nx = %d, Ny = %d\n", Nx, Ny );

// Create images from point lists

    vector<double>	i1, i2;

    ImageFromValuesAndPoints( i1, Nx, Ny, iv1, ip1, B1.L, B1.B );
    ImageFromValuesAndPoints( i2, Nx, Ny, iv2, ip2, B2.L, B2.B );

// FFTs and lags

    vector<double>	rslt;
    vector<CD>		fft1;

#ifdef ALN_USE_MKL

    vector<CD>	_fft2;

    FFT_2D( _fft2, i2, Nx, Ny, false, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    for( int i = 0; i < M; ++i )
        fft1[i] = _fft2[i] * conj( fft1[i] );

#else

    FFT_2D( fft2, i2, Nx, Ny, true, flog );
    FFT_2D( fft1, i1, Nx, Ny, false, flog );

    pthread_mutex_lock( &mutex_fft );

    for( int i = 0; i < M; ++i )
        fft1[i] = fft2[i] * conj( fft1[i] );

    pthread_mutex_unlock( &mutex_fft );

#endif

    IFT_2D( rslt, fft1, Nx, Ny, flog );

// Prepare correlation calculator

    RCalc	calc;

    calc.Initialize( i1, w1, h1, i2, w2, h2, Nx, Ny );

// Reorganize valid entries of rslt image so that (dx,dy)=(0,0)
// is at the image center and (cx,cy)+(dx,dy) indexes all pixels
// of any sign.

    R.resize( nR );
    A.resize( nR );

    double	vmin = 1e7, vmax = -1e7;

    for( int y = -nnegy; y < nposy; ++y ) {

        int	iy = Nx * (y >= 0 ? y : Ny + y);

        for( int x = -nnegx; x < nposx; ++x ) {

            int	ix = (x >= 0 ? x : Nx + x);
            int	ir = cx+x + wR*(cy+y);

            A[ir] = calc.Valid( LegalRgn, arglr,
                        LegalCnt, arglc, x, y );

            R[ir] = calc.CalcR( rslt[ix+iy] );

            if( R[ir] < vmin )
                vmin = R[ir];

            if( R[ir] > vmax )
                vmax = R[ir];
        }
    }

    if( dbgCor )
        fprintf( flog, "Corr: Center = (%d %d).\n", cx, cy );

    if( verbose ) {
        fprintf( flog,
        "Corr: R image range [%11.6f %11.6f].\n", vmin, vmax );
    }

    MaskA( A, Ox, Oy, Rx, Ry );

// Create S =========================================

    for( int i = 0; i < M; ++i ) {

        double	mag = abs( fft1[i] );

        if( mag > 1e-10 )
            fft1[i] /= sqrt( mag );
    }

    IFT_2D( rslt, fft1, Nx, Ny, flog );

    S.resize( nR );

    vmin = 1e7, vmax = -1e7;

    for( int y = -nnegy; y < nposy; ++y ) {

        int	iy = Nx * (y >= 0 ? y : Ny + y);

        for( int x = -nnegx; x < nposx; ++x ) {

            int	ix = (x >= 0 ? x : Nx + x);
            int	ir = cx+x + wR*(cy+y);

            S[ir] = calc.CalcS( rslt[ix+iy] );

            if( S[ir] < vmin )
                vmin = S[ir];

            if( S[ir] > vmax )
                vmax = S[ir];
        }
    }

    if( verbose ) {
        fprintf( flog,
        "Corr: S image range [%11.6f %11.6f].\n", vmin, vmax );
    }

// ==================================================

    if( dbgCor ) {

        char	simg[32];

        sprintf( simg, "thmA_%d.tif", _dbg_simgidx );
        CorrThmToTif8( simg, i1, Nx, w1, h1, flog );

        sprintf( simg, "thmB_%d.tif", _dbg_simgidx );
        CorrThmToTif8( simg, i2, Nx, w2, h2, flog );

        sprintf( simg, "corr_A_%d.tif", _dbg_simgidx );
        Raster8ToTif8( simg, &A[0], wR, hR, flog );

        sprintf( simg, "corr_R_%d.tif", _dbg_simgidx );
        RasterDblToTifFlt( simg, &R[0], wR, hR, flog );

        sprintf( simg, "corr_S_%d.tif", _dbg_simgidx );
        RasterDblToTifFlt( simg, &S[0], wR, hR, flog );
    }
}

/* --------------------------------------------------------------- */
/* CCorImg::MakeF ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Some interesting filters...

// 0, -1,  0,
//-1,  5, -1,
// 0, -1,  0};

//-1, -1, -1,
//-1,  9, -1,
//-1, -1, -1};

//-1, -1, -1,
//-1, 10, -1,
//-1, -1, -1};

//-1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1,
//-1, -1, 24, -1, -1,
//-1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1};

//-1, -1, -1, -1, -1,
//-1, -4, -4, -4, -1,
//-1, -4, 48, -4, -1,
//-1, -4, -4, -4, -1,
//-1, -1, -1, -1, -1};

//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, 81, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1,
//-1, -1, -1, -1, -1, -1, -1, -1, -1};

void CCorImg::MakeF(
    vector<double>			&F,
    vector<uint8>			&A,
    const vector<double>	&R )
{
    vector<CD>	kfft;
    double		K[] = {
                -1, -1, -1,
                -1,  9, -1,
                -1, -1, -1};
    double		vmax;
    int			ksize	= (int)sqrt( sizeof(K) / sizeof(double) );
    int			grd		= ksize / 2;

    Convolve( F, R, wR, hR, K, ksize, ksize, true, true, kfft, flog );

// Zero invalid border of width = grd (guard band)

    if( grd < 2 )
        grd = 2;

    // top & bottom
    {
        int	nz = grd * wR;

        memset( &A[0], 0, nz * sizeof(uint8) );
        memset( &A[nR - nz], 0, nz * sizeof(uint8) );

        memset( &F[0], 0, nz * sizeof(double) );
        memset( &F[nR - nz], 0, nz * sizeof(double) );
    }

    // left & right
    for( int y = grd; y < hR - grd; ++y ) {

        for( int x = 0; x < grd; ++x ) {

            int	i = x + wR*y;

            A[i] = 0;
            F[i] = 0.0;

            i = wR - 1 - x + wR*y;

            A[i] = 0;
            F[i] = 0.0;
        }
    }

// Normalize the filtered image, keep positive only

    vmax = 1e-7;

    for( int i = 0; i < nR; ++i ) {

        if( F[i] > vmax )
            vmax = F[i];
    }

    for( int i = 0; i < nR; ++i ) {

        if( F[i] > 0.0 )
            F[i] /= vmax;
        else
            F[i] = 0.0;
    }

    if( verbose ) {
        fprintf( flog,
        "Corr: F image range [%11.6f %11.6f] before renorm.\n",
        0.0, vmax );
    }

    if( dbgCor ) {

        char	simg[32];

        sprintf( simg, "corr_F_%d.tif", _dbg_simgidx );
        RasterDblToTifFlt( simg, &F[0], wR, hR, flog );
    }
}

/* --------------------------------------------------------------- */
/* CCorImg::OrderF ----------------------------------------------- */
/* --------------------------------------------------------------- */

class CSort_F_Dec {
public:
    const vector<double>	&F;
public:
    CSort_F_Dec( const vector<double> &F )
        : F(F) {};
    bool operator() ( int a, int b )
        {return F[a] > F[b];};
};


bool CCorImg::OrderF(
    vector<int>				&forder,
    const vector<double>	&F,
    const vector<uint8>		&A,
    const vector<double>	&R,
    double					mincor )
{
// List F pixels in A

    forder.reserve( nR );

    if( mincor > 0.0 ) {

        for( int i = 0; i < nR; ++i ) {
            if( A[i] && R[i] > mincor )
                forder.push_back( i );
        }
    }
    else {

        for( int i = 0; i < nR; ++i ) {
            if( A[i] )
                forder.push_back( i );
        }
    }

// Check non-empty

    if( !forder.size() ) {

        if( verbose )
            fprintf( flog, "Corr: No peak candidates.\n" );

        return false;
    }

// Sort by decreasing F

    CSort_F_Dec	sorter( F );

    sort( forder.begin(), forder.end(), sorter );

    if( dbgCor ) {

        int	ndbg = forder.size();
        if( ndbg > 40 )
            ndbg = 40;

        fprintf( flog, "top %d by F:\n", ndbg );

        for( int i = 0; i < ndbg; ++i ) {
            int	k	= forder[i];
            int	y0	= k / wR;
            int	x0	= k - wR * y0;
            fprintf( flog, "F %.3f R %.3f (x,y)=(%d,%d)\n",
                F[k], R[k], x0, y0 );
        }
    }

    return true;
}

/* --------------------------------------------------------------- */
/* CCorImg::FPeak ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Select highest F that is well isolated from neighbors.
//
// Return true if success.
//
#if 1
bool CCorImg::FPeak(
    int						&rx,
    int						&ry,
    const vector<int>		&forder,
    const vector<double>	&F,
    const vector<double>	&R,
    double					nbmaxht )
{
/* ---------------------------------------------------- */
/* For efficiency, mask innermost peak areas as visited */
/* ---------------------------------------------------- */

    vector<uint8>	mask( nR, 0 );

/* --------------------------------------- */
/* For each peak candidate (highest first) */
/* --------------------------------------- */

    int	nO = forder.size();

    for( int i = 0; i < nO; ++i ) {

        int		ri	= BIG,
                k	= forder[i],
                ro, rm;

        if( mask[k] )
            continue;

        /* ----------------------- */
        /* Candidate {rx, ry, fpk} */
        /* ----------------------- */

        double	fpk	= F[k],
                tol = nbmaxht * fpk;

        ry	= k / wR;
        rx	= k - wR * ry;

        /* ------------------------------------------- */
        /* Scan for inner radius ri (pk ht down by 2x) */
        /* ------------------------------------------- */

        // We select the minimum ri from each of four major
        // scan directions. Originally we used the average
        // over the four directions and this had a problem
        // when the peak (even though compact) rides atop an
        // extended ridge (from a fold) or extended patch
        // (due to resin). In such cases the rate of falloff
        // could be very slow along the ridge, biasing the
        // peak size estimate.

//up:
        for( int y = ry - 1; y >= 0; --y ) {
            if( F[rx + wR*y] <= 0.5*fpk ) {
                ri = min( ri, ry - y );
                goto down;
            }
        }

        goto ri_close;

down:
        for( int y = ry + 1; y < hR; ++y ) {
            if( F[rx + wR*y] <= 0.5*fpk ) {
                ri = min( ri, y - ry );
                goto left;
            }
        }

        goto ri_close;

left:
        for( int x = rx - 1; x >= 0; --x ) {
            if( F[x + wR*ry] <= 0.5*fpk ) {
                ri = min( ri, rx - x );
                goto right;
            }
        }

        goto ri_close;

right:
        for( int x = rx + 1; x < wR; ++x ) {
            if( F[x + wR*ry] <= 0.5*fpk ) {
                ri = min( ri, x - rx );
                goto set_ri;
            }
        }

ri_close:
        if( verbose ) {
            fprintf( flog,
            "Ri_edge F=%.3f R=%.3f (%4d,%4d)\n",
            fpk, R[k], rx, ry );
        }

        continue;

set_ri:
        // noise
        if( ri <= 1 && R[k] < 0.1 ) {

            if( verbose ) {
                fprintf( flog,
                "Micropk F=%.3f R=%.3f (%4d,%4d)\n",
                fpk, R[k], rx, ry );
            }

            continue;
        }

        /* ----------------------- */
        /* Set outer guard band ro */
        /* ----------------------- */

        ro = 3 * ri;

        if( ry - ro < 0 || ry + ro >= hR ||
            rx - ro < 0 || rx + ro >= wR ) {

            if( verbose ) {
                fprintf( flog,
                "Ro_edge F=%.3f R=%.3f (%4d,%4d):"
                " (ri,ro)=(%4d,%4d)\n",
                fpk, R[k], rx, ry,
                ri, ro );
            }

            continue;
        }

        /* ------------------------- */
        /* Mask out core peak region */
        /* ------------------------- */

        if( (rm = ri/3) < 1 )
            rm = 1;

        for( int y = ry - rm; y <= ry + rm; ++y ) {
            for( int x = rx - rm; x <= rx + rm; ++x )
                mask[x + wR*y] = 1;
        }

        /* ----------------------------------------- */
        /* Any pixel > tol in region between ri, ro? */
        /* ----------------------------------------- */

        for( int y = ry - ro; y <= ry + ro; ++y ) {

            for( int x = rx - ro; x <= rx + ro; ++x ) {

                // ignore inner region
                if( y >= ry - ri && y <= ry + ri &&
                    x >= rx - ri && x <= rx + ri ) {

                    continue;
                }

                if( F[x + wR*y] > tol ) {

                    if( verbose ) {

                        fprintf( flog,
                        "Hi_neib F=%.3f R=%.3f (%4d,%4d):"
                        " (ri,ro)=(%4d,%4d)"
                        " neib%%=%.3f @ (%4d,%4d)\n",
                        fpk, R[k], rx, ry,
                        ri, ro,
                        F[x + wR*y]/fpk, x, y );
                    }

                    goto next_i;
                }
            }
        }

        // tested all - good.
        return true;

next_i:;
    }

    if( verbose )
        fprintf( flog, "Corr: All peaks rejected.\n" );

    return false;
}
#endif

#if 0
bool CCorImg::FPeak(
    int						&rx,
    int						&ry,
    const vector<int>		&forder,
    const vector<double>	&F,
    const vector<double>	&R,
    double					nbmaxht )
{
/* ---------------------------------------------------- */
/* For efficiency, mask innermost peak areas as visited */
/* ---------------------------------------------------- */

    vector<uint8>	mask( nR, 0 );

/* --------------------------------------- */
/* For each peak candidate (highest first) */
/* --------------------------------------- */

    int	nO = forder.size();

    for( int i = 0; i < nO; ++i ) {

        int		ri	= 0,
                k	= forder[i],
                ro, rm, bl, br, bb, bt;

        if( mask[k] )
            continue;

        /* ----------------------- */
        /* Candidate {rx, ry, fpk} */
        /* ----------------------- */

        double	fpk	= F[k],
                tol = nbmaxht * fpk;

        ry	= k / wR;
        rx	= k - wR * ry;

        /* ------------------------------------------- */
        /* Scan for inner radius ri (pk ht down by 2x) */
        /* ------------------------------------------- */

up:
        for( int y = ry - 1; y >= 0; --y ) {
            if( F[rx + wR*y] <= 0.5*fpk ) {
                ri += ry - y;
                goto down;
            }
        }

        continue;

down:
        for( int y = ry + 1; y < hR; ++y ) {
            if( F[rx + wR*y] <= 0.5*fpk ) {
                ri += y - ry;
                goto left;
            }
        }

        continue;

left:
        for( int x = rx - 1; x >= 0; --x ) {
            if( F[x + wR*ry] <= 0.5*fpk ) {
                ri += rx - x;
                goto right;
            }
        }

        continue;

right:
        for( int x = rx + 1; x < wR; ++x ) {
            if( F[x + wR*ry] <= 0.5*fpk ) {
                ri += x - rx;
                goto set_ri;
            }
        }

        continue;

set_ri:
        ri /= 4;

        // noise
        if( ri <= 1 && R[k] < 0.1 )
            continue;

        /* ------------------------- */
        /* Mask out core peak region */
        /* ------------------------- */

        if( (rm = ri/3) < 1 )
            rm = 1;

        bb = max( ry - rm, 0 );
        bt = min( ry + rm, hR - 1 );
        bl = max( rx - rm, 0 );
        br = min( rx + rm, wR - 1 );

        for( int y = bb; y <= bt; ++y ) {
            for( int x = bl; x <= br; ++x )
                mask[x + wR*y] = 1;
        }

        /* ----------------------- */
        /* Set outer guard band ro */
        /* ----------------------- */

        ro = 3 * ri;

        bb = max( ry - ro, 0 );
        bt = min( ry + ro, hR - 1 );
        bl = max( rx - ro, 0 );
        br = min( rx + ro, wR - 1 );

        /* ----------------------------------------- */
        /* Any pixel > tol in region between ri, ro? */
        /* ----------------------------------------- */

        for( int y = bb; y <= bt; ++y ) {

            for( int x = bl; x <= br; ++x ) {

                // ignore inner region
                if( y >= ry - ri && y <= ry + ri &&
                    x >= rx - ri && x <= rx + ri ) {

                    continue;
                }

                if( F[x + wR*y] > tol ) {

                    if( verbose ) {

                        fprintf( flog,
                        "Reject F=%.3f R=%.3f (%4d,%4d):"
                        " (ri,ro)=(%4d,%4d)"
                        " neib%%=%.3f @ (%4d,%4d)\n",
                        fpk, R[k], rx, ry,
                        ri, ro,
                        F[x + wR*y]/fpk, x, y );
                    }

                    goto next_i;
                }
            }
        }

        // tested all - good.
        return true;

next_i:;
    }

    if( verbose )
        fprintf( flog, "Corr: All peaks rejected.\n" );

    return false;
}
#endif

/* --------------------------------------------------------------- */
/* CCorImg::SimpleMax -------------------------------------------- */
/* --------------------------------------------------------------- */

void CCorImg::SimpleMax(
    int						&rx,
    int						&ry,
    const vector<int>		&forder )
{
    int	k = forder[0];

    ry	= k / wR;
    rx	= k - wR * ry;
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
    int				w,
    int				h )
{
    int		ix = (int)xpk,
            iy = (int)ypk;
    double	z1 = I[ix + w*iy],
            z0,
            z2;

    if( ix - d >= 0 &&
        ix + d <  w &&
        (z0 = I[ix-d + w*iy]) > 0.0 &&
        (z2 = I[ix+d + w*iy]) > 0.0 ) {

        z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

        if( fabs( z0 ) <= 8 * d )
            xpk += z0;
    }

    if( iy - d >= 0 &&
        iy + d <  h &&
        (z0 = I[ix + w*(iy-d)]) > 0.0 &&
        (z2 = I[ix + w*(iy+d)]) > 0.0 ) {

        z0 = d * (z0 - z2) / (2 * (z0 + z2 - z1 - z1));

        if( fabs( z0 ) <= 8 * d )
            ypk += z0;
    }
}

/* --------------------------------------------------------------- */
/* CCorImg::ReturnR ---------------------------------------------- */
/* --------------------------------------------------------------- */

double CCorImg::ReturnR(
    double					&dx,
    double					&dy,
    int						rx,
    int						ry,
    const vector<double>	&Rreport,
    const vector<double>	&Rtweak )
{
    double	bigR = Rreport[rx + wR*ry];

    if( verbose ) {
        fprintf( flog,
        "Corr: Max corr %f at [dx,dy]=[%d,%d] (x,y)=(%d,%d).\n",
        bigR, rx - cx, ry - cy, rx, ry );
    }

    dx = rx;
    dy = ry;
    ParabPeak( dx, dy, 1, &Rtweak[0], wR, hR );

    dx += B2.L - B1.L - cx;
    dy += B2.B - B1.B - cy;

    if( dbgCor )
         ++_dbg_simgidx;

    return bigR;
}

/* --------------------------------------------------------------- */
/* CorrImagesF --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// mincor: If non-zero, used to prescreen F values during the
// peak-hunting phase.
//
// nbmaxht: A candidate F-peak is rejected if its guard band
// contains another pixel with F > nbmaxht*peak.
//
// {Ox,Oy,Rx,Ry}: If Rx > 0 and Ry > 0, search narrowed to oval
// with origin (Ox,Oy) and semimajor axes Rx,Ry.
//
// fft2: Cache of image2 FFT. On entry, if fft2 has the
// correct size it is used. Otherwise recomputed here.
//
// Version using F and well isolated F peak.
//
double CorrImagesF(
    FILE					*flog,
    int						verbose,
    double					&dx,
    double					&dy,
    const vector<Point>		&ip1,
    const vector<double>	&iv1,
    const vector<Point>		&ip2,
    const vector<double>	&iv2,
    EvalType				LegalRgn,
    void*					arglr,
    EvalType				LegalCnt,
    void*					arglc,
    double					mincor,
    double					nbmaxht,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry,
    vector<CD>				&fft2 )
{
    CCorImg			cc;
    vector<double>	R;
    vector<uint8>	A;
    vector<double>	F;
    vector<int>		forder;
    int				rx;
    int				ry;

    if( dbgCor )
        verbose = true;

    if( !cc.SetDims( flog, verbose, ip1, ip2 ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    cc.MakeRandA( R, A, ip1, iv1, ip2, iv2,
        LegalRgn, arglr, LegalCnt, arglc,
        Ox, Oy, Rx, Ry, fft2 );

    cc.MakeF( F, A, R );

    if( !cc.OrderF( forder, F, A, R, mincor ) ||
        !cc.FPeak( rx, ry, forder, F, R, nbmaxht ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    return cc.ReturnR( dx, dy, rx, ry, R, R );
}

/* --------------------------------------------------------------- */
/* CorrImagesR --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// mincor: If non-zero, used to prescreen R values during the
// peak-hunting phase.
//
// nbmaxht: Dummy slot for compatibility with F-version.
//
// {Ox,Oy,Rx,Ry}: If Rx > 0 and Ry > 0, search narrowed to oval
// with origin (Ox,Oy) and semimajor axes Rx,Ry.
//
// fft2: Cache of image2 FFT. On entry, if fft2 has the
// correct size it is used. Otherwise recomputed here.
//
// Version using straight max R.
//
double CorrImagesR(
    FILE					*flog,
    int						verbose,
    double					&dx,
    double					&dy,
    const vector<Point>		&ip1,
    const vector<double>	&iv1,
    const vector<Point>		&ip2,
    const vector<double>	&iv2,
    EvalType				LegalRgn,
    void*					arglr,
    EvalType				LegalCnt,
    void*					arglc,
    double					mincor,
    double					nbmaxht,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry,
    vector<CD>				&fft2 )
{
    CCorImg			cc;
    vector<double>	R;
    vector<uint8>	A;
    vector<int>		order;
    int				rx;
    int				ry;

    if( dbgCor )
        verbose = true;

    if( !cc.SetDims( flog, verbose, ip1, ip2 ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    cc.MakeRandA( R, A, ip1, iv1, ip2, iv2,
        LegalRgn, arglr, LegalCnt, arglc,
        Ox, Oy, Rx, Ry, fft2 );

    // actually orders R, here
    if( !cc.OrderF( order, R, A, R, mincor ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    cc.SimpleMax( rx, ry, order );

    return cc.ReturnR( dx, dy, rx, ry, R, R );
}

/* --------------------------------------------------------------- */
/* CorrImagesS --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return cross-correlation and additive displacement
// (dx, dy) that places set ip1 into bounding box of ip2.
//
// mincor: If non-zero, used to prescreen S values during the
// peak-hunting phase.
//
// nbmaxht: Dummy slot for compatibility with F-version.
//
// {Ox,Oy,Rx,Ry}: If Rx > 0 and Ry > 0, search narrowed to oval
// with origin (Ox,Oy) and semimajor axes Rx,Ry.
//
// fft2: Cache of image2 FFT. On entry, if fft2 has the
// correct size it is used. Otherwise recomputed here.
//
// Version using FFT power spectrum filtering as prescribed
// by Art Wetzel, followed by straight max S.
//
double CorrImagesS(
    FILE					*flog,
    int						verbose,
    double					&dx,
    double					&dy,
    const vector<Point>		&ip1,
    const vector<double>	&iv1,
    const vector<Point>		&ip2,
    const vector<double>	&iv2,
    EvalType				LegalRgn,
    void*					arglr,
    EvalType				LegalCnt,
    void*					arglc,
    double					mincor,
    double					nbmaxht,
    int						Ox,
    int						Oy,
    int						Rx,
    int						Ry,
    vector<CD>				&fft2 )
{
    CCorImg			cc;
    vector<double>	R;
    vector<uint8>	A;
    vector<double>	S;
    vector<int>		order;
    int				rx;
    int				ry;

    if( dbgCor )
        verbose = true;

    if( !cc.SetDims( flog, verbose, ip1, ip2 ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    cc.MakeSandRandA( S, R, A, ip1, iv1, ip2, iv2,
        LegalRgn, arglr, LegalCnt, arglc,
        Ox, Oy, Rx, Ry, fft2 );

    // actually orders S, here
    if( !cc.OrderF( order, S, A, R, mincor ) ) {

        dx	= 0.0;
        dy	= 0.0;

        return 0.0;
    }

    cc.SimpleMax( rx, ry, order );

// S is really better than R for comparing different cases,
// say sweeps for strips or blocks. So we return S directly.
// To scale S roughly like R, we create the following hack:
// To get better content independence, form snr = (r-mn)/sd.
// Then, to account somewhat for image size, divide by linear
// size estimator sqrt( area ). Finally, arbitrary factor 10
// scales result "roughly" into range [-1, 1].

    double	mn, sd, r = cc.ReturnR( dx, dy, rx, ry, S, S );

    Stats( S, mn, sd );

    return 10.0 * (r - mn) / (sd * sqrt( S.size() ));
}


