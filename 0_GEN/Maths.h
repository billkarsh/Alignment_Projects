

#pragma once


#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Inlines ------------------------------------------------------- */
/* --------------------------------------------------------------- */

inline int iabs( int x )		{return x >= 0 ? x : -x;};
inline int min( int x, int y )	{return x < y ? x : y;};
inline int max( int x, int y )	{return x > y ? x : y;};
inline int ROUND( double x )	{return x >= 0.0 ? int(x+0.5) : int(x-0.5);}
inline int RND( double x )		{return x >  0.0 ? int(x+0.5) : 0;}

inline double sgn( double a )	{return a > 0.0 ? 1.0 : (a < 0.0 ? -1.0 : 0.0);};

/* --------------------------------------------------------------- */
/* Hash ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

uint32 SuperFastHash( const char * data, int len );

/* --------------------------------------------------------------- */
/* Numerics ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int CeilPow2( int n );

/* --------------------------------------------------------------- */
/* Statistics ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class MeanStd {	// for computing mean and standard deviation

private:
    double	sum, sum2;
    int		n;

public:
    MeanStd()					{Reset();};
    void Reset()				{sum = sum2 = 0.0; n=0;};
    void Element( double a )	{sum += a; sum2 += a*a; ++n;};
    void Run( double a, int k ) {sum += a*k; sum2 += a*a*k; n += k;};
    int  HowMany()				{return n;};
    void Stats( double &avg, double &std );
};


void Stats( const vector<double> &v, double &avg, double &std );
void Stats( const uint8 *v, int N, double &avg, double &std );

void StatsRasterNonZeros(
    const uint8*	raster,
    int				npix,
    double			&avg,
    double			&std );

double Normalize( vector<double> &v );
double Normalize( double *a, int n );

double NormalizeNonZeros( vector<double> &v );

double CoNormalize( vector<double> &v1, vector<double> &v2 );

void CoUpperTail(
    vector<double>	&v1,
    vector<double>	&v2,
    double			nsigma );

void CoExcludeMiddle(
    vector<double>	&v1,
    vector<double>	&v2,
    double			nsigma );

int FirstNonzero( const double *v, int n );

int IndexOfMaxVal( const double *v, int n );

double MedianVal( vector<double> &V );

void LineFit(
    double			*icpt,
    double			*slope,
    double			*lincor,
    const double	*x,
    const double	*y,
    int				i0,
    int				ilim );

void Histogram(
    double			&uflo,
    double			&oflo,
    double			*bins,
    int				nbins,
    double			datamin,
    double			datamax,
    const uint16	*data,
    int				ndata,
    bool			breset );

int PercentileBin(
    const double	*bins,
    int				nbins,
    double			frac );

double IsoDataThresh(
    const double	*bins,
    int				nbins,
    double			datamin,
    double			datamax );

double MinSepThresh(
    const double	*bins,
    int				nbins,
    double			datamin,
    double			datamax );

double OtsuThresh(
    const double	*bins,
    int				nbins,
    double			datamin,
    double			datamax );

/* --------------------------------------------------------------- */
/* Matrices ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void PrintVectorAsMat( FILE* f, const vector<double> &v, int w );

void Print3x3Matrix( FILE* f, const double a[3][3] );

double Invert3x3Rowlist( double i[9], const double a[9] );

double Invert3x3Matrix( double i[3][3], const double a[3][3] );
double Invert4x4Matrix( double i[4][4], const double a[4][4] );

/* --------------------------------------------------------------- */
/* Legendre polys ------------------------------------------------ */
/* --------------------------------------------------------------- */

void LegPolyCreate(
    vector<vector<double> >	&L,
    int						maxOrder,
    int						nSamples );

void LegPolyFlatten(
    vector<double>		&vals,
    const uint8*		src,
    int					w,
    int					h,
    int					maxOrder );

void LegPolyFlatten(
    vector<double>		&vals,
    const uint16*		src,
    int					w,
    int					h,
    int					maxOrder,
    int					offset );

void LegPolyFlatten(
    vector<double>		&vals,
    const vector<Point>	&pts,
    const uint8*		src,
    int					w,
    int					h,
    int					maxOrder );

/* --------------------------------------------------------------- */
/* Images -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CopyRaster(
    double*			D,
    int				wD,
    const double*	S,
    int				wS,
    int				wCopy,
    int				hCopy );

void DecimateVector(
    vector<Point>	&p,
    vector<double>	&v,
    int				w,
    int				h,
    int				lbin );

double InterpolatePixel(
    double					x,
    double					y,
    const uint8*			image,
    int						w );

double InterpolatePixel(
    double					x,
    double					y,
    const vector<double>	&image,
    int						w );

double SafeInterp(
    double					x,
    double					y,
    const uint8*			image,
    int						w,
    int						h );

double SafeInterp(
    double					x,
    double					y,
    const uint16*			image,
    int						w,
    int						h );

double SafeInterp(
    double					x,
    double					y,
    const double*			image,
    int						w,
    int						h );

void DistributePixel(
    double			x,
    double			y,
    double			val,
    float*			image,
    int				w,
    int				h );

inline void DistributePixel(
    double			x,
    double			y,
    double			val,
    vector<double>	&image,
    int				w,
    int				h )
{
    int		ix		= (int)floor( x );
    int		iy		= (int)floor( y );
    int		i		= ix + w * iy;
    double	alpha	= x - ix;
    double	beta	= y - iy;

    if( ix >= 0 && ix < w ) {

        if( iy >= 0 && iy < h )
            image[i] += (1.0-alpha)*(1.0-beta)*val;

        ++iy;

        if( iy >= 0 && iy < h )
            image[i + w] += (1.0-alpha)*beta*val;

        --iy;
    }

    ++ix;
    ++i;

    if( ix >= 0 && ix < w ) {

        if( iy >= 0 && iy < h )
            image[i] += alpha*(1.0-beta)*val;

        ++iy;

        if( iy >= 0 && iy < h )
            image[i + w] += alpha*beta*val;
    }
}

double BiCubicInterp( const double* image, int w, Point p );
double BiCubicInterp( const uint8* image, int w, Point p );

void ImageGradients(
    double				&slopex,
    double				&slopey,
    const vector<uint8>	&image,
    int					w,
    const vector<Point>	&pts,
    FILE				*flog );

bool IsLowContrast(
    const vector<double>	&I,
    double					std,
    FILE					*flog = stdout );

void EmbedExtended8(
    vector<uint8>	&dst,
    const uint8		*src,
    int				w,
    int				h,
    int				r );

void ExtractEmbedded8(
    uint8			*dst,
    const uint8		*src,
    int				w,
    int				h,
    int				r );

void Downsample8(
    vector<uint8>	&dst,
    const uint8		*src,
    int				&w,
    int				&h,
    int				scl );

void Upsize8(
    vector<uint8>	&dst,
    const uint8		*src,
    int				wd,
    int				hd,
    int				ws,
    int				hs );

void Sobel8(
    uint8			*dst,
    const uint8		*src,
    int				w,
    int				h );

void Median8(
    uint8			*dst,
    const uint8		*src,
    int				w,
    int				h,
    int				r );

void ResinMask8(
    vector<uint8>	&msk,
    const uint8		*src,
    int				w,
    int				h,
    bool			samelayer );


