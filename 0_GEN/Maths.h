

#pragma once


#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stdio.h>

#include	<vector>
using namespace std;


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
	int  HowMany()				{return n;};
	void Stats( double &avg, double &std );
};


void Stats( const vector<double> &v, double &avg, double &std );

void StatsRasterNonZeros(
	const uint8*	raster,
	int				npix,
	double			&avg,
	double			&std );

void Normalize( vector<double> &v );
void Normalize( double *a, int n );

void NormalizeNonZeros( vector<double> &v );

void CoNormalize( vector<double> &v1, vector<double> &v2 );

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

void DistributePixel(
	double			x,
	double			y,
	double			val,
	float*			image,
	int				w );

inline void DistributePixel(
	double			x,
	double			y,
	double			val,
	vector<double>	&image,
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

double BiCubicInterp( const double* image, int w, Point p );
double BiCubicInterp( const uint8* image, int w, Point p );

void ImageGradients(
	double				&slopex,
	double				&slopey,
	const vector<uint8>	&image,
	int					w,
	const vector<Point>	&pts,
	FILE				*flog );

bool IsLowContrast( vector<double> &I, double std );


