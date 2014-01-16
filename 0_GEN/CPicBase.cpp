

#include	"CPicBase.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"

#include	<stdlib.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* PicBase constructor ------------------------------------------- */
/* --------------------------------------------------------------- */

PicBase::PicBase()
{
	raster		= NULL;
	original	= NULL;
	external	= NULL;
	z			= 0;
	scale		= 1;
}


PicBase::~PicBase()
{
	if( raster == original ) {

		if( raster && raster != external )
			RasterFree( raster );
	}
	else {

		if( raster && raster != external )
			RasterFree( raster );

		if( original && original != external )
			RasterFree( original );
	}
}

/* --------------------------------------------------------------- */
/* PicBase::LoadOriginal ----------------------------------------- */
/* --------------------------------------------------------------- */

void PicBase::LoadOriginal(
	const char*	name,
	FILE*		flog,
	bool		transpose )
{
	fname = name;

	original =
	Raster8FromAny( name, w, h, flog, transpose );

	CopyOriginal();
}

/* --------------------------------------------------------------- */
/* PicBase::SetExternal ------------------------------------------ */
/* --------------------------------------------------------------- */

void PicBase::SetExternal( const uint8* in_raster, uint32 w, uint32 h )
{
	external	= original = (uint8*)in_raster;
	this->w		= w;
	this->h		= h;

	CopyOriginal();
}

/* --------------------------------------------------------------- */
/* PicBase::CopyOriginal ----------------------------------------- */
/* --------------------------------------------------------------- */

void PicBase::CopyOriginal()
{
	if( raster && raster != original && raster != external )
		RasterFree( raster );

	scale	= 1;
	raster	= original;
}

/* --------------------------------------------------------------- */
/* PicBase::DownsampleIfNeeded ----------------------------------- */
/* --------------------------------------------------------------- */

// If the original is bigger than 2K,
// create a downsampled copy (but keep the original)
//
void PicBase::DownsampleIfNeeded( FILE* flog )
{
	CopyOriginal();

	if( w <= 2048 && h <= 2048 )
		return;

// Find scale factor that will make each dimension <= 2048

	int		origw = w, origh = h;

	do {
		w		/= 2;
		h		/= 2;
		scale	*= 2;
	} while( w > 2048 || h > 2048 );

	fprintf( flog, "DownsampleIfNeeded: Will scale by %d\n", scale );

	if( w * scale != origw || h * scale != origh ) {

		fprintf( flog,
		"DownsampleIfNeeded: Image does not divide evenly!\n" );
		exit( 45 );
	}

	int		npq = scale * scale;

	raster = (uint8*)RasterAlloc( w * h * sizeof(uint8) );

	for( int iy = 0; iy < h; ++iy ) {

		for( int ix = 0; ix < w; ++ix ) {

			int		sum = 0;

			for( int dy = 0; dy < scale; ++dy ) {

				for( int dx = 0; dx < scale; ++dx )
					sum += original[ix*scale+dx + origw*(iy*scale+dy)];
			}

			raster[ix+w*iy] = (sum + npq/2)/npq;	// rounding
		}
	}

	tr.MulXY( 1.0 / scale );
}

/* --------------------------------------------------------------- */
/* PicBase::MakeFFTExist ----------------------------------------- */
/* --------------------------------------------------------------- */

// Calculate whole-frame FFT (fft_of_frame).
//
// Arg (i) is just for printing.
//
void PicBase::MakeFFTExist( int i )
{
	printf(
	"MakeFFTExist: called with %d,"
	" fft_of_frame.size=%ld\n", i, fft_of_frame.size() );

	if( fft_of_frame.size() )	// already exists
		return;

	int N	= 4096;
	int M	= N*(N/2+1);	// count of complex in FFT of 2D real
	int np	= w * h;

// Create vector of all the pixels in the frame,
// so we can normalize their values.

	vector<double>	px( np );
	double			avg, std;

	for( int k = 0; k < np; ++k ) {

		int	y = k / w;
		int	x = k - w*y;
		px[k] = raster[k];
	}

	Stats( px, avg, std );
	printf( "MakeFFTExist: Frame %d: mean %f, std %f\n", i, avg, std );

// Recompute the mean and standard dev. using only 'sensible' pixels

	px.clear();

	for( int k = 0; k < np; ++k ) {

		int		y	= k / w;
		int		x	= k - w*y;
		double	val	= (raster[k] - avg)/std;

		if( fabs( val ) <= 1.5 )
			px.push_back( raster[k] );
	}

	Stats( px, avg, std );
	printf( "MakeFFTExist: Frame %d: mean %f, std %f\n", i, avg, std );

// Now make frame of normalized values.

	vector<double>	fr( N * N, 0.0 );

	for( int k = 0; k < np; ++k ) {

		int		y	= k / w;
		int		x	= k - w*y;
		double	val	= (raster[k] - avg)/std;

		// replace outlying values with 0 (blank).
		// These are normally artifacts producing false matches.

		if( fabs( val ) > 1.5 )
			val = 0.0;

		fr[x + N*y] = val;
	}

	FFT_2D( fft_of_frame, fr, N, N, false );

#if 0
// Now zero out the first few terms.
// These are distractions since we are looking for small slices.
// We should do this differently for X and Y matches, but just do
// the first few for now.

	int	stride	= N/2+1;
	int	NR		= 3;	// remove through this component
	CD	cz( 0.0, 0.0 );

	for( int j = 0; j <= NR; ++j ) {

		for( int k = 0; k <= NR; ++k )
			fft_of_frame[k + j*stride] = cz;
	}
#endif
}

/* --------------------------------------------------------------- */
/* PicBase::MakeDoGExist ----------------------------------------- */
/* --------------------------------------------------------------- */

// Compute a Difference-of-Gaussian raster, if not done previously.
//
// This is an edge enhancement filter. Each Gaussian is a blurring
// filter that passes frequencies longer than its (inverse) standard
// dev (here, specified as a radius). The difference, then, passes
// frequencies between the two (inverse) radii: r1 < r2. The length
// scale of interest is the number of pixels over which there is
// rapid intensity variation, usually, near object edges. This is
// usually a few pixels, so choose {r1, r2} around {3, 6}.
//
// Note: We calc filter out to length = fltsz = 3 sigma = 3 x r2.
//
void PicBase::MakeDoGExist( vector<CD> &filter, int r1, int r2 )
{
	int		npix	= w * h,
			fltsz	= 3 * r2;
	int		i, x, y;

// Early exit?

	if( DoG.size() )
		return;

	DoG.resize( npix );

	if( r2 <= 0.0 ) {
		memcpy( &DoG[0], &raster[0], npix );
		return;
	}

// Set FFT sizing

	int	N = CeilPow2( max( w, h ) + fltsz );

	printf( "MakeDoGExist: Using size %d\n", N );

	int M = N * (N/2 + 1);	// size of array of complex coefficients

// Construct filter

	if( !filter.size() ) {

		vector<double>	image( N * N, 0.0 );
		double			sum1 = 0.0,	// for normalizing.
						sum2 = 0.0;	// each integral should = 2*pi

		// they are not quite the same due to
		// quantization and removal of tails.

		for( x = -fltsz; x <= fltsz; ++x ) {

			for( y = -fltsz; y <= fltsz; ++y ) {

				double	rad2 = x*x + y*y;

				sum1 += exp( -( rad2/(2*r1*r1) ) );
				sum2 += exp( -( rad2/(2*r2*r2) ) );
			}
		}

		printf( "MakeDoGExist: Check %f equals %f equals 2pi.\n",
			sum1/(r1*r1), sum2/(r2*r2) );

		for( x = -fltsz; x <= fltsz; ++x ) {

			for( y = -fltsz; y <= fltsz; ++y ) {

				double	rad2	= x*x + y*y;
				double	v1		= exp( -( rad2/(2*r1*r1) ) );
				double	v2		= exp( -( rad2/(2*r2*r2) ) );
				int		ix		= (x >= 0 ? x : N + x);
				int		iy		= (y >= 0 ? y : N + y);

				image[ix + N*iy] = v1/sum1 - v2/sum2;
			}
		}

		// Now fft it
		FFT_2D( filter, image, N, N, false );
	}

// Normalize image values

	vector<double>	v( npix );

	for( i = 0; i < npix; ++i )
		v[i] = raster[i];

	Normalize( v );

// Copy values to the large image, then FFT

	vector<double>	imag( N * N, 0.0 );
	vector<CD>		freq;

	CopyRaster( &imag[0], N, &v[0], w, w, h );

	FFT_2D( freq, imag, N, N, false );

// Convolve with filter

	for( i = 0; i < M; ++i )
		freq[i] *= filter[i];

// FFT back

	IFT_2D( imag, freq, N, N );

// Copy back to v for normalization

	CopyRaster( &v[0], w, &imag[0], N, w, h );

	Normalize( v );

// Copy v to DoG

	for( i = 0; i < npix; ++i ) {

		int pix	= 127 + int(32 * v[i]);

		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;

		DoG[i] = pix;
	}
}


