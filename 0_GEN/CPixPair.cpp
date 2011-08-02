

#include	"CPixPair.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Timer.h"






/* --------------------------------------------------------------- */
/* PixPair::DoGSetFilter ----------------------------------------- */
/* --------------------------------------------------------------- */

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
void PixPair::DoGSetFilter( DoGdata &D, int r1, int r2, FILE* flog )
{
	int	N, fltsz = 3 * r2;

// Set FFT sizing

	D.N = N = CeilPow2( max( wf, hf ) + fltsz );
	D.M = N * (N/2 + 1);

	fprintf( flog, "DoGSetFilter: Using size %d.\n", N );

// Integral normalizations
//
// Normalization sums should ideally equal 2pi. They will differ
// slightly due to quantization and the finite tail length.

	double	sum1 = 0.0, sum2 = 0.0;

	D.I.resize( N * N, 0.0 );

	for( int x = -fltsz; x <= fltsz; ++x ) {

		for( int y = -fltsz; y <= fltsz; ++y ) {

			double	R = x*x + y*y;

			sum1 += exp( -( R/(2*r1*r1) ) );
			sum2 += exp( -( R/(2*r2*r2) ) );
		}
	}

	fprintf( flog,
	"DoGSetFilter: Check %f equals %f equals 2pi.\n",
	sum1/(r1*r1), sum2/(r2*r2) );

// Set filter values in workspace

	for( int x = -fltsz; x <= fltsz; ++x ) {

		for( int y = -fltsz; y <= fltsz; ++y ) {

			double	R	= x*x + y*y;
			double	v1	= exp( -( R/(2*r1*r1) ) );
			double	v2	= exp( -( R/(2*r2*r2) ) );
			int		ix	= (x >= 0 ? x : N + x);
			int		iy	= (y >= 0 ? y : N + y);

			D.I[ix + N*iy] = v1/sum1 - v2/sum2;
		}
	}

// FFT of filter

	FFT_2D( D.flt, D.I, N, N );
}

/* --------------------------------------------------------------- */
/* PixPair::DoGApply --------------------------------------------- */
/* --------------------------------------------------------------- */

// Assumes values are already normalized on entry.
//
void PixPair::DoGApply(
	vector<double>			&dst,
	const vector<double>	&src,
	DoGdata					&D )
{
// Copy values to workspace then FFT

	vector<CD>	fft;
	int			N = D.N;

	memset( &D.I[0], 0, N * N * sizeof(double) );
	CopyRaster( &D.I[0], N, &src[0], wf, wf, hf );

	FFT_2D( fft, D.I, N, N );

// Convolve with filter

	for( int i = 0; i < D.M; ++i )
		fft[i] *= D.flt[i];

// FFT back

	IFT_2D( D.I, fft, N, N );

// Copy to dst and normalize

	dst.resize( wf * hf );

	CopyRaster( &dst[0], wf, &D.I[0], N, wf, hf );

	Normalize( dst );
}

/* --------------------------------------------------------------- */
/* PixPair::Downsample ------------------------------------------- */
/* --------------------------------------------------------------- */

void PixPair::Downsample(
	vector<double>			&dst,
	const vector<double>	&src )
{
	int	n = scl * scl;

	dst.resize( ws * hs );

	for( int iy = 0; iy < hs; ++iy ) {

		for( int ix = 0; ix < ws; ++ix ) {

			double	sum = 0.0;

			for( int dy = 0; dy < scl; ++dy ) {

				for( int dx = 0; dx < scl; ++dx )
					sum += src[ix*scl+dx + wf*(iy*scl+dy)];
			}

			dst[ix+ws*iy] = sum / n;
		}
	}
}

/* --------------------------------------------------------------- */
/* PixPair::Load ------------------------------------------------- */
/* --------------------------------------------------------------- */

bool PixPair::Load(
	const char	*apath,
	const char	*bpath,
	int			order,
	int			bDoG,
	int			r1,
	int			r2,
	FILE*		flog )
{
	printf( "\n---- Image loading ----\n" );

	clock_t		t0 = StartTiming();

/* ----------------------------- */
/* Load and sanity check rasters */
/* ----------------------------- */

	uint8	*aras, *bras;
	uint32	wa, ha, wb, hb;
	int		ok = false;

	aras = Raster8FromAny( apath, wa, ha, flog, false );
	bras = Raster8FromAny( bpath, wb, hb, flog, false );

	if( !aras || !bras ) {
		fprintf( flog,
		"PixPair: Picture load failure.\n" );
		goto exit;
	}

	if( wa != wb || ha != hb ) {
		fprintf( flog,
		"PixPair: Nonmatching picture dimensions.\n" );
		goto exit;
	}

	ok		= true;
	wf		= wa;
	hf		= ha;
	ws		= wa;
	hs		= ha;
	scl		= 1;

/* ------- */
/* Flatten */
/* ------- */

	LegPolyFlatten( _avf, aras, wf, hf, order );
	RasterFree( aras );
	aras = NULL;

	LegPolyFlatten( _bvf, bras, wf, hf, order );
	RasterFree( bras );
	bras = NULL;

	avs_vfy	= avs_aln = avf_vfy	= avf_aln = &_avf;
	bvs_vfy	= bvs_aln = bvf_vfy	= bvf_aln = &_bvf;

/* ------------- */
/* Apply filters */
/* ------------- */

	if( bDoG ) {

		DoGdata		D;

		DoGSetFilter( D, r1, r2, flog );
		DoGApply( _avfflt, _avf, D );
		DoGApply( _bvfflt, _bvf, D );

		avs_aln = avf_aln = &_avfflt;
		bvs_aln = bvf_aln = &_bvfflt;
	}

/* --------------------- */
/* Downsample all images */
/* --------------------- */

	if( ws > 2048 || hs >= 2048 ) {

		do {
			ws		/= 2;
			hs		/= 2;
			scl		*= 2;
		} while( ws > 2048 || hs > 2048 );

		fprintf( flog, "PixPair: Scaling by %d.\n", scl );

		if( ws * scl != wf || hs * scl != hf ) {

			fprintf( flog,
			"PixPair: Dimensions not multiple of scale!\n" );
			goto exit;
		}

		Downsample( _avs, _avf );
		Downsample( _bvs, _bvf );

		avs_vfy = avs_aln = &_avs;
		bvs_vfy = bvs_aln = &_bvs;

		if( bDoG ) {

			Downsample( _avsflt, _avfflt );
			Downsample( _bvsflt, _bvfflt );

			avs_aln = &_avsflt;
			bvs_aln = &_bvsflt;
		}
	}
	else
		fprintf( flog, "PixPair: Using image scale=1.\n" );

/* ------------------------------ */
/* Write DoG images for debugging */
/* ------------------------------ */

#if 0
	if( bDoG ) {
		VectorDblToTif8( "DoGa.tif", avs_aln, ws, hs );
		VectorDblToTif8( "DoGb.tif", bvs_aln, ws, hs );
	}
#endif

/* -------- */
/* Clean up */
/* -------- */

exit:
	if( aras )
		RasterFree( aras );

	if( bras )
		RasterFree( bras );

	StopTiming( flog, "Image conditioning", t0 );

	return ok;
}


