

#include	"CPixPair.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Timer.h"






/* --------------------------------------------------------------- */
/* MakeDoGKernel ------------------------------------------------- */
/* --------------------------------------------------------------- */

// This is an edge enhancement filter. Each Gaussian is a blurring
// filter that passes frequencies longer than its (inverse) standard
// dev (here, specified as a radius). The difference, then, passes
// frequencies between the two (inverse) radii: r1 < r2. The length
// scale of interest is the number of pixels over which there is
// rapid intensity variation, usually, near object edges. This is
// usually a few pixels, so choose {r1, r2} around {3, 6}.
//
// Note: Radius of kernel: 3 sigma = 3 x r2.
//
static int MakeDoGKernel(
	vector<double>	&DoG,
	int				r1,
	int				r2,
	FILE*			flog )
{
	int	d = 3 * r2,
		D = 2 * d + 1;

// Normalization sums should ideally equal 2pi. They will differ
// slightly due to quantization and the finite tail length.

	double	sum1 = 0.0, sum2 = 0.0;

	for( int y = -d; y <= d; ++y ) {

		for( int x = -d; x <= d; ++x ) {

			double	R = x*x + y*y;

			sum1 += exp( -( R/(2*r1*r1) ) );
			sum2 += exp( -( R/(2*r2*r2) ) );
		}
	}

	fprintf( flog,
	"MakeDoGKernel: Check %f equals %f equals 2pi.\n",
	sum1/(r1*r1), sum2/(r2*r2) );

// Store values

	DoG.resize( D * D );

	for( int y = -d; y <= d; ++y ) {

		for( int x = -d; x <= d; ++x ) {

			double	R	= x*x + y*y;
			double	v1	= exp( -( R/(2*r1*r1) ) );
			double	v2	= exp( -( R/(2*r2*r2) ) );

			DoG[x+d + D*(y+d)] = v1/sum1 - v2/sum2;
		}
	}

	return D;
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

	LegPolyFlatten( _bvf, bras, wf, hf, order );
	RasterFree( bras );

	avs_vfy	= avs_aln = avf_vfy	= avf_aln = &_avf;
	bvs_vfy	= bvs_aln = bvf_vfy	= bvf_aln = &_bvf;

/* ------------- */
/* Apply filters */
/* ------------- */

	if( bDoG ) {

		vector<double>	DoG;
		vector<CD>		kfft;
		int				dim = MakeDoGKernel( DoG, r1, r2, flog );

		Convolve( _avfflt, _avf, wf, hf,
			&DoG[0], dim, dim, true, true, kfft );
		Normalize( _avfflt );

		Convolve( _bvfflt, _bvf, wf, hf,
			&DoG[0], dim, dim, true, true, kfft );
		Normalize( _bvfflt );

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

			if( _avfflt.size() ) {
				Downsample( _avsflt, _avfflt );
				avs_aln = &_avsflt;
			}

			if( _bvfflt.size() ) {
				Downsample( _bvsflt, _bvfflt );
				bvs_aln = &_bvsflt;
			}
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

	ok = true;

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


