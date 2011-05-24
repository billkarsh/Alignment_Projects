

#include	"Inspect.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"






/* --------------------------------------------------------------- */
/* ABOverlay ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write out a diagnostic overlay image.
//
// A large black image is created to show A and B relative
// to each other. A is shown unmodified in GREEN. Regions
// of B that get mapped from A, or that fall outside A are
// shown in RED. Pixels in B that correspond to A but did
// not get a mapping are shown in BLUE.
//
void ABOverlay(
	const PixPair	&px,
	const uint16*	rmap,
	int				Ntrans,
	const TForm*	tfs,
	const TForm*	ifs )
{
	int		w		= px.wf,
			h		= px.hf,
			npix	= w * h;
	int		xmin	= -1000,			// overlay img limits...
			xmax	= w + 1000,
			ymin	= -1000,
			ymax	= h + 1000;
	int		w2		= xmax - xmin + 1,	// ...and its dims
			h2		= ymax - ymin + 1,
			npix2	= w2 * h2;
	int		bgaps	= 0;

// init overlay images

	vector<uint32>	raster2( npix2, 0xFF000000 );	// set alpha chan
	vector<float>	bpix( npix2, 0.0 );

// In this section...
// (a) Accumulate centroid data for each mapping region.
// (b) For each point in A that maps onto B:
//		(1) Store its tfs number in bmap.
//		(2) Store the targeted B pixel value in bpix.

	vector<uint16>	bmap( npix, 0 );
	vector<Point>	ctr( Ntrans, Point(0.0, 0.0) );
	vector<int>		ncpts( Ntrans, 0 );
	MeanStd			m;
	double			mean, std;

	for( int y = 0; y < h; ++y ) {

		for( int x = 0; x < w; ++x ) {

			int		mv = rmap[x + w*y] - 10;

			if( mv >= Ntrans ) {

				printf(
				"ABOverlay: ERROR - tform idx=%d, Ntrans=%d.\n",
				mv, Ntrans );

				return;
			}

			if( mv >= 0 ) {

				ctr[mv].x += x;
				ctr[mv].y += y;
				++ncpts[mv];

				Point p( x, y );
				tfs[mv].Transform( p );

				int	ix = (int)p.x;
				int	iy = (int)p.y;

				if( ix >= 0 && ix < w && iy >= 0 && iy < h )
					bmap[ix + w*iy] = mv + 10;

				// slightly stricter test for interpolation
				if( p.x >= 0.0 && p.x < w-1 &&
					p.y >= 0.0 && p.y < h-1 ) {

					double	pix =
					InterpolatePixel( p.x, p.y, *px.bvf_vfy, w );

					bpix[(x-xmin) + w2*(y-ymin)] = pix;
					m.Element( pix );
				}
			}
		}
	}

// Finish centroids and transform them to B-coords.

	for( int i = 0; i < Ntrans; ++i ) {

		ctr[i].x /= ncpts[i];
		ctr[i].y /= ncpts[i];

		printf(
			"ABOverlay: region %d, center (%f %f) in A,"
			" npixels=%d.\n",
			i, ctr[i].x, ctr[i].y, ncpts[i] );

		tfs[i].Transform( ctr[i] );

		printf( "ABOverlay: Maps to (%f %f) in image B.\n",
		ctr[i].x, ctr[i].y );
	}

// Now we scan the bmap and characterize/colorize the regions.
//
// (1) If the bmap has a defined entry, that pixel is already
//		included in the bpix RED channel, so skip it.
//
// (2) If the closest (by centroid) ifs maps back to A we add
//		it to the BLUE channel.
//
// (3) If outside A we add it to the RED channel to see both
//		A and B in the common global composite.

	for( int i = 0; i < npix; ++i ) {

		int		iy = i / w;
		int		ix = i - w * iy;

		// only work on undefined pixels

		if( bmap[ix + w*iy] < 10 ) {

			// if no transforms, can't even guess
			if( Ntrans <= 0 )
				continue;

			double	dbest	= 1.0E30;
			int		best	= -1;

			for( int q = 0; q < Ntrans; ++q ) {

				double	dx	= ix - ctr[q].x;
				double	dy	= iy - ctr[q].y;
				double	d	= dx*dx + dy*dy;

				if( d < dbest ) {
					dbest	= d;
					best	= q;
				}
			}

			// using best tform, inverse map back to A

			Point	p( ix, iy );
			ifs[best].Transform( p );

			if( p.x >= 0.0 && p.x < w && p.y >= 0.0 && p.y < h ) {

				// record gap pixels in BLUE
				ix = (int)p.x - xmin;
				iy = (int)p.y - ymin;
				raster2[ix + w2*iy] |= (255 << 16);
				++bgaps;
			}
			else {

				const vector<double>&	bv = *px.bvf_vfy;

				// Subtract (xmin, ymin) to place
				// in larger composite.

				p.x -= xmin;
				p.y -= ymin;

				if( p.x >= 0 && p.x < w2-1 &&
					p.y >= 0 && p.y < h2-1 ) {

					double	pix = bv[ix + w*iy];

					DistributePixel( p.x, p.y, pix, &bpix[0], w2 );
					m.Element( pix );
				}
			}
		}
	}

// Report bgaps = how many pixels in B are not in the image
// (range) of some tfs.

	printf(
	"ABOverlay: %d pixels in A-B overlap not from A.\n", bgaps );

// Transfer bpix to the RED channel.

	m.Stats( mean, std );

	if( std ) {

		for( int i = 0; i < npix2; ++i ) {

			if( bpix[i] != 0.0 ) {

				int	pix = 127 + int((bpix[i] - mean)/std*40.0);

				if( pix < 0 )
					pix = 0;
				else if( pix > 255 )
					pix = 255;

				raster2[i] |= pix;
			}
		}
	}

// Now make copy of A in the GREEN channel.

	const vector<double>&	av = *px.avf_vfy;

	for( int i = 0; i < npix; ++i ) {

		int	iy	= i / w;
		int	ix	= i - w * iy;
		int	pix	= 127 + int(40 * av[i]);

		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;

		ix -= xmin;
		iy -= ymin;

		if( ix >= 0 && ix < w2 && iy >= 0 && iy < h2 )
			raster2[ix + w2*iy] |= (pix << 8);
	}

// Write it out
//	Raster32ToTifRGBA( "comp.tif", &raster2[0], w2, h2 );
	Raster32ToPngRGBA( "comp.png", &raster2[0], w2, h2 );
}

/* --------------------------------------------------------------- */
/* Legal funs for CorrView --------------------------------------- */
/* --------------------------------------------------------------- */

// Is it a legal region?
//
static bool lr31x31( int sx, int sy, void *arg )
{
	return sx == 31 && sy == 31;
}


// Is the pixel count legal?
//
static bool lcTRUE( int c1, int c2, void *arg )
{
	return true;
}

/* --------------------------------------------------------------- */
/* CorrView ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Check results by visualizing correlation from
// aligning several randomly chosen small patches.
//
// Result saved as 'qual.tif'.
//
static void CorrView(
	const double*	apix,
	const uint8*	bpix,
	uint32			w,
	uint32			h,
	const uint16*	rmap )
{
	const int	PATCH	= 15;	// radius of patch
	const int	LOOK	= 47;	// how far to look (radius)

// Create qual, a map that tells where the quality is OK

	int		npix	= w * h;
	uint8*	qual	= (uint8*)RasterAlloc( npix * sizeof(uint8) );

	memset( qual, 0, npix );

// Do trial alignments of several small random patches

	for( int ntries = 0; ntries < 10000; ++ntries ) {

		// generate a random point in the range [LOOK, dim-1-LOOK]

		int x = LOOK + int((w-1-2*LOOK) * random() / RAND_MAX);
		int y = LOOK + int((h-1-2*LOOK) * random() / RAND_MAX);

		if( rmap[x + w*y] < 10 )
			continue;	// if not mapped, skip it

		// collect the data from the two images.

		int				nv = 0;	// # pixels with 'reasonable values'
		vector<Point>	pts1;
		vector<double>	vals1;
		MeanStd			m;

		for( int ix = x-PATCH; ix <= x+PATCH; ++ix ) {

			for( int iy = y-PATCH; iy <= y+PATCH; ++iy ) {

				uint8	v = 127 + int(40 * apix[ix + w*iy]);

				vals1.push_back( v );
				pts1.push_back( Point( ix-x, iy-y ) );
				nv += (v > 12);	// 12 is sort of arbitrary
				m.Element( v );
			}
		}

		int		side = 2*PATCH+1;
		double	frac = nv / double(side*side);

		if( frac < 0.1 ) {
			printf(
			"CorrView: Too many dark pixels"
			" (%d light of %d, %f) - skipped.\n",
			nv, side*side, frac );
			continue;
		}

		// now the same for the other image

		nv = 0;		// # pixels with 'reasonable values'
		vector<Point>	pts2;
		vector<double>	vals2;

		for( int ix = x-LOOK; ix <= x+LOOK; ++ix ) {

			for( int iy = y-LOOK; iy <= y+LOOK; ++iy ) {

				uint8	v = bpix[ix + w*iy];

				vals2.push_back( v );
				pts2.push_back( Point( ix-x, iy-y ) );
				nv += (v > 5);	// looking for folds here
			}
		}

		side = 2*LOOK+1;
		frac = nv / double(side*side);

		if( frac < 0.2 ) {
			printf( "CorrView: Too many dark pixels(2)"
			" (%d light of %d, %f) - skipped.\n",
			nv, side*side, frac );
			continue;
		}

		// now see how they align

		double		dx, dy;
		vector<CD>	ftc;	// cache for Fourier transform

		double co	= CorrPatches(
						stdout, false, dx, dy,
						pts1, vals1, pts2, vals2, 0, 0, 4000,
						lr31x31, NULL, lcTRUE, NULL, ftc );

		double	d	= sqrt( dx*dx + dy*dy );
		double	mean, std;
		m.Stats( mean, std );

		printf( "CorrView: pt %6d %6d ---> dx, dy= %6.2f, %6.2f,"
		" d=%6.2f, corr %6.3f, std %6.2f.\n\n",
		x, y, dx, dy, d, co, std );

		// where 'qual' is good, write a block into the 'qual' map

#if 0	// original 3-shade version with preset thresholds

		uint8	rslt = 0x7F;	// 127 = gray if we even tried it

		if( d < 10.0 && co > 0.5 )
			rslt = 0xC0;		// 192 = brighter gray for OK

		if( d < 6.0 && co > 0.6 )
			rslt = 0xFF;		// 255 = white if good match

#else	// continuous auto-shading

		int	rslt = 10;

		if( d < 10.0 && co > 0.1 )
			rslt = int(co * 100.0);

		if( d < 6.0 && co > 0.1 )
			rslt = int(100.0 + co * 100.0);

		if( rslt > 255 )
			rslt = 255;

#endif

		for( int ix = x-LOOK; ix <= x+LOOK; ++ix ) {

			for( int iy = y-LOOK; iy <= y+LOOK; ++iy ) {

				int	n = ix + w*iy;

				qual[n] = max( qual[n], rslt );
			}
		}
	}

	Raster8ToTif8( "qual.tif", qual, w, h );

	RasterFree( qual );
}

/* --------------------------------------------------------------- */
/* RunCorrView --------------------------------------------------- */
/* --------------------------------------------------------------- */

// (1) Build copy of the B-image and save that as 'registered.tif'.
//		If 'registered.tif' had already existed, the current B-image
//		is added into that and then saved.
//
// (2) The current A-image and new B-image are sent to CorrView.
//
void RunCorrView(
	const PixPair	&px,
	const uint16*	rmap,
	const TForm*	tfs,
	FILE*			flog )
{
	FILE*	f;
	uint8*	bpix;
	int		k, w, h, npix;

	w		= px.wf;
	h		= px.hf;
	npix	= w * h;

// Read existing file...

	if( f = fopen( "registered.tif", "r" ) ) {

		uint32	w2, h2;

		fclose( f );

		printf(
		"RunTripleChk: Reading existing 'registered.tif'.\n" );

		bpix = Raster8FromTif( "registered.tif", w2, h2, flog );
	}
	else {

		// ... or create new black raster

		printf(
		"RunTripleChk: No existing 'registered.tif'"
		" - creating one.\n" );

		bpix = (uint8*)RasterAlloc( npix * sizeof(uint8) );

		memset( bpix, 0, npix );
	}

// For each black pixel in bpix, if there is a mapping there
// from A to B, fill pixel with B-data.

	for( k = 0; k < npix; ++k ) {

		if( !bpix[k] ) {

			int		iy	= k / w;
			int		ix	= k - w * iy;
			int		mv	= rmap[ix + w*iy] - 10;

			if( mv >= 0 ) {	// a transformation exists

				Point	p( ix, iy );
				tfs[mv].Transform( p );

				if( p.x >= 0.0 && p.x < w-1 &&
					p.y >= 0.0 && p.y < h-1 ) {

					int	pix = 127 + int(40 *
					InterpolatePixel( p.x, p.y, *px.bvf_vfy, w ));

					if( pix < 0 )
						pix = 0;
					else if( pix > 255 )
						pix = 255;

					bpix[ix + w*iy] = pix;
				}
			}
		}
	}

	Raster8ToTif8( "registered.tif", bpix, w, h );

// Triple check, by picking regions at random and
// making sure they match.

	CorrView( &(*px.avf_vfy)[0], bpix, w, h, rmap );

	RasterFree( bpix );
}


