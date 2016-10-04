

#include	"Inspect.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"

#include	<stdlib.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* YellowView ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void YellowView(
    const PixPair	&px,
    const TAffine	&T,
    FILE			*flog )
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

    vector<uint32>			raster2( npix2, 0xFF000000 );
    const vector<double>&	av = *px.avf_vfy;

// A in green

    for( int i = 0; i < npix; ++i ) {

        int	ay	= i / w;
        int	ax	= i - w * ay;
        int	pix	= 127 + int(40 * av[i]);

        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;

        ax -= xmin;
        ay -= ymin;

        raster2[ax + w2*ay] |= (pix << 8);
    }

// B in red

    for( int i = 0; i < npix2; ++i ) {

        int		ry	= i / w2;
        int		rx	= i - w2 * ry;
        Point	p( rx + xmin, ry + ymin );

        T.Transform( p );

        if( p.x >= 0.0 && p.x < w-1 &&
            p.y >= 0.0 && p.y < h-1 ) {

            double	dpix =
            SafeInterp( p.x, p.y, &(*px.bvf_vfy)[0], w, h );

            int	pix	= 127 + int(40 * dpix);

            if( pix < 0 )
                pix = 0;
            else if( pix > 255 )
                pix = 255;

            raster2[i] |= pix;
        }
    }

    Raster32ToPngRGBA( "Ylw_Aff.png", &raster2[0], w2, h2, flog );
}

/* --------------------------------------------------------------- */
/* YellowView ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void YellowView(
    const PixPair	&px,
    const THmgphy	&T,
    FILE			*flog )
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

    vector<uint32>			raster2( npix2, 0xFF000000 );
    const vector<double>&	av = *px.avf_vfy;

// A in green

    for( int i = 0; i < npix; ++i ) {

        int	ay	= i / w;
        int	ax	= i - w * ay;
        int	pix	= 127 + int(40 * av[i]);

        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;

        ax -= xmin;
        ay -= ymin;

        raster2[ax + w2*ay] |= (pix << 8);
    }

// B in red

    for( int i = 0; i < npix2; ++i ) {

        int		ry	= i / w2;
        int		rx	= i - w2 * ry;
        Point	p( rx + xmin, ry + ymin );

        T.Transform( p );

        if( p.x >= 0.0 && p.x < w-1 &&
            p.y >= 0.0 && p.y < h-1 ) {

            double	dpix =
            SafeInterp( p.x, p.y, &(*px.bvf_vfy)[0], w, h );

            int	pix	= 127 + int(40 * dpix);

            if( pix < 0 )
                pix = 0;
            else if( pix > 255 )
                pix = 255;

            raster2[i] |= pix;
        }
    }

    Raster32ToPngRGBA( "Ylw_Hmg.png", &raster2[0], w2, h2, flog );
}

/* --------------------------------------------------------------- */
/* ABOverlay ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write out a diagnostic overlay png (comp_png).
//
// If comp_png omitted, default is './comp.png'.
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
    const TAffine*	tfs,
    const TAffine*	ifs,
    const char		*comp_png,
    FILE			*flog )
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

    if( !comp_png )
        comp_png = "comp.png";

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

                fprintf( flog,
                "ABOverlay: ERROR - tform idx=%d, Ntrans=%d\n",
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

        fprintf( flog,
            "ABOverlay: region %d, center (%f %f) in A,"
            " npixels=%d\n",
            i, ctr[i].x, ctr[i].y, ncpts[i] );

        tfs[i].Transform( ctr[i] );

        fprintf( flog, "ABOverlay: Maps to (%f %f) in image B.\n",
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

                double	pix = bv[ix + w*iy];

                DistributePixel( p.x, p.y, pix, &bpix[0], w2, h2 );
                m.Element( pix );
            }
        }
    }

// Report bgaps = how many pixels in B are not in the image
// (range) of some tfs.

    fprintf( flog,
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
    Raster32ToPngRGBA( comp_png, &raster2[0], w2, h2, flog );
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
    const uint16*	rmap,
    FILE			*flog )
{
    const int	PATCH	= 15;	// radius of patch
    const int	LOOK	= 47;	// how far to look (radius)

// Create qual, a map that tells where the quality is OK

    vector<uint8> qual( w * h, 0 );

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
            fprintf( flog,
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
            fprintf( flog,
            "CorrView: Too many dark pixels(2)"
            " (%d light of %d, %f) - skipped.\n",
            nv, side*side, frac );
            continue;
        }

        // now see how they align

        double		dx, dy;
        vector<CD>	ftc;	// cache for Fourier transform

        double co	= CorrPatches(
                        flog, false, dx, dy,
                        pts1, vals1, pts2, vals2, 0, 0, 4000,
                        lr31x31, NULL, lcTRUE, NULL, ftc );

        double	d	= sqrt( dx*dx + dy*dy );
        double	mean, std;
        m.Stats( mean, std );

        fprintf( flog,
        "CorrView: pt %6d %6d ---> dx, dy= %6.2f, %6.2f,"
        " d=%6.2f, corr %6.3f, std %6.2f\n\n",
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

    Raster8ToTif8( "qual.tif", &qual[0], w, h, flog );
}

/* --------------------------------------------------------------- */
/* RunCorrView --------------------------------------------------- */
/* --------------------------------------------------------------- */

// (1) Build copy of the B-image and save that as png registered_png.
//		If registered_png had already existed, the current B-image
//		is added into that and then saved.
//
// (2) If registered_png omitted, default is './registered.png'.
//
// (3) If heatmap is true, the current A-image and new B-image
//		are sent to CorrView.
//
void RunCorrView(
    const PixPair	&px,
    const uint16*	rmap,
    const TAffine*	tfs,
    bool			heatmap,
    const char		*registered_png,
    FILE			*flog )
{
    FILE*	f;
    uint8*	bpix;
    int		k, w, h, npix;

    if( !registered_png )
        registered_png = "registered.png";

    fprintf( flog, "\n" );

    w		= px.wf;
    h		= px.hf;
    npix	= w * h;

// Read existing file...

    if( f = fopen( registered_png, "r" ) ) {

        uint32	w2, h2;

        fclose( f );

        fprintf( flog,
        "RunTripleChk: Reading existing registered_png [%s].\n",
        registered_png );

        bpix = Raster8FromPng( registered_png, w2, h2, flog );
    }
    else {

        // ... or create new black raster

        fprintf( flog,
        "RunTripleChk: No previous registered file"
        " - creating one.\n" );

        bpix = (uint8*)RasterAlloc( npix * sizeof(uint8) );

        memset( bpix, 0, npix );
    }

// For each black pixel in bpix, if there is a mapping there
// from A to B, fill pixel with B-data. Note that we defer
// writing into bpix until we can normalize the new data.

    vector<double>	new_val;
    vector<int>		new_idx;

    for( k = 0; k < npix; ++k ) {

        if( !bpix[k] ) {

            int		iy	= k / w;
            int		ix	= k - w * iy;
            int		mv	= rmap[ix + w*iy] - 10;

            if( mv >= 0 ) {	// a transformation exists

                Point	p( ix, iy );
                tfs[mv].Transform( p );

#if 1	// using bicubic
                if( p.x >= 1.0 && p.x < w-2 &&
                    p.y >= 1.0 && p.y < h-2 ) {

                    new_val.push_back(
                    BiCubicInterp( &(*px.bvf_vfy)[0], w, p ) );

                    new_idx.push_back( ix + w*iy );
                }
#else	// bilinear
                if( p.x >= 0.0 && p.x < w-1 &&
                    p.y >= 0.0 && p.y < h-1 ) {

                    new_val.push_back(
                    InterpolatePixel( p.x, p.y, *px.bvf_vfy, w ) );

                    new_idx.push_back( ix + w*iy );
                }
#endif
            }
        }
    }

// Normalize and store new values

    int	np = new_val.size();

    if( np ) {

        Normalize( new_val );

        for( int k = 0; k < np; ++k ) {

            double	pix	= 127 + 40.0 * new_val[k];

            if( pix < 0 )
                pix = 0;
            else if( pix > 255 )
                pix = 255;

            bpix[new_idx[k]] = int(pix + 0.5);
        }
    }

    Raster8ToPng8( registered_png, bpix, w, h, flog );

// Triple check, by picking regions at random and
// making sure they match.

    if( heatmap )
        CorrView( &(*px.avf_vfy)[0], bpix, w, h, rmap, flog );

    RasterFree( bpix );

    fprintf( flog, "\n" );
}


