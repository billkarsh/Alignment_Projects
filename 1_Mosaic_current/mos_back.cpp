
// Make a mosaic


#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"LinEqu.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"Draw.h"
#include	"Memory.h"

#include	<string.h>

#include	<sys/resource.h>
#include	<sys/stat.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	kCorrRadius	200

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class glob_spot {
public:
	vector<int>		which;	// which picture
	vector<int>		patch;	// which patch
	vector<Point>	where;	// local coords
	Point			gpt;	// point in global coords spawning these
	double			moved;	// distance it moved
};


class image {
public:
	uint8*					raster;		// gray scale image
	uint8*					foldmap;	// fold map
	uint8*					bmap;		// boundary map
	uint32					w, h;
	int						layer;		// layer number
	int						FirstGlobalPoint;
	int						FirstTriangle;
	int						spbase;		// should be added to all sp ids in this image, for uniqueness
	char					*rname;		// raster filename
	char					*fname;		// foldmap filename
	char					*spname;	// super-pixel filename
	char					*bname;		// boundary-map filename
	uint16					*spmap;
	vector<int>				SPmapping;  // tells what original SP numbers are mapped to
	vector<TAffine>			tf;			// image to global space, one for each patch (0 unused)
	vector<TAffine>			inv;		// inverse transform
	vector<vector<TAffine> >sectors;	// forward transform for each sector of each patch
	vector<vector<TAffine> >sinvs;		// sector inverses
	vector<uint8>			FoldmapRenumber;
};

// Class describing a triple of numbers that in turn describe
// where the point is that maps to a given point in global space.
// It comes from a particular image, a particular patch within
// the image, and a sector within the patch (the sectors are
// used to deform the patch so the edges line up better.)
//
class Triple {
public:
	uint16 image;
	uint8  patch;
	uint8  sector;
public:
	Triple()
	: image(0), patch(0), sector(0) {};

	Triple( uint16 i, uint8 p, uint8 s )
	: image(i), patch(p), sector(s) {};
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static FILE *of;
static bool dp = false;				// debug print
static int nwarn_edge_interp = 0;	// interpolation cannot be done since max is on edge of region
static int nwarn_bad_corr = 0;		// could not find a good correlation
static int nwarn_good_on_edge = 0;	// worst case - good correlation, but on edge






/* --------------------------------------------------------------- */
/* FileNameToLayerNumber ----------------------------------------- */
/* --------------------------------------------------------------- */

// Find the layer number. We have a list of sub-strings; if the
// file name contains the sub-string then that defines the layer.
// Sub-strings and layer numbers stored in two parallel vectors.
//
// BK added: If fname is an nmrc file, then the name has the form
// 'nmrc_z_tile.png' and we decode the z from that.
//
static int FileNameToLayerNumber(
	const vector<char*>	&dnames,
	const vector<int>	&lnums,
	const char			*fname )
{
/* ----- */
/* nmrc? */
/* ----- */

	const char	*s = FileNamePtr( fname );

	if( !strncmp( s, "nmrc_", 5 ) )
		return atoi( s + 5 );

/* ------------------ */
/* Otherwise use ldir */
/* ------------------ */

int k;
for(k=0; k<dnames.size(); k++)
    if( strstr(fname, dnames[k]) != NULL )
	break;
if( k >= dnames.size() ) {
    printf("Oops - '%s' not in layer directory.\n", fname);
    exit( 42 );
    }
return lnums[k];
}

/* --------------------------------------------------------------- */
/* PrintTransAndInv ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintTransAndInv( FILE *of, const TAffine &tr )
{
	TAffine	in;

	in.InverseOf( tr );

	fprintf( of,
	"Fwd:%11.8f %11.8f %10.4f   %11.8f %11.8f %10.4f\n",
	tr.t[0], tr.t[1], tr.t[2],
	tr.t[3], tr.t[4], tr.t[5] );

	fprintf( of,
	"Rev:%11.8f %11.8f %10.4f   %11.8f %11.8f %10.4f\n",
	in.t[0], in.t[1], in.t[2],
	in.t[3], in.t[4], in.t[5] );
}

/* --------------------------------------------------------------- */
/* GlobalBounds -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Find bounding box in global space.
//
static void GlobalBounds(
	double				&xmin,
	double				&xmax,
	double				&ymin,
	double				&ymax,
	const vector<image>	&images,
	int					w,
	int					h )
{
	if( xmin > BIG/2 ) {

		int	ni = images.size();

		for( int i = 0; i < ni; ++i ) {

			int	nt = images[i].tf.size();

			for( int k = 1; k < nt; ++k ) {

				const TAffine&	T = images[i].tf[k];

				if( !T.det() )
					continue;

				vector<Point>	pts;

				pts.push_back( Point( 0.0, 0.0 ) );
				pts.push_back( Point( w-1, 0.0 ) );
				pts.push_back( Point( w-1, h-1 ) );
				pts.push_back( Point( 0.0, h-1 ) );

				for( int j = 0; j < 4; ++j ) {

					T.Transform( pts[j] );

					xmin = fmin( xmin, pts[j].x );
					xmax = fmax( xmax, pts[j].x );
					ymin = fmin( ymin, pts[j].y );
					ymax = fmax( ymax, pts[j].y );
				}
			}
		}
	}

	printf( "Bounds of global image: x=[%f %f] y=[%f %f].\n",
	xmin, xmax, ymin, ymax );
}

/* --------------------------------------------------------------- */
/* OffsetCoords -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Modify bounds and transformations so the min is (0,0)
// then apply user's x0, y0.
//
// We change by an integer amount to make layer-to-layer
// alignment easier.
//
static void OffsetCoords(
	double			&xmin,
	double			&xmax,
	double			&ymin,
	double			&ymax,
	vector<image>	&images,
	vector<double>	&x1,
	vector<double>	&x2,
	vector<double>	&y1,
	vector<double>	&y2,
	int				x0,
	int				y0,
	int				xsize,
	int				ysize )
{
	int xfl = int(floor( xmin ));
	int yfl = int(floor( ymin ));

// bounds

	xmin -= xfl;
	xmax -= xfl;
	ymin -= yfl;
	ymax -= yfl;

	printf( "Bounds of global image: x=[%f %f] y=[%f %f].\n",
	xmin, xmax, ymin, ymax );

	if( xsize != -1 ) {
		xmax = xsize;
		ymax = ysize;
	}

// transforms

	xfl = -(xfl + x0);
	yfl = -(yfl + y0);

	int	ni = images.size();

	for( int i = 0; i < ni; ++i ) {

		int	nt = images[i].tf.size();

		for( int k = 1; k < nt; ++k ) {

			TAffine&	T = images[i].tf[k];

			if( !T.det() )
				continue;

			T.AddXY( xfl, yfl );
			images[i].inv[k].InverseOf( T );
		}
	}

// debugging - print transforms

#if 0
	for( int i = 0; i < ni; ++i ) {

		int	nt = images[i].tf.size();

		for( int k = 1; k < nt; ++k ) {

			const TAffine&	T = images[i].tf[k];

			if( !T.det() )
				continue;

			printf(
			"Post offset: image %d, layer %d tf[%d]"
			" = %f %f %f   %f %f %f\n",
			i, images[i].layer, k,
			T.t[0], T.t[1], T.t[2],
			T.t[3], T.t[4], T.t[5] );
		}
	}
#endif

// correspondence points

	int	n = x1.size();

	for( int i = 0; i < n; ++i ) {

		x1[i] += xfl;
		x2[i] += xfl;
		y1[i] += yfl;
		y2[i] += yfl;
	}
}

/* --------------------------------------------------------------- */
/* LoadNormImg --------------------------------------------------- */
/* --------------------------------------------------------------- */

static uint8* LoadNormImg( const char *name, uint32 &w, uint32 &h )
{
	uint8*	ras = Raster8FromAny( name, w, h );

	if( !strstr( name, ".mrc" ) ) {

		double	mean, std;
		int		np  = w * h;
		MeanStd	m;

		for( int i = 0; i < np; ++i )
			m.Element( double(ras[i]) );

		m.Stats( mean, std );

		for( int i = 0; i < np; ++i ) {

			int pix = RND( 127 + (ras[i] - mean) / std * 30.0 );

			if( pix > 255 )
				pix = 255;

			ras[i] = pix;
		}
	}

	return ras;
}

/* --------------------------------------------------------------- */
/* PrintMagnitude ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void PrintMagnitude( const vector<double> &X )
{
	int	k = X.size() - 6;

	if( k >= 0 ) {

		double	mag	= sqrt( X[k]*X[k] + X[k+1]*X[k+1] );

		printf( "Final magnitude is %g\n", mag );
	}
}

/* --------------------------------------------------------------- */
/* lr ------------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Legal region for correlation?
//
static bool lr( int sx, int sy, void *arg )
{
	return sx == 21 && sy == 21;
}

/* --------------------------------------------------------------- */
/* lc ------------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Legal count for correlation?
//
static bool lc( int c1, int c2, void *arg )
{
	return true;
}

/* --------------------------------------------------------------- */
/* CommonPoint --------------------------------------------------- */
/* --------------------------------------------------------------- */

// In the set of 'images', point pt (in global space) occurs
// in more than one image. Find out which, create a spot, and
// add it to the vector of spots. We know for sure the point
// is part of the image i1, patch patch1.
//
static void CommonPoint(
	vector<glob_spot>		&spots,
	const vector<image>		&images,
	Point					pt,
	int						i1,
	int						patch1,
	int						out_layer )
{
Point p1 = pt;
printf("image %d, global x,y= %f %f\n", i1, p1.x, p1.y);
images[i1].inv[patch1].Transform( p1 );
const int PATCH=10;  // radius of patch
//const int LOOK=25;   // radius of region to look.  For Marta, must be smaller
const int LOOK=52;   // radius of region to look
// first, make sure we have enough pixels in the image
int p1x = RND(p1.x);
int p1y = RND(p1.y);
if( p1x < PATCH || p1x > images[i1].w-1-PATCH ) return;
if( p1y < PATCH || p1y > images[i1].h-1-PATCH ) return;

// Surrounding pixels exist, so go find the data
vector<Point> pts1;
vector<double> vals1;
int nv = 0;   // number of pixels with 'reasonable values'
for(int ix = p1x-PATCH; ix <= p1x+PATCH; ix++) {
    for(int iy = p1y-PATCH; iy <= p1y+PATCH; iy++) {
        uint8 v = images[i1].raster[ix + images[i1].w * iy];
        vals1.push_back(v);
        pts1.push_back(Point(ix-p1x, iy-p1y));
        nv += (v > 12);  // 12 is sort of arbitrary
        }
    }
int side = 2*PATCH+1;
double frac = nv/double(side*side);
if( frac < 0.5 ) {
    printf("Too many dark pixels (%d light of %d, %f) - skipped\n", nv, side*side, frac);
    return;
    }

// so now we have the master copy.   Find all images in which it occurs, including the
// one we started with (clearly, this one should end up with a correlation very near 1.0,
// and an offset of very close to (0,0). )
glob_spot s;   // this will contain all the detailed info
s.gpt = pt;
for(int i=0; i<images.size(); i++) {
    if( images[i].layer != out_layer )
	continue;
    // need to go through all patches.  Just because the inverse transform maps into the image,
    // it's not enough.  It also needs to land in the corresponding patch
    for(int patch=1; patch < images[i].tf.size(); patch++) {
        if( images[i].tf[patch].det() == 0.0 )
	    continue;
        Point p1 = pt;
        images[i].inv[patch].Transform( p1 );
        //printf(" image %d, x,y= %f %f\n", i, p1.x, p1.y);
        // first, make sure it lands in the image, and we have enough pixels in the image
        int p1x = RND(p1.x);
        int p1y = RND(p1.y);
        if( p1x < LOOK || p1x > images[i].w-1-LOOK ) continue;
        if( p1y < LOOK || p1y > images[i].h-1-LOOK ) continue;

        bool nip = false;   // stands for Not In Patch
        vector<Point> pts2;
        vector<double> vals2;
        for(int ix = p1x-LOOK; ix <= p1x+LOOK; ix++) {
            for(int iy = p1y-LOOK; iy <= p1y+LOOK; iy++) {
                int loc = ix + images[i].w * iy;
                nip = nip || (images[i].foldmap[loc] != patch);
                uint8 v = images[i].raster[loc];
                vals2.push_back(v);
                pts2.push_back(Point(ix-p1x, iy-p1y));
                }
            }
        // if any of the pixels were not in the relevant patch, we should not use this
        if( nip )
	    continue;

        double dx, dy;
        vector<CD> ftc;  // cache for Fourier transform
        int save = nwarn_edge_interp;
        // dp = (pt.x == 5570.0 && pt.y == 10356);
        // dp = false;

		double	co = CorrPatches(
						stdout, false, dx, dy,
						pts1, vals1, pts2, vals2,
						0, 0, kCorrRadius,
						lr, NULL, lc, NULL, ftc );

        printf(" tx=%f, ---> dx, dy= %f, %f, corr %f\n\n", pt.x, dx, dy, co);
        if( co < 0.7 ) {
             printf("WARNING - poor image correlation %f - ignored\n", co);
             nwarn_bad_corr++;
             nwarn_edge_interp = save;  // reset this since we don't care since we are not using the point
             continue;
             }
        // OK, got a matching point in this picture.  Add to lists.
        p1.x += dx;
        p1.y += dy;
        s.which.push_back(i);      // image i
        s.patch.push_back(patch);  // patch 'patch'
        s.where.push_back(p1);     // and location in local coordinates
        }
    }
// OK, now s contains all the matches. If there's more than one, we've got a stitch point
if( s.which.size() > 1 ) {
    printf("More than one image for global point (%f %f)\n", pt.x, pt.y);
    // and add all the matching spots to the list for each image
    for(int i=0; i<s.which.size(); i++) {
	printf("  Image %d, patch %d, coords (%f %f)\n", s.which[i], s.patch[i], s.where[i].x, s.where[i].y);
	}
    spots.push_back(s);  // so we can re-evaluate residuals later.
    }
}

/* --------------------------------------------------------------- */
/* DrawText ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Super simple routine for drawing text. Uses simple 5x7 font.
// Font by Imran Nazar under BSD Public License.
//
static uint8 font5x7[][7]={
/*   */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00},
/* ! */ {0x04,0x04,0x04,0x04,0x00,0x00,0x04},
/* " */ {0x0A,0x0A,0x0A,0x00,0x00,0x00,0x00},
/* # */ {0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A},
/* $ */ {0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04},
/* % */ {0x18,0x19,0x02,0x04,0x08,0x13,0x03},
/* & */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00},
/* ' */ {0x0C,0x04,0x08,0x00,0x00,0x00,0x00},
/* ( */ {0x02,0x04,0x08,0x08,0x08,0x04,0x02},
/* ) */ {0x08,0x04,0x02,0x02,0x02,0x04,0x08},
/* * */ {0x00,0x04,0x15,0x0E,0x15,0x04,0x00},
/* + */ {0x00,0x04,0x04,0x1F,0x04,0x04,0x00},
/* , */ {0x00,0x00,0x00,0x00,0x0C,0x04,0x08},
/* - */ {0x00,0x00,0x00,0x1F,0x00,0x00,0x00},
/* . */ {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C},
/* / */ {0x00,0x01,0x02,0x04,0x08,0x10,0x00},
/* 0 */ {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E},
/* 1 */ {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E},
/* 2 */ {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F},
/* 3 */ {0x1F,0x02,0x04,0x02,0x01,0x11,0x0E},
/* 4 */ {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02},
/* 5 */ {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E},
/* 6 */ {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E},
/* 7 */ {0x1F,0x01,0x02,0x04,0x08,0x08,0x08},
/* 8 */ {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E},
/* 9 */ {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C},
/* : */ {0x00,0x0C,0x0C,0x00,0x0C,0x0C,0x00},
/* ; */ {0x00,0x0C,0x0C,0x00,0x0C,0x04,0x08},
/* < */ {0x02,0x04,0x08,0x10,0x08,0x04,0x02},
/* = */ {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00},
/* > */ {0x08,0x04,0x02,0x01,0x02,0x04,0x08},
/* ? */ {0x0E,0x11,0x01,0x02,0x04,0x00,0x04},
/* @ */ {0x0E,0x11,0x01,0x0D,0x15,0x15,0x0E},
/* A */ {0x0E,0x11,0x11,0x11,0x1F,0x11,0x11},
/* B */ {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E},
/* C */ {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E},
/* D */ {0x1C,0x12,0x11,0x11,0x11,0x12,0x1C},
/* E */ {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F},
/* F */ {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10},
/* G */ {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F},
/* H */ {0x11,0x11,0x11,0x1F,0x11,0x11,0x11},
/* I */ {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E},
/* J */ {0x07,0x02,0x02,0x02,0x02,0x12,0x0C},
/* K */ {0x11,0x12,0x14,0x18,0x14,0x12,0x11},
/* L */ {0x10,0x10,0x10,0x10,0x10,0x10,0x1F},
/* M */ {0x11,0x1B,0x15,0x15,0x11,0x11,0x11},
/* N */ {0x11,0x11,0x19,0x15,0x13,0x11,0x11},
/* O */ {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E},
/* P */ {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10},
/* Q */ {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D},
/* R */ {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11},
/* S */ {0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E},
/* T */ {0x1F,0x04,0x04,0x04,0x04,0x04,0x04},
/* U */ {0x11,0x11,0x11,0x11,0x11,0x11,0x0E},
/* V */ {0x11,0x11,0x11,0x11,0x11,0x0A,0x04},
/* W */ {0x11,0x11,0x11,0x15,0x15,0x15,0x0A},
/* X */ {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11},
/* Y */ {0x11,0x11,0x11,0x0A,0x04,0x04,0x04},
/* Z */ {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F},
/* [ */ {0x0E,0x08,0x08,0x08,0x08,0x08,0x0E},
/* \ */ {0x00,0x10,0x08,0x04,0x02,0x01,0x00},
/* ] */ {0x0E,0x02,0x02,0x02,0x02,0x02,0x0E},
/* ^ */ {0x04,0x0A,0x11,0x00,0x00,0x00,0x00},
/* _ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x1F},
/* ` */ {0x08,0x04,0x02,0x00,0x00,0x00,0x00},
/* a */ {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F},
/* b */ {0x10,0x10,0x10,0x16,0x19,0x11,0x1E},
/* c */ {0x00,0x00,0x0E,0x10,0x10,0x11,0x0E},
/* d */ {0x01,0x01,0x01,0x0D,0x13,0x11,0x0F},
/* e */ {0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E},
/* f */ {0x06,0x09,0x08,0x1C,0x08,0x08,0x08},
/* g */ {0x00,0x0F,0x11,0x11,0x0F,0x01,0x0E},
/* h */ {0x10,0x10,0x16,0x19,0x11,0x11,0x11},
/* i */ {0x00,0x04,0x00,0x04,0x04,0x04,0x04},
/* j */ {0x02,0x00,0x06,0x02,0x02,0x12,0x0C},
/* k */ {0x10,0x10,0x12,0x14,0x18,0x14,0x12},
/* l */ {0x04,0x04,0x04,0x04,0x04,0x04,0x0F},
/* m */ {0x00,0x00,0x1A,0x15,0x15,0x11,0x11},
/* n */ {0x00,0x00,0x16,0x19,0x11,0x11,0x11},
/* o */ {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E},
/* p */ {0x00,0x00,0x1E,0x11,0x1E,0x10,0x10},
/* q */ {0x00,0x00,0x0D,0x13,0x0F,0x01,0x01},
/* r */ {0x00,0x00,0x16,0x19,0x10,0x10,0x10},
/* s */ {0x00,0x00,0x0E,0x10,0x0E,0x01,0x1E},
/* t */ {0x08,0x08,0x1C,0x08,0x08,0x09,0x06},
/* u */ {0x00,0x00,0x11,0x11,0x11,0x13,0x0D},
/* v */ {0x00,0x00,0x11,0x11,0x11,0x0A,0x04},
/* w */ {0x00,0x00,0x11,0x11,0x15,0x15,0x0A},
/* x */ {0x00,0x00,0x11,0x0A,0x04,0x0A,0x11},
/* y */ {0x00,0x00,0x11,0x11,0x0F,0x01,0x0E},
/* z */ {0x00,0x00,0x1F,0x02,0x04,0x08,0x1F},
/* { */ {0x02,0x04,0x04,0x08,0x04,0x04,0x02},
/* | */ {0x04,0x04,0x04,0x04,0x04,0x04,0x04},
/* } */ {0x08,0x04,0x04,0x02,0x04,0x04,0x08},
/* ~ */ {0x00,0x00,0x04,0x15,0x02,0x00,0x00}};


static void DrawText(
	uint8		*raster,
	int			w,
	int			h,
	const char	*str,
	int			x,
	int			y )
{
	for( const char *p = str; *p; p++ ) {	// go through the string

		int	ch = *p;

		if( ch < 32 || ch > 126 )
			continue;

		// draw as N times original size
		//printf("draw character '%c' at %d, %d\n", ch, x, y);

		ch -= 32;

		int	N = 3;

		for( int dy = 0; dy < 7; ++dy ) {

			for( int dx = 0; dx < 5; ++dx ) {

				bool	bit = (font5x7[ch][dy] & (1 << (4-dx)));

				if( !bit )
					continue;

				int		bx = x + dx*N;
				int		by = y + dy*N;

				// draw an N x N box at (bx, by)

				for( int j = 0; j < N; ++j ) {

					for( int k = 0; k < N; ++k ) {

						int	m = bx+j + w*(by+k);

						if( m < w*h )
							raster[m] = 0xFF;
					}
				}
			}
		}

		x += 6*N;
	}
}

/* --------------------------------------------------------------- */
/* PseudoAngle --------------------------------------------------- */
/* --------------------------------------------------------------- */

#if 0

// An angle like atan2(y,x), but runs from -4 to +4,
// and is continuous, but not uniform.
//
static double PseudoAngle( double y, double x )
{
bool xish = fabs(y) <= fabs(x);  // closer to x than y axis
if( x > 0 && xish )
	return y/x;
if( y > 0 ) {
    if( xish )
	return 4 + y/x;
    else
        return 2 - x/y;
    }
// now we know y < 0
if( xish )
    return -4 +y/x;
return -2 - x/y;
}

#endif

/* --------------------------------------------------------------- */
/* FillInHolesInFoldmap ------------------------------------------ */
/* --------------------------------------------------------------- */

// We consider every 0 in the foldmap, and decide whether to
// fill it in. We will change each 0 to 255 as we examine it
// so we only see it once. After we see each hole, we will
// either fill it in, or add it to the list of 0s to be
// restored at the end. Basic idea is to expand out through
// each 0 area, and see what we hit. If we hit only one thing
// then it's a hole. If one thing and the boundary, also fill
// it in. Two different things, then it's a gap between two
// patches; leave it a hole.
//
static void FillInHolesInFoldmap( uint8 *map, uint32 w, uint32 h )
{
uint32 np = w*h;  // number of pixels
int start = 0;
vector<int> setback;  // pixels we will set back to 0 when we are done
int nholes = 0;
for(int i=start; i < np; i++) {
    if( map[i] == 0 ) {
	start = i;  // so next time resume here
	stack<int> st;
	st.push(i);
	vector<int> boundary( 257, 0 );	// histogram of what we hit.
        vector<int> pixels; // 0 valued pixels in this area
	while( !st.empty() ) {
	    int j = st.top(); st.pop();
            if( map[j] == 0 ) {
	        map[j] = 255;  // so we won't get it again
                pixels.push_back(j);
	        int y = j / w;
	        int x = j - w * y;
	        if( x-1 >= 0 ) st.push(j-1); else boundary[256]++;
            if( x+1 <  w ) st.push(j+1); else boundary[256]++;
            if( y-1 >= 0 ) st.push(j-w); else boundary[256]++;
            if( y+1 <  h ) st.push(j+w); else boundary[256]++;
		}
            else {// map value is not 0, add what we ran into to the set of boundary values
		if( map[j] != 255 ) boundary[map[j]]++;
		}
	    }
        // how many different things did we hit?  And what was one of them?
        int howmany = 0, what = 0;
        for( int k = 0; k < 256; ++k ) {
	    if( boundary[k] != 0 ) {
			howmany++;
            what = k;
		}
	    }
        if( howmany == 1 ) {	// then we have a hole, or a gap starting at a boundary
	    for( int k = 0; k < pixels.size(); ++k )
			map[pixels[k]] = what;
        nholes++;
	    }
	else { // not a pure hole.  record pixel values to set back later
            for(int k=0; k<pixels.size(); k++)
		setback.push_back(pixels[k]);
	     }
	}
    }
// OK, reset all those we decided not to fill in
for(int k=0; k<setback.size(); k++)
    map[setback[k]] = 0;
printf("Filled in %d holes\n", nholes);
}

/* --------------------------------------------------------------- */
/* cot ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

static double cot( double x )
{
	return 1.0 / tan( x );
}

/* --------------------------------------------------------------- */
/* Inside -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Inside( triangle &tri, vector<Point> gvtx, Point &p )
{
Point p1(gvtx[tri.v[0]]);
Point p2(gvtx[tri.v[1]]);
Point p3(gvtx[tri.v[2]]);
return LeftSide(p1,p2,p) && LeftSide(p2,p3,p) && LeftSide(p3,p1,p);
}

/* --------------------------------------------------------------- */
/* FindTriangle -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Find the triangle where a point resides. This calculation
// can be done equally well in global or local coordinates.
// We will use global coords. Last two args are returned.
// Which 3 variables, and which 3 proportions.
//
static void FindTriangle(
	vector<image>		&images,
	int					im,
	int					patch,
	Point				p1,
	vector<triangle>	&tris,
	vector<Point>		&gvtx,
	int					N,
	int					*which,
	double				*amts )
{
Point p = p1;
images[im].tf[patch].Transform( p );
int first_tri = images[im].FirstTriangle + (patch-1)*N;  // since N triangle per patch, and
					                // first one is #1
for(int i=first_tri; i < first_tri+N; i++) {
    if( Inside(tris[i], gvtx, p) ) {
	printf("Got one - Point %f %f inside triangle (%f %f) (%f %f) (%f %f)\n",
	 p.x, p.y,
	 gvtx[tris[i].v[0]].x, gvtx[tris[i].v[0]].y,
	 gvtx[tris[i].v[1]].x, gvtx[tris[i].v[1]].y,
	 gvtx[tris[i].v[2]].x, gvtx[tris[i].v[2]].y );
        // find barycentric coordinates (linear combo of the three vertices
        amts[0] = tris[i].a[0][0]*p.x + tris[i].a[0][1]*p.y + tris[i].a[0][2]*1.0;
        amts[1] = tris[i].a[1][0]*p.x + tris[i].a[1][1]*p.y + tris[i].a[1][2]*1.0;
        amts[2] = tris[i].a[2][0]*p.x + tris[i].a[2][1]*p.y + tris[i].a[2][2]*1.0;
        printf("Barycentric coords: %f %f %f\n", amts[0], amts[1], amts[2] );
        for(int k=0; k<3; k++)
	    which[k] = tris[i].v[k];
	return;
	}
    }
printf("Oops - no triangle for point (%f %f)\n", p.x, p.y);
exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteTriangles ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Write triangles out to a text file in a format gnuplot can read.
//
static void WriteTriangles(
	const char			*name,
	vector<triangle>	&tris,
	vector<Point>		&pts )
{
FILE *ft = FileOpenOrDie( name, "w" );

for(int i=0; i<tris.size(); i++) {
    Point ps[3];
    for(int k=0; k<3; k++)
	ps[k] = pts[tris[i].v[k]];
    for(int k=0; k<4; k++)
	fprintf(ft, "%f %f\n", ps[k%3].x, ps[k%3].y);
    fprintf(ft, "\n");
    }
fclose(ft);
}

/* --------------------------------------------------------------- */
/* FoldmapRenumber ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void FoldmapRenumber( image &im )
{
// if not computed yet, make the renumbering
if( im.FoldmapRenumber.size() == 0 ) {
    // compute it
    im.FoldmapRenumber.push_back(0);  // 0 always maps to 0
    int k=1;
    for(int i=1; i<im.tf.size(); i++) {
	if( im.tf[i].det() != 0.0 ) {           // used
	    im.tf[k]  = im.tf[i];
            im.inv[k] = im.inv[i];
            im.FoldmapRenumber.push_back(k++);
            }
        else                                   // not used - map to 0
            im.FoldmapRenumber.push_back(0);
        }
    im.tf.erase(im.tf.begin()+k,   im.tf.end());   // remove any unneeded transforms
    im.inv.erase(im.inv.begin()+k, im.inv.end());  // and inverses
    }
uint32 np = im.w * im.h;
uint8 m = im.FoldmapRenumber.size();
printf("Renumber %d pixels:", np);
for(int i=0; i<im.FoldmapRenumber.size(); i++)
    printf(" %d", im.FoldmapRenumber[i]);
printf("\n");
for(int i=0; i<np; i++) {
    uint8 pix = im.foldmap[i];
    im.foldmap[i] = (pix >= m) ? 0 : im.FoldmapRenumber[pix];
    }
}

/* --------------------------------------------------------------- */
/* RemapSuperPixelsOneImage -------------------------------------- */
/* --------------------------------------------------------------- */

// Modify a super-pixel map so numbers are assigned consectively
// from MaxSPUsed+1 on up, and update MaxSPUsed accordingly.
// Fill the vector SPmapping to record what was done.
//
// So if input was 0 0 0 0 0 0 1 2 4 0 0 0 0 and MAXSPUsed = 3000
// The output is   0 0 0 0 0 0 1 2 3 0 0 0 0
// and SPMapping[0] = 0
//     SPMapping[1] = 3001
//     SPMapping[2] = 3002
//     SPMapping[4] = 3003
// on exit, MAXSPused = 3003
//
static void RemapSuperPixelsOneImage(
	uint16		*test,
	int			w,
	int			h,
	int			&MaxSPUsed,
	vector<int>	&SPmapping )
{
int biggest = -1;
vector<int> vals(65536,0);
for(int i=0; i<w*h; i++) {
    vals[test[i]]++;
    biggest = max(int(test[i]), biggest);
    }
printf("Biggest value in 16 bit map is %d\n", biggest);
SPmapping.resize(biggest+1,0);  // make the remapping vector the right size
int base = MaxSPUsed;  // will add this to all
int n=0;
for(int i=1; i<65536; i++) {  // 0 always maps to 0
    if( vals[i] != 0 ) { // this value was used
	MaxSPUsed++;
        SPmapping[i] = MaxSPUsed;
	n++;
	}
    }
printf("--- %d different values were used\n", n);
// Now do the remapping
for(int i=0; i<w*h; i++)
    if( test[i] != 0 )
        test[i] = SPmapping[test[i]] - base;  // cannot overflow, since cannot have more than 65535 new values
}

/* --------------------------------------------------------------- */
/* RemapSuperPixelsOneImage -------------------------------------- */
/* --------------------------------------------------------------- */

// Same thing, but with ints. Could templatize this.
//
static void RemapSuperPixelsOneImage(
	uint32		*test,
	int			w,
	int			h,
	int			&MaxSPUsed,
	vector<int>	&SPmapping )
{
// First find the biggest
int biggest = -1;
for(int i=0; i<w*h; i++)
    biggest = max(int(test[i]), biggest);
printf("Biggest value in 32 bit map is %d\n", biggest);
//for(int i=0; i<w*h; i++) {
    //if( int(test[i]) == biggest )
	//printf("Occurs: test[%d] = %d\n", i, test[i]);
    //}

// Find which numbers were used.
vector<int> vals(biggest+1,0);
for(int i=0; i<w*h; i++)
    vals[test[i]]++;

SPmapping.resize(biggest+1,0);  // make the remapping vector the right size
int base = MaxSPUsed;  // will add this to all
int n=0;
for(int i=1; i<=biggest; i++) {  // 0 always maps to 0
    if( vals[i] != 0 ) { // this value was used
	MaxSPUsed++;
        SPmapping[i] = MaxSPUsed;
	n++;
	}
    }
printf("--- %d different values were used\n", n);
// Now do the remapping.  'Base' not currently used for 32 bit version
if( base != 0 ) {
    printf("Expected base to be 0, not %d\n", base);
    exit( 42 );
    }
for(int i=0; i<w*h; i++)
    test[i] = SPmapping[test[i]] - base;
}

/* --------------------------------------------------------------- */
/* RuleOutIntersection ------------------------------------------- */
/* --------------------------------------------------------------- */

// Can we rule out that an image from (0,0) to (w-1,h-1),
// when transformed by any of its 'tfs', will hit the image
// (xmin,ymin) to (xmax,ymax)? Try to find a separating line...
//
static bool RuleOutIntersection(
	uint32			w,
	uint32			h,
	vector<TAffine>	&tfs,
	int				xmin,
	int				ymin,
	int				xmax,
	int				ymax )
{
// make the four corners
vector<Point> corners;
corners.push_back(Point(0.0, 0.0));
corners.push_back(Point(w-1, 0.0));
corners.push_back(Point(w-1, h-1));
corners.push_back(Point(0.0, h-1));
// Make a list of all corners, transformed by all transforms
vector<Point> all;
printf("xmin %d, ymin %d, xmax %d, ymax %d, w %d, h %d\n", xmin, ymin, xmax, ymax, w, h);
for(int j=1; j<tfs.size(); j++) {  // patches start at 1
    tfs[j].TPrint();
    for(int i=0; i < corners.size(); i++) {
        Point p = corners[i];
        tfs[j].Transform( p );
        printf("Now (%f %f)\n", p.x, p.y);
        all.push_back(p);
	}
    }
// Now check against each edge of the viewing region
bool AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x < xmin);
if( AllOut ) return true;          // all on left of view
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x > xmax);
if( AllOut ) return true;          // all on right
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].y < ymin);
if( AllOut ) return true;          // all on bottom
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].y > ymax);
if( AllOut ) return true;          // all off the top
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x + all[i].y > xmax + ymax);
if( AllOut ) return true;          // all off northeast
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x + all[i].y < xmin + ymin);
if( AllOut ) return true;          // all off southwest
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x - all[i].y < xmin - ymax);
if( AllOut ) return true;          // all off northwest
AllOut = true;
for(int i=0; i<all.size(); i++)
    AllOut &= (all[i].x - all[i].y > xmax - ymin);
if( AllOut ) return true;          // all off southeast
return false;
}

/* --------------------------------------------------------------- */
/* ImageResize --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Routines to re-size an image down by a factor of 'scale';
// ww and hh are the original size; scale is the integer factor.
//
static void ImageResize(
	uint8	*&raster,
	int		ww,
	int		hh,
	int		scale,
	bool	sampling = true )
{
if( scale == 1 )
    return;
int w = ww/scale;
int h = hh/scale;
uint8* new_ras = (uint8 *)malloc(w*h*sizeof(uint8));
if( sampling ) { // just pick one of the pixels
    for(int y=0; y<h; y++)
	for(int x=0; x<w; x++)
	    new_ras[x + w*y] = raster[scale*(x + ww*y)];
    }
else {
    for(int y=0; y<h; y++) {
	int y0 = y*scale;
	for(int x=0; x<w; x++) {
	    int x0 = x*scale;
	    int sum = 0;
	    for(int i=0; i<scale; i++) {
		for(int j=0; j < scale; j++)
		    sum += raster[x0+x + ww*(y0+y)];
		}
	    new_ras[x + w*y] = sum/(scale*scale);
	    }
	}
    }
free(raster);
raster = new_ras;
}

/* --------------------------------------------------------------- */
/* ImageResize16bits --------------------------------------------- */
/* --------------------------------------------------------------- */

// Only used for super-pixels, so only sampling and not averaging.
//
static void ImageResize16bits(
	uint16	*&raster,
	int		ww,
	int		hh,
	int		scale )
{
if( scale == 1 )
    return;
int w = ww/scale;
int h = hh/scale;
uint16* new_ras = (uint16 *)malloc(w*h*sizeof(uint16));
for(int y=0; y<h; y++)
    for(int x=0; x<w; x++)
	new_ras[x + w*y] = raster[scale*(x + ww*y)];
free(raster);
raster = new_ras;
}

/* --------------------------------------------------------------- */
/* PointsInRing -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Produce the points at radius R (a square) around the center,
// in steps < 1 in size. Might include duplicates.
//
static void PointsInRing(
	vector<Point>	&pts,
	const Point		&center,
	double			r,
	double			delta,
	uint32			w,
	uint32			h )
{
	double	lx = center.x - r;
	double	rx = center.x + r;
	double	by = center.y - r;
	double	ty = center.y + r;
	vector<Point>	ppts; // potential points

	for( double d = -r; d <= r; d += delta ) {

		ppts.push_back( Point( lx, center.y + d ) );
		ppts.push_back( Point( rx, center.y + d ) );
		ppts.push_back( Point( center.x + d, by ) );
		ppts.push_back( Point( center.x + d, ty ) );
	}

// Now push back upper right corner, which might not be covered.

	ppts.push_back( Point( rx, ty ) );

// Now keep only those within the image
// (can happen if image is a rectangle, not a square)

	int	n = ppts.size();

	for( int i = 0; i < n; ++i ) {

		if( ppts[i].x >= 0.0 &&
			ppts[i].x <= w-1 &&
			ppts[i].y >= 0.0 &&
			ppts[i].y <= h-1 ) {

			pts.push_back( ppts[i] );
		}
	}
}

/* --------------------------------------------------------------- */
/* FillInFromCenterOut ------------------------------------------- */
/* --------------------------------------------------------------- */

// Routine to fill in the map, starting with the center of each
// image and filling outward. This sets each pixel based on how
// close it is the center, avoiding the edges where the biggest
// distortions are.
//
static void FillInFromCenterOut(
	vector<uint16>	&imap,
	vector<uint8>	&PatchFrom,
	int				nx,
	int				ny,
	int				w,
	int				h,					// size of images coming in
	double			delta_image_space,	// smallest size that maps to each pixel on output
	vector<image>	&images,			// the images
	vector<int>		&relevant_images )	// the ones that are used among them
{
// find the center of each picture.  May not be integer
Point center((w-1)/2.0, (h-1)/2.0);
int report = 100;
for(double r=0.0; r <= max(center.x,center.y)+0.0001; r += delta_image_space) {
    if( r >= report ) {
	printf("Starting radius %f\n", r);
	report += 100;
	}
    //compute a ring of radius R
    vector<Point> pts;
    PointsInRing(pts, center, r, delta_image_space, w, h);
    // Now transform this by all relevant images
    for(int k=0; k<relevant_images.size(); k++) {      // for each picture
	int i = relevant_images[k];
	for(int j=0; j<pts.size(); j++) {
	    Point p(pts[j]);
	    int ix = int(p.x);  // PointsInRing only returns legal points
	    int iy = int(p.y);
	    int patch = images[i].foldmap[ix + w*iy];
	    // There can be small patches which were not used, so need to check for legality here
	    // No longer needed since we compress patches and foldmaps on input
	    if( patch == 0 )  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
		continue;
	    images[i].tf[patch].Transform( p );  // change to global coordinates
	    ix = ROUND(p.x);
	    iy = ROUND(p.y);
	    if( ix < 0 || ix >= nx || iy < 0 || iy >= ny )
		continue;  // outside the image
	    //printf("%f i=%d ix,iy=%d %d\n", r, i, ix, iy);
	    uint32 nn = ix + uint32(iy)*uint32(nx);   // index into array
	    if( imap[nn] == 0 ) {
		imap[nn] = k+1;       // this pixel will be set by the kth relevant picture
		PatchFrom[nn] = patch;
		}
	    else { // already set, but we still might be better (because of finite
		   // grid in source, rounding to int, ordering ).
		Point pt_us = Point(ix,iy);
		Point pt_ot = pt_us;
		// transform back into original frames.
		int oi = relevant_images[imap[nn]-1];     // other image
		int op = PatchFrom[nn]; // other patch
		images[oi].inv[op].Transform( pt_ot );
		images[i].inv[patch].Transform( pt_us );

		double d_us = max(abs(pt_us.x-center.x), abs(pt_us.y-center.y));  // L-infinity norm
		double d_ot = max(abs(pt_ot.x-center.x), abs(pt_ot.y-center.y));  // L-infinity norm
		if( d_us < d_ot ) {
		    //printf("It CAN happen... d_us= %f, d_ot= %f\n", d_us, d_ot);
		    imap[nn] = k+1;       // this pixel will be set by the ith picture
		    PatchFrom[nn] = patch;
		    }
		}
	    }
	}
    }
}

/* --------------------------------------------------------------- */
/* NameSorter ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// An alternate to the above routine, this one fills in as the
// Matlab version did, to help with re-preducing earlier Matlab
// results for ease of back-annotation of synapse and other
// annotations. Here each tile is completely rendered in alphabetical
// order (emperically true of the XML files matlab used for ordering).
//
class NameSorter{
public:
	char	*name;
	int		num;
public:
	bool operator < ( const NameSorter &rhs ) const
		{return strcmp( name, rhs.name ) < 0;};
};

/* --------------------------------------------------------------- */
/* NameSortFn ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool NameSortFn( NameSorter a, NameSorter b )
{
	return strcmp( a.name, b.name ) > 0;
}

/* --------------------------------------------------------------- */
/* FillInMatlabOrder --------------------------------------------- */
/* --------------------------------------------------------------- */

static void FillInMatlabOrder(
	vector<uint16>	&imap,
	vector<uint8>	&PatchFrom,
	int				nx,
	int				ny,
	int				w,
	int				h,					// size of images coming in
	double			delta_image_space,	// smallest size that maps to each pixel on output
	vector<image>	&images,			// the images
	vector<int>		&relevant_images)	// the ones that are used among them
{
// find the center of each picture.  May not be integer
Point center((w-1)/2.0, (h-1)/2.0);
// sort the relevent images (null for now)
vector<NameSorter> temp(relevant_images.size());
for(int k=0; k<relevant_images.size(); k++) {      // for each picture
    int i = relevant_images[k];
    temp[k].num  = k;                              // where this name is in the relevant picture vector
    temp[k].name = images[i].rname;
    }
sort(temp.begin(), temp.end(), NameSortFn);    // sort in decreasing order

// OK, render in order
for(int kk=0; kk<temp.size(); kk++) {      // for each picture
    printf("Start image %s\n", temp[kk].name);
    int k = temp[kk].num;
    int i = relevant_images[k];
    int report = 500;
    for(double r=0.0; r <= max(center.x,center.y)+0.0001; r += delta_image_space) {
	if( r >= report ) {
	    printf("Starting radius %f\n", r);
	    report += 500;
	    }
	//compute a ring of radius R
	vector<Point> pts;
	PointsInRing(pts, center, r, delta_image_space, w, h);
	// Now transform this by all relevant images
	for(int j=0; j<pts.size(); j++) {
	    Point p(pts[j]);
	    int ix = int(p.x);  // PointsInRing only returns legal points
	    int iy = int(p.y);
	    int patch = images[i].foldmap[ix + w*iy];
	    // There can be small patches which were not used, so need to check for legality here
	    // No longer needed since we compress patches and foldmaps on input
	    if( patch == 0 )  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
		continue;
	    images[i].tf[patch].Transform( p );  // change to global coordinates
	    ix = ROUND(p.x);
	    iy = ROUND(p.y);
	    if( ix < 0 || ix >= nx || iy < 0 || iy >= ny )
		continue;  // outside the image
	    //printf("%f i=%d ix,iy=%d %d\n", r, i, ix, iy);
	    uint32 nn = ix + uint32(iy)*uint32(nx);   // index into array
	    if( imap[nn] == 0 ) {
		imap[nn] = k+1;       // this pixel will be set by the kth relevant picture
		PatchFrom[nn] = patch;
		}
	    }
	}
    }
}

/* --------------------------------------------------------------- */
/* MakeDirExist -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Makes the directory exist.
// If it already exists and is not a directory, or cannot be created,
// print an error message and exit with an error code.
// If it does exist and is a directory, or if we can create it,
// return normally.
//
static void MakeDirExist( const char *name, set<string> &dir_cache )
{
// first, see if we already created this directory.
if( dir_cache.find(string(name)) != dir_cache.end() )
    return;

// OK, have not created it already.  One pass occasionally fails when running many
// parallel jobs.  I suspect the stat fails (dir does not exist) but by the time
// the makedir() is called, it *does* exist, and the mkdir fails.  So try PASS times
const int MAX_PASS = 3;
int pass;
for(pass=1; pass<=MAX_PASS; pass++) {
    struct stat buf;
    if( stat(name, &buf) ) {
        // error, assume directory does not exist
        if( mkdir(name,0775) )
	    printf("mkdir('%s') failed.\n", name);
        else // it worked, stop trying
	    break;
        }
    else {  // was able to stat the file; see if it's a directory
        if( !S_ISDIR(buf.st_mode) ) {
	    printf("File '%s' exists, but is not a directory.\n", name);
	    exit( 42 );
	    }
        else
	    break;  // file exists and is a directory
	}
    }
if( pass > MAX_PASS )
    exit( 42 );
// if we get here, either the directory already exists, or we created it.  In
// either case, no need to make it again.
dir_cache.insert(string(name));
}

/* --------------------------------------------------------------- */
/* CopyRaster ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// For images...
// Copy and condense the raster at src into the location in dest
// (dest, w, h) = describes dest
// (dx, dy)     = where to write
// (src, sw, sh)= describes the source
// Only problem is that sw, or sh, or both could be odd,
// so each pixel can get a contribution from 1, 2, or 4 pixels.
//
static void CopyRaster(
	uint8	*dest,
	int		w,
	int		h,
	int		dx,
	int		dy,
	uint8	*src,
	int		sw,
	int		sh )
{
// Do this a brute force way.  Go through all the src pixels, mapping each to dest,
// and counting how many map to each.  Then average.
// Note added later:  This is now too general, as Phil says the tile size always
// rounds down, so don't worry about the half pixels.  Could now make make more efficient.
if( src == NULL )
    return;
vector<int> N(w*h,0);
vector<int> d(w*h,0);
int swt = sw & (~1);  // create copies that are rounded down to nearest even
int sht = sh & (~1);  // stands for "source height truncated"
for(int y=0; y<sht; y++) {
    for(int x=0; x<swt; x++) {
	int xx = dx + x/2;
        int yy = dy + y/2;
        d[xx + w*yy] += src[x + sw*y];
        N[xx + w*yy]++;
	}
    }
// now, for every pixel that was written at least once, average it
for(int x = 0; x<w; x++)
    for(int y=0; y<h; y++)
	if( N[x + w*y] != 0 )
	    dest[x + w*y] = d[x +w*y]/N[x+w*y];
}

/* --------------------------------------------------------------- */
/* CopyRasterSP -------------------------------------------------- */
/* --------------------------------------------------------------- */

// For SuperPixels...
// Copy and condense the raster at src into the location in dest
// (dest, w, h) = describes dest
// (dx, dy)     = where to write
// (src, sw, sh)= describes the source
// Only problem is that sw, or sh, or both could be odd,
// so each pixel an get a contribution from 1, 2, or 4 pixels.
// Do NOT want interpolation; instead take first one.
//
static void CopyRasterSP(
	uint32			*dest,
	int				w,
	int				h,
	int				dx,
	int				dy,
	const uint32	*src,
	int				sw,
	int				sh )
{
// Do this a brute force way.  Go through all the src pixels, mapping each to dest,
// and counting how many map to each.  Then average.
// Note added later:  This is now too general, as Phil says the tile size always
// rounds down, so don't worry about the half pixels.  Could now make make more efficient.
// Bizzare convoluted logic since copied from image, which needed to average.
if( src == NULL )
    return;
vector<int> N(w*h,0);
vector<int> d(w*h,0);
int swt = sw & (~1);  // create copies that are rounded down to nearest even
int sht = sh & (~1);  // stands for "source height truncated"
for(int y=0; y<sht; y++) {
    for(int x=0; x<swt; x++) {
	int xx = dx + x/2;
        int yy = dy + y/2;
        d[xx + w*yy] = src[x + sw*y];  // just replace; no average
        N[xx + w*yy]++;
	}
    }
// now, for every pixel that was written at least once, set it
for(int x = 0; x<w; x++)
    for(int y=0; y<h; y++)
	if( N[x + w*y] != 0 )
	    dest[x + w*y] = d[x +w*y];
}

/* --------------------------------------------------------------- */
/* WriteSummaryTile ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Condenses 4 tiles into 1, for images.
//
static int WriteSummaryTile(
	string		rav_name,
	int			section,
	int			level,
	int			row,
	int			col,
	set<string>	&dir_cache )
{
// open the 4 relevant tiles - (row,col), (row+1,col), (row,col+1), and (row+1,col+1)
uint32 w1, h1, w2, h2, w3, h3, w4, h4;
char fname[2048];
int dir = section/1000;
char ds[32];                     // string for the final directory
strcpy(ds,"");                   // for layers 0-999, it's blank (no directory)
if( dir != 0 )                    // but otherwise, one per thousand layers
    sprintf(ds, "%d/", dir*1000);// for example, layer 1761 would be in directory "1000/"

//Tiles are named after compass directions on the screen (sw= southwest, for example)
sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%s%03d.png", rav_name.c_str(), level-1, row, col, ds, section);
uint8* sw = Raster8FromPng(fname, w1, h1);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%s%03d.png", rav_name.c_str(), level-1, row+1, col, ds, section);
uint8* nw = Raster8FromPng(fname, w2, h2);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%s%03d.png", rav_name.c_str(), level-1, row, col+1, ds, section);
uint8* se = Raster8FromPng(fname, w3, h3);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%s%03d.png", rav_name.c_str(), level-1, row+1, col+1, ds, section);
uint8* ne = Raster8FromPng(fname, w4, h4);
if( sw == NULL )  // if even the base tile does not exist, don't write anything.
    return 0;
// if we are at the origin, and one tile holds everything, we don't need a tile either
if( nw == NULL && se == NULL && row == 0 & col == 0 )
    return 0;

// otherwise, create the tile...
int new_w = w1 + (se != NULL ? w3 : 0);
int new_h = h1 + (nw != NULL ? h2 : 0);
// now we need to divide these by 2, rounding down according to spec
int w = new_w/2;
int h = new_h/2;
uint8* raster = (uint8 *)malloc(w*h*sizeof(uint8));
// copy each of the rasters into its spot.
int half = 1024/2;
CopyRaster(raster, w, h, 0,    (new_h-h1)/2, sw,   w1, h1);
CopyRaster(raster, w, h, 0,    0, nw , w2, h2);
CopyRaster(raster, w, h, half, (new_h-h1)/2, se , w3, h3);
CopyRaster(raster, w, h, half, 0, ne, w4, h4);
// Make sure the level, the row, and the colum exist
sprintf(fname, "%s/tiles/1024/%d", rav_name.c_str(), level);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d", rav_name.c_str(), level, row/2);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d/%d", rav_name.c_str(), level, row/2, col/2);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/g", rav_name.c_str(), level, row/2, col/2);
MakeDirExist(fname, dir_cache);
if( dir == 0 ) // first layers are written in the root directory
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%03d.png", rav_name.c_str(), level, row/2, col/2, section);
else { // need another layer of directories
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%d", rav_name.c_str(), level, row/2, col/2, dir*1000);
    MakeDirExist(fname, dir_cache);
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/g/%d/%03d.png", rav_name.c_str(), level, row/2, col/2, dir*1000, section);
    }
printf("writing tile '%s', %d by %d\n", fname, w, h);
Raster8ToPng8( fname, raster, w, h );
free(raster);
return 1;
}

/* --------------------------------------------------------------- */
/* WriteSummaryTileSP -------------------------------------------- */
/* --------------------------------------------------------------- */

// Condenses 4 tiles into 1, for SuperPixels.
//
static int WriteSummaryTileSP(
	string		rav_name,
	int			section,
	int			level,
	int			row,
	int			col,
	set<string>	&dir_cache )
{
// open the 4 relevant tiles - (row,col), (row+1,col), (row,col+1), and (row+1,col+1)
uint32 w1, h1, w2, h2, w3, h3, w4, h4;
char fname[2048];
int dir = section/1000;
char ds[32];                     // string for the final directory
strcpy(ds,"");                   // for layers 0-999, it's blank (no directory)
if( dir != 0 )                    // but otherwise, one per thousand layers
    sprintf(ds, "%d/", dir*1000);// for example, layer 1761 would be in directory "1000/"

//Tiles are named after compass directions on the screen (sw= southwest, for example)
sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%s%03d.png", rav_name.c_str(), level-1, row, col, ds, section);
uint32* sw = Raster32FromPng(fname, w1, h1);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%s%03d.png", rav_name.c_str(), level-1, row+1, col, ds, section);
uint32* nw = Raster32FromPng(fname, w2, h2);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%s%03d.png", rav_name.c_str(), level-1, row, col+1, ds, section);
uint32* se = Raster32FromPng(fname, w3, h3);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%s%03d.png", rav_name.c_str(), level-1, row+1, col+1, ds, section);
uint32* ne = Raster32FromPng(fname, w4, h4);
if( sw == NULL )  // if even the base tile does not exist, don't write anything.
    return 0;
// if we are at the origin, and one tile holds everything, we don't need a tile either
if( nw == NULL && se == NULL && row == 0 & col == 0 )
    return 0;

// otherwise, create the tile...
int new_w = w1 + (se != NULL ? w3 : 0);
int new_h = h1 + (nw != NULL ? h2 : 0);
// now we need to divide these by 2, rounding down according to spec
int w = new_w/2;
int h = new_h/2;
uint32* raster = (uint32 *)malloc(w*h*sizeof(uint32));
// copy each of the rasters into its spot.
int half = 1024/2;
CopyRasterSP(raster, w, h, 0,    (new_h-h1)/2, sw,   w1, h1);
CopyRasterSP(raster, w, h, 0,    0, nw , w2, h2);
CopyRasterSP(raster, w, h, half, (new_h-h1)/2, se , w3, h3);
CopyRasterSP(raster, w, h, half, 0, ne, w4, h4);
// Make sure the level, the row, and the column exist
sprintf(fname, "%s/tiles/1024/%d", rav_name.c_str(), level);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d", rav_name.c_str(), level, row/2);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d/%d", rav_name.c_str(), level, row/2, col/2);
MakeDirExist(fname, dir_cache);
sprintf(fname, "%s/tiles/1024/%d/%d/%d/s", rav_name.c_str(), level, row/2, col/2);
MakeDirExist(fname, dir_cache);
if( dir == 0 ) // first layers are written in the root directory
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%03d.png", rav_name.c_str(), level, row/2, col/2, section);
else { // need another layer of directories
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%d", rav_name.c_str(), level, row/2, col/2, dir*1000);
    MakeDirExist(fname, dir_cache);
    sprintf(fname, "%s/tiles/1024/%d/%d/%d/s/%d/%03d.png", rav_name.c_str(), level, row/2, col/2, dir*1000, section);
    }
printf("SP:writing tile '%s', %d by %d\n", fname, w, h);
Raster32ToPngRGBA( fname, raster, w, h );
free(raster);
return 1;
}

/* --------------------------------------------------------------- */
/* CreateLevelNTiles --------------------------------------------- */
/* --------------------------------------------------------------- */

//Create smaller tiles, for images.
//
static void CreateLevelNTiles(
	string		rav_name,
	int			section,
	set<string>	&dir_cache )
{
int at_this_level = 1;  // so we loop at least once
for(int level=1; level < 20 && at_this_level > 0; level++) {  // enough for trillions of tiles
    at_this_level = 0;
    int wrote = 1;
    // now do this level
    for(int row=0; wrote > 0; row += 2) {
        wrote = 0;
        int stat = 1;
        for(int col = 0; stat == 1; col += 2) {
            printf("Trying to write tile for level %d, row %d, col %d\n", level, row, col);
            stat = WriteSummaryTile(rav_name, section, level, row, col, dir_cache);
            wrote += stat;
            }
        at_this_level += wrote;
	}
    printf("At level %d, wrote %d tiles\n", level, at_this_level);
    }
}

/* --------------------------------------------------------------- */
/* CreateLevelNTilesSP ------------------------------------------- */
/* --------------------------------------------------------------- */

//Create smaller tiles, for SuperPixels.
//
void CreateLevelNTilesSP(
	string		rav_name,
	int			section,
	set<string>	&dir_cache )
{
int at_this_level = 1;  // so we loop at least once
for(int level=1; level < 20 && at_this_level > 0; level++) {  // enough for trillions of tiles
    at_this_level = 0;
    int wrote = 1;
    // now do this level
    for(int row=0; wrote > 0; row += 2) {
        wrote = 0;
        int stat = 1;
        for(int col = 0; stat == 1; col += 2) {
            printf("SP:Trying to write tile for level %d, row %d, col %d\n", level, row, col);
            stat = WriteSummaryTileSP(rav_name, section, level, row, col, dir_cache);
            wrote += stat;
            }
        at_this_level += wrote;
	}
    printf("SP:At level %d, wrote %d tiles\n", level, at_this_level);
    }
}

/* --------------------------------------------------------------- */
/* UpdateMetaData ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Updates the meta-data file. For now, just writes.
//
static void UpdateMetaData(
	const char	*name,
	int			section,
	int			w,
	int			h )
{
	int		zmin, zmax;
	FILE	*fd = fopen( name, "r" );

	if( fd ) {

		CLineScan	LS;
		ssize_t		read;

		while( (read = LS.Get( fd )) > 0 ) {

			//printf( "Retrieved line of length %zu :\n", read );

			if( !strncmp( LS.line, "zmin=", 5 ) )
				zmin = min( section, atoi(LS.line+5) );
			if( !strncmp( LS.line, "zmax=", 5 ) )
				zmax = max( section, atoi(LS.line+5) );
		}

		fclose( fd );
	}
	else {	// no file yet; make one with one layer
		zmin = section;
		zmax = section;
	}

// Now write the new file

	fd = FileOpenOrDie( name, "w" );

	fprintf( fd, "version=1\n" );
	fprintf( fd, "width=%d\n", w );
	fprintf( fd, "height=%d\n", h );
	fprintf( fd, "zmin=%d\n", zmin );
	fprintf( fd, "zmax=%d\n", zmax );
	fprintf( fd, "superpixel-format=RGBA\n" );
	fprintf( fd, "channels=\"g;s\"\n" );

	fclose( fd );
}

/* --------------------------------------------------------------- */
/* WriteImageTiles ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Writes the lowest level image tiles.
//
static int WriteImageTiles(
	const string	rav_name,
	int				section,
	const uint8		*image,
	int				w,
	int				h )
{
set<string> dir_cache;  // names of directories already created
string base_name = rav_name;
base_name += "/tiles";
MakeDirExist(base_name.c_str(), dir_cache);    // create the directory 'tiles'
char fname[2048];
sprintf(fname, "%s/metadata.txt", base_name.c_str());
UpdateMetaData(fname, section, w, h);
base_name += "/1024";
MakeDirExist(base_name.c_str(), dir_cache);    // we will always make tiles 1024 on a side
base_name += "/0";
MakeDirExist(base_name.c_str(), dir_cache);    // and always level 0 (raw data)
uint8* raster = (uint8 *)malloc(1024*1024*sizeof(uint8));  // make a raster big enough for even the biggest tile
for(int row = 0; row*1024 < h; row++) {
    int ymax = h-1-1024*row;  // top scan line
    int ymin = max(ymax - 1024+1, 0);             // bottom scan line
    sprintf(fname, "%s/%d", base_name.c_str(), row);
    MakeDirExist(fname, dir_cache);
    for(int col = 0; col*1024 < w; col++) {
	int xmin = col*1024;             // bottom scan line
	int xmax = min(xmin+1023, w-1);  // top scan line
        sprintf(fname, "%s/%d/%d", base_name.c_str(), row, col);      // make the dir for the column
	MakeDirExist(fname, dir_cache);
        sprintf(fname, "%s/%d/%d/g", base_name.c_str(), row, col);    // and the grey scale underneath
	MakeDirExist(fname, dir_cache);
        // now, copy the portion of the image from xmin..xmax to ymin..ymax;
        int new_w = xmax - xmin + 1;
        int n = 0;
	for(int y=ymin; y<=ymax; y++) {
            uint32 m = uint32(xmin) + uint32(y)*uint32(w); // location of first byte on scan line y of 'image'
	    for(int x=0; x < new_w; x++)
		raster[n++] = image[m++];
	    }
        int dir = section/1000;
        if( dir == 0 ) // first layers are written in the root directory
	    sprintf(fname,"%s/%d/%d/g/%03d.png", base_name.c_str(), row, col, section);
        else { // need another layer of directories
	    sprintf(fname,"%s/%d/%d/g/%d", base_name.c_str(), row, col, dir*1000);
	    MakeDirExist(fname, dir_cache);
	    sprintf(fname,"%s/%d/%d/g/%d/%03d.png", base_name.c_str(), row, col, dir*1000, section);
	    }
        Raster8ToPng8( fname, raster, new_w, ymax-ymin+1 );
	}
    }
free(raster);

// Now create the higher level tiles.
CreateLevelNTiles(rav_name, section, dir_cache);
return 0;
}

/* --------------------------------------------------------------- */
/* WriteSPTiles -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Writes the lowest level SP tiles.
//
static int WriteSPTiles(
	const string	rav_name,
	int				section,
	uint32			*image,
	int				w,
	int				h )
{
set<string> dir_cache;  // names of directories already created
string base_name = rav_name;
base_name += "/tiles";
MakeDirExist(base_name.c_str(), dir_cache);    // create the directory 'tiles'
char fname[2048];
sprintf(fname, "%s/metadata.txt", base_name.c_str());
//UpdateMetaData(fname, section, w, h);   image tiles already did this
base_name += "/1024";
MakeDirExist(base_name.c_str(), dir_cache);    // we will always make tiles 1024 on a side
base_name += "/0";
MakeDirExist(base_name.c_str(), dir_cache);    // and always level 0 (raw data)
uint32* raster = (uint32 *)malloc(1024*1024*sizeof(uint32));  // make a raster big enough for even the biggest tile
for(int row = 0; row*1024 < h; row++) {
    int ymax = h-1-1024*row;  // top scan line
    int ymin = max(ymax - 1024+1, 0);             // bottom scan line
    sprintf(fname, "%s/%d", base_name.c_str(), row);
    MakeDirExist(fname, dir_cache);
    for(int col = 0; col*1024 < w; col++) {
	int xmin = col*1024;             // bottom scan line
	int xmax = min(xmin+1023, w-1);  // top scan line
        sprintf(fname, "%s/%d/%d", base_name.c_str(), row, col);      // make the dir for the column
	MakeDirExist(fname, dir_cache);
        sprintf(fname, "%s/%d/%d/s", base_name.c_str(), row, col);    // and the sp map underneath
	MakeDirExist(fname, dir_cache);
        // now, copy the portion of the image from xmin..xmax to ymin..ymax;
        int new_w = xmax - xmin + 1;
        int n = 0;
	for(int y=ymin; y<=ymax; y++) {
            uint32 m = uint32(xmin) + uint32(y)*uint32(w); // location of first byte on scan line y of 'image'
	    for(int x=0; x < new_w; x++)
		raster[n++] = image[m++];
	    }
        int dir = section/1000;
        if (dir == 0) // first layers are written in the root directory
	    sprintf(fname,"%s/%d/%d/s/%03d.png", base_name.c_str(), row, col, section);
        else { // need another layer of directories
	    sprintf(fname,"%s/%d/%d/s/%d", base_name.c_str(), row, col, dir*1000);
	    MakeDirExist(fname, dir_cache);
	    sprintf(fname,"%s/%d/%d/s/%d/%03d.png", base_name.c_str(), row, col, dir*1000, section);
	    }
        Raster32ToPngRGBA( fname, raster, new_w, ymax-ymin+1 );
	}
    }
free(raster);

// Now create the higher level tiles.
CreateLevelNTilesSP(rav_name, section, dir_cache);
return 0;
}

/* --------------------------------------------------------------- */
/* Annotation Inversion ------------------------------------------ */
/* --------------------------------------------------------------- */

#include	"json/json.h"

class Bookmark {
public:
	string	text;
	Point	pt;
	int		z;
	int		body_id;
	int		image_number;
public:
	bool operator < ( const Bookmark &rhs ) const
		{return image_number < rhs.image_number;};
};


class Partner {
public:
	double	confidence;
	Point	pt;
	int		body_id;
	int		z;
};


class TBar {
public:
	string			status;
	double			confidence;
	Point			pt;
	int				z;
	vector<Partner>	partners;
	int				body_id;
	int				image_number;
public:
	bool operator < ( const TBar &rhs ) const
		{return image_number < rhs.image_number;};
};

/* --------------------------------------------------------------- */
/* ReadBookmarks ------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool ReadBookmarks(
	char				*fname,
	int					section,
	vector<Bookmark>	&bookmarks )
{
Json::Value root;   // will contains the root value after parsing.
Json::Reader reader;
FILE *fd = fopen( fname, "r" );
if( fd == NULL ) {
    printf("Cannot open '%s' for read.\n", fname);
    return false;
    }
const int N=1000000;
char huge[N];
fread(huge, 1, N, fd);
fclose(fd);
bool parsingSuccessful = reader.parse(huge, root );
if( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    exit( 42 );
    }

printf("parsed file...\n");
Json::Value md = root["metadata"];
string des = md["description"].asString();
int         ver = md["version"].asInt();

printf("Description: '%s', version %d\n", des.c_str(), ver );
// Go through all the bookmarks, keeping only those on the specified layers
Json::Value data = root["data"];
int	nd = data.size();
for( unsigned int i = 0; i < nd; ++i ) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only is known
    Bookmark b;
    b.text    = v["text"].asString();
    b.body_id = v["body ID"].isInt();
    Json::Value lo = v["location"];
    b.z = lo[2u].asInt();
    if( b.z == section ) {
	b.pt.x = lo[0u].asDouble();
	b.pt.y = lo[1u].asDouble();
        bookmarks.push_back(b);
	}
    }
return true;
}

/* --------------------------------------------------------------- */
/* ReadSynAnnot -------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool ReadSynAnnot(
	char			*fname,
	int				section,
	vector<TBar>	&tbs )
{
Json::Value root;   // will contains the root value after parsing.
Json::Reader reader;
FILE *fd = fopen( fname, "r" );
if( fd == NULL ) {
    printf("Cannot open '%s' for read.\n", fname);
    return false;
    }
const int N=1000000;
char huge[N];
fread(huge, 1, N, fd);
fclose(fd);
bool parsingSuccessful = reader.parse(huge, root );
if( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    exit( 42 );
    }

printf("parsed file...\n");
Json::Value md = root["metadata"];
string des = md["description"].asString();
int         ver = md["version"].asInt();

printf("Description: '%s', version %d\n", des.c_str(), ver );
// Go through all the T-bars, keeping only those on the specified layers
Json::Value data = root["data"];
int	nd = data.size();
for( unsigned int i = 0; i < nd; ++i ) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only {"T-bar","partner"} is known
    Json::Value t = v["T-bar"];
    Json::Value part = v["partners"];
    if( !t.isNull() ) {
        Json::Value lo = t["location"];
        TBar tb;
        tb.z = lo[2u].asInt();
        if( tb.z == section ) {
	    tb.status      = t["status"].asString();
	    tb.confidence  = t["confidence"].asDouble();
	    tb.body_id     = t["body ID"].asInt();
	    tb.pt.x = lo[0u].asDouble();
	    tb.pt.y = lo[1u].asDouble();
	    //printf("status is %s, confidence %f, loc %.1f,%lf,%d,  %d partners\n",
	    //tb.status.c_str(), tb.confidence, tb.pt.x, tb.pt.y, tb.z, part.size() );
	    for(int j=0; j<part.size(); j++) {
		Partner prt;
		prt.confidence = part[j]["confidence"].asDouble();
		prt.body_id    = part[j]["body ID"   ].asInt();
		Json::Value lo = part[j]["location"];
		prt.pt.x = lo[0u].asDouble();
		prt.pt.y = lo[1u].asDouble();
		prt.z = lo[2u].asInt();
		//printf("   Partner %d, confidence %f. loc=%f %f %d\n", j, prt.confidence, prt.pt.x, prt.pt.y, prt.z);
		tb.partners.push_back(prt);
		}
	    tbs.push_back(tb);
	    }
	}
    }
return true;
}

/* --------------------------------------------------------------- */
/* LookForAndTransform ------------------------------------------- */
/* --------------------------------------------------------------- */

// Read a file of synapse annotations. Find all the ones that
// map to the final image, transform them into final image
// coordinates, then append them to the vector result[].
//
static bool LookForAndTransform(
	const char		*name,
	char			*root,
	int				section,
	int				image_no,
	vector<image>	&images,
	bool			Warp,
	vector<uint8>	&tmap,
	vector<TBar>	&result )
{
uint32 w = images[image_no].w;
uint32 h = images[image_no].h;
char lname[2048];
sprintf(lname, "%s/%s", root, name);
FILE *fp = fopen( lname, "r" );
if( fp == NULL ) {
    printf("No file '%s'\n", lname);
    return false;
    }
printf("Found '%s', copying annotation.\n", lname);
fclose(fp);
if( strstr(name, "annotations-synapse") != NULL ) {
    vector<TBar> tbs;
    ReadSynAnnot(lname, section, tbs);
    // Use the location to find the patch
    for(int j=0; j<tbs.size(); j++) {
        // OK, find the consistent transform that maps into final image
        Point p = tbs[j].pt;
        int ix = int(p.x);
        int iy = int(p.y);
	// make sure it's legal
	if( ix < 0 ||  ix >= w || iy < 0 || iy >= h ) {
	    printf("Synapse outside diagram?? %d %d %d %d\n", ix, w, iy, h);
	    continue;
	    }
        int patch = images[image_no].foldmap[ix + w*iy];
        // There can be small patches which were not used, so need to check for legality here
        //             // No longer needed since we compress patches and foldmaps on input
        if( patch == 0 )  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
            continue;
        TAffine t;
        if( !Warp )  // just one transform per image
            t = images[image_no].tf[patch];  // change to global coordinates
        else { //we are warping
            int sector = tmap[ix + w*iy];
            t = images[image_no].sectors[patch][sector];
	    }
        t.Transform( p );
        tbs[j].pt = p;
	for(int m=0; m<tbs[j].partners.size(); m++)
	    t.Transform( tbs[j].partners[m].pt );
	result.push_back(tbs[j]);
	}
    }
return true;
}

/* --------------------------------------------------------------- */
/* LookForAndTransform ------------------------------------------- */
/* --------------------------------------------------------------- */

// Read a file of synapse annotations. Find all the ones that
// map to the final image, transform them into final image
// coordinates, then append them to the vector result[].
//
static bool LookForAndTransform(
	const char			*name,
	char				*root,
	int					section,
	int					image_no,
	vector<image>		&images,
	bool				Warp,
	vector<uint8>		&tmap,
	vector<Bookmark>	&result )
{
uint32 w = images[image_no].w;
uint32 h = images[image_no].h;
char lname[2048];
sprintf(lname, "%s/%s", root, name);
FILE *fp = fopen( lname, "r" );
if( fp == NULL ) {
    printf("No file '%s'\n", lname);
    return false;
    }
printf("Found '%s', copying bookmarks.\n", lname);
fclose(fp);
if( strstr(name, "annotations-bookmarks") != NULL ) {
    vector<Bookmark> bookmarks;
    ReadBookmarks(lname, section, bookmarks);
    // Since each image is divided into patches (often 3 or more), and each patch
    // may be divided into sectors, there are many very similar transformations
    // mapping this image into global space.  However, there should be just one
    // point in global space that maps to this point in this image.  So try each
    // one, then look up the reverse transformation from global space, and find
    // the one that maps here.
    for(int j=0; j<bookmarks.size(); j++) {
	Point p = bookmarks[j].pt;
        int ix = int(p.x);
        int iy = int(p.y);
        // make sure it's legal
        if( ix < 0 ||  ix >= w || iy < 0 || iy >= h ) {
             printf("Synapse outside diagram?? %d %d %d %d\n", ix, w, iy, h);
             continue;
             }
	int patch = images[image_no].foldmap[ix + iy*w];
        // There can be small patches which were not used, so need to check for legality here
        //             // No longer needed since we compress patches and foldmaps on input
        if( patch == 0 )  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0)
            continue;
        TAffine	t;
	if( !Warp )  // just one transform per image
	    t = images[image_no].tf[patch];  // change to global coordinates
        else { //we are warping
	    int sector = tmap[ix+iy*w];
            t = images[image_no].sectors[patch][sector];
            }
        t.Transform( p );
        bookmarks[j].pt = p;
	result.push_back(bookmarks[j]);
	}
    }
return true;
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Main program for creating mosaics of grayscale, superpixel,
// and boundary maps from individual em images. Needs a file
// containing transforms and a mapping of image names to sections
// (layers).

int main( int argc, char **argv )
{

bool Debug = false;
bool Warp = false;      // generate warped images, that better align at seams?
bool FoldMasks = true;  // by default, use fold masks
double DontMoveStrength = 0.01;  // default value
bool Annotate = false;
bool make_tiles = false;
bool make_flat  = true;       // create a flat output image
bool make_map = true;         // create a map saying how image was generated?
bool matlab_order = false;    // render tiles in matlab order, not by closeness
bool RenumberSuperPixels = true; // by default, renumber these
string fold_dir   = ".";      // pre-pend this to fold mask locations
string region_dir = ".";      // where results go
string gray_dir   = ".";      // where gray scale images go
string sp_dir     = ".";      // where super-pixel maps go
string inv_dir    = ".";      // where information needed to invert transformations
                              // (map back to tiles) goes.
string bmap_dir   = ".";      // boundary maps go here
string rav_dir    = ".";      // Where to write the raveler tiles.
int scale = 1;                // scale everything down by a factor of 'scale'

vector<char *>noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    // process arguments here
    if( argv[i][0] != '-' )
	noa.push_back(argv[i]);
    else if( strcmp(argv[i],"-d") == 0 )
	Debug = true;
    else if( strncmp(argv[i],"-warp", 5) == 0 )
	Warp = true;
    else if( strcmp(argv[i],"-tiles") == 0 )
	make_tiles = true;
    else if( strcmp(argv[i],"-noflat") == 0 )
	make_flat = false;
    else if( strcmp(argv[i],"-nomap") == 0 )
	make_map = false;
    else if( strcmp(argv[i], "-nf") == 0 )
	FoldMasks = false;
    else if( strcmp(argv[i], "-a") == 0 )
	Annotate = true;
    else if( strcmp(argv[i], "-matlab") == 0 )
	matlab_order = true;
    else if( strcmp(argv[i], "-drn") == 0 )
	RenumberSuperPixels = false;
    else if( strncmp(argv[i], "-dms=",5) == 0 ) {
	DontMoveStrength = atof(argv[i]+5);
	}
    else if( strncmp(argv[i], "-fold_dir=",10) == 0 )
        fold_dir = string(argv[i]+10);
    else if( strncmp(argv[i], "-region_dir=",12) == 0 )
        region_dir = string(argv[i]+12);
    else if( strncmp(argv[i], "-gray_dir=",10) == 0 )
        gray_dir = string(argv[i]+10);
    else if( strncmp(argv[i], "-grey_dir=",10) == 0 )   // both spellings of grey, just in case
        gray_dir = string(argv[i]+10);
    else if( strncmp(argv[i], "-sp_dir=",8) == 0 )
        sp_dir = string(argv[i]+8);
    else if( strncmp(argv[i], "-inv_dir=",9) == 0 )
        inv_dir = string(argv[i]+9);
    else if( strncmp(argv[i], "-rav_dir=",9) == 0 )
        rav_dir = string(argv[i]+9);
    else if( strncmp(argv[i], "-bmap_dir=",10) == 0 )
        inv_dir = string(argv[i]+10);
    else if( strncmp(argv[i], "-s=",3) == 0 )
        scale = atoi(argv[i]+3);
    else
	printf("Ignored option '%s'\n", argv[i]);
    }
if( Warp )
    printf("Don't Move Strength = %f\n", DontMoveStrength);
if( noa.size() < 1 ) {
    printf("Usage: mos <simple file> [region xmin,ymin,dx,dy] [minlayer,maxlayer] [options]\n");
    exit( 42 );
    }

vector<char*>	dnames;   // directory names
vector<int>		lnums;    // parallel vector of layer numbers
int highest = 0; // highest and lowest layer numbers encountered
int lowest  = 1000000000;
vector<image> images;
map<string,int> imap;  // map of image names
uint32 w=0,h=0;        // size of each sub-image

int x0=0,y0=0;           // where to start output image
int xsize=-1,ysize=-1;  // size of output image, if specified
if( noa.size() >= 2 ) {
    if( sscanf(noa[1],"%d,%d,%d,%d", &x0, &y0, &xsize, &ysize) != 4 ) {
	printf("Expected x0,y0,xsize,ysize, got '%s'\n", noa[1]);
	exit( 42 );
	}
    }
int lspec1=-1, lspec2;   //layers specified by the user
if( noa.size() >= 3 ) {
    if( sscanf(noa[2],"%d,%d", &lspec1, &lspec2) != 2 ) {
	printf("Expected min,max layers, got '%s'\n", noa[2]);
	exit( 42 );
	}
    }
printf( "layer range = [%d, %d].\n", lspec1, lspec2 );
printf( "Fold mask   directory = '%s'\n",   fold_dir.c_str() );
printf( "Region      directory = '%s'\n", region_dir.c_str() );
printf( "Gray scale  directory = '%s'\n",   gray_dir.c_str() );
printf( "Super-pixel directory = '%s'\n",     sp_dir.c_str() );
printf( "Inverse-map directory = '%s'\n",    inv_dir.c_str() );
printf( "Raveler     directory = '%s'\n",    rav_dir.c_str() );

// read the file.  Here we assume all images are the same size, so if w and/or h are non-zero, we don't really need to read the
// images, since we already read at least one.

double	xmin =  BIG, ymin =  BIG;
double	xmax = -BIG, ymax = -BIG;
vector<double>	x1, y1, x2, y2;
vector<int>		z1, z2;
FILE			*fp = FileOpenOrDie( noa[0], "r" );
CLineScan		*ls = new CLineScan;

	for(;;) {

		if( ls->Get( fp ) <= 0 )
			break;

		if( !strncmp( ls->line, "DIR", 3 ) ) {

			// parallel vectors
			char *z = strtok( ls->line + 4, " " );
			char *n = strtok( NULL, " \n" );

			lnums.push_back( atoi( z ) );
			dnames.push_back( strdup( n ) );
		}
		else if( !strncmp( ls->line, "FOLDMAP", 7 ) ) {

			char	*tname = strtok( ls->line + 7, " '\n" );
			char	*mname = strtok( NULL, " '\n" );
			int		id_dum;

			if( !tname || !mname ) {
				printf( "Bad FOLDMAP statement '%s'.", ls->line );
				exit( 42 );
			}

			map<string,int>::iterator it = imap.find(string(tname));

			if( it != imap.end() )
				continue;  // already seen this one

			// Get image dimensions
			if( !w ) {
				uint8	*r = Raster8FromAny( tname, w, h );
				RasterFree( r );
				w /= scale;
				h /= scale;
			}

			image ii;

			ii.raster	= NULL;
			ii.foldmap	= NULL;
			ii.bmap		= NULL;
			ii.w		= w;
			ii.h		= h;
			ii.spbase	= 0;
			ii.rname	= strdup( tname );
			ii.fname	= strdup( mname );
			ii.spname	= NULL;
			ii.bname	= NULL;
			ii.spmap	= NULL;

			ZIDFromFMPath( ii.layer, id_dum, mname );

			if( lspec1 < 0 ||
				(lspec1 <= ii.layer && ii.layer <= lspec2) ) {

				// Point at drawing (d-suffixed) foldmap if exists

				char	*suf = strstr( mname, ".tif" );

				if( suf ) {

					char	buf[2048];

					sprintf( buf, "%.*sd.tif", suf - mname, mname );

					if( DskExists( buf ) ) {
						printf( "Using drawing file '%s'.\n", buf );
						free( ii.fname );
						ii.fname = strdup( buf );
					}
				}

				// Finally, add to collection
				imap[string(tname)] = images.size();
				images.push_back( ii );

				lowest  = min( lowest, ii.layer );
				highest = max( highest, ii.layer );
			}
		}
		else if( !strncmp( ls->line, "TRANSFORM", 9 ) ) {

			char	name[2048];
			double	a, b, c, d, e, f;
			int		nread;

			// newer format has '' delimited name
			// older just space separated

			if( ls->line[10] == '\'' ) {
				nread = sscanf(
				ls->line + 10, "'%[^']' %lf %lf %lf %lf %lf %lf",
				name, &a, &b, &c, &d, &e, &f );
			}
			else {
				nread = sscanf(
				ls->line + 10, "%s %lf %lf %lf %lf %lf %lf",
				name, &a, &b, &c, &d, &e, &f );
			}

			if( 7 != nread ) {
				printf( "Bad TRANSFORM statement '%s'.", ls->line );
				exit( 42 );
			}

			TAffine	tf( a, b, c/scale, d, e, f/scale );
			char	*fname = strtok( name, " ':" );

			//printf("File '%s'\n", fname);

			map<string,int>::iterator imit = imap.find(string(fname));

			if( imit == imap.end() ) {

				// This is normal with single layer image generation.

				//printf( "File in TRANSFORM has no FOLDMAP - ignored.\n" );
				continue;
			}

			int	k = imit->second;

			// get the rest of the string, now of the form ::123
			int	patch = atoi( strtok( NULL, " :'\n" ) );

			printf( "ImageIdx::patch = %d::%d\n", k, patch );

			// Make new slot?
			if( images[k].tf.size() <= patch ) {

				// initialize to an illegal transform
				images[k].tf.resize( patch+1, TAffine(0,0,0,0,0,0) );
				images[k].inv.resize( patch+1 );
			}

			images[k].tf[patch] = tf;
		}
		else if( !strncmp( ls->line, "SPMAP", 5 ) ) {

			char	name[2048], where[2048];
			int		z, nprm;

			nprm = sscanf( ls->line + 5, "%s %s %d", name, where, &z );

			if( nprm < 2 ) {
				printf( "Bad SPMAP statement '%s'.", ls->line );
				exit( 42 );
			}

			char	*fname = strtok( name, " ':" );

			//printf("File '%s'\n", fname);

			map<string,int>::iterator imit = imap.find(string(fname));

			if( imit == imap.end() ) {

				// This is normal with single layer image generation.
				// If it's a layer we care about, print a message.
				// Otherwise silently ignore.

				int	layer;

				if( nprm > 2 )
					layer = z;
				else
					layer = FileNameToLayerNumber( dnames, lnums, fname );

				if( lspec1 < 0 ||
					(lspec1 <= layer && layer <= lspec2) ) {

					printf(
					"File in SPMAP has no image - ignored.\n" );
				}

				continue;
			}

			int	k = imit->second;

			printf( "ImageIdx = %d\n", k );
			images[k].spname = strdup( where );
		}
		else if( !strncmp( ls->line, "BOUNDARYMAP", 11 ) ) {

			char	name[2048], where[2048];
			int		z, nprm;

			nprm = sscanf( ls->line + 11, "%s %s %d", name, where, &z );

			if( nprm < 2 ) {
				printf( "Bad BOUNDARYMAP statement '%s'.", ls->line );
				exit( 42 );
			}

			char	*fname = strtok( name, " ':" );

			//printf("File '%s'\n", fname);

			map<string,int>::iterator imit = imap.find(string(fname));

			if( imit == imap.end() ) {

				// This is normal with single layer image generation.
				// If it's a layer we care about, print a message.
				// Otherwise silently ignore.

				int	layer;

				if( nprm > 2 )
					layer = z;
				else
					layer = FileNameToLayerNumber( dnames, lnums, fname );

				if( lspec1 < 0 ||
					(lspec1 <= layer && layer <= lspec2) ) {

					printf(
					"File in BOUNDARYMAP has no image - ignored.\n" );
				}

				continue;
			}

			int	k = imit->second;

			printf( "ImageIdx = %d\n", k );
			images[k].bname = strdup( where );
		}
		else if( !strncmp( ls->line, "MPOINTS", 7 ) ) {

			double	a, b, c, d;
			int		za, zc;

			if( 6 != sscanf( ls->line + 7,
				"%d %lf %lf %d %lf %lf",
				&za, &a, &b, &zc, &c, &d ) ) {

				printf( "Bad MPOINTS statement '%s'.", ls->line );
				exit( 42 );
			}

			if( lspec1 < 0 ||
				(lspec1-1 <= za && za <= lspec2+1) ||
				(lspec1-1 <= zc && zc <= lspec2+1) ) {

				z1.push_back( za );
				x1.push_back( a/scale );
				y1.push_back( b/scale );

				z2.push_back( zc );
				x2.push_back( c/scale );
				y2.push_back( d/scale );
			}
		}
		else if( !strncmp( ls->line, "BBOX", 4 ) ) {

			if( 4 != sscanf( ls->line + 4,
				"%lf %lf %lf %lf", &xmin, &ymin, &xmax, &ymax ) ) {

				printf( "Bad BBOX statement '%s'.\n", ls->line );
				exit( 42 );
			}

			xmin /= scale;
			ymin /= scale;
			xmax /= scale;
			ymax /= scale;
		}
		else if( !strncmp( ls->line, "IMAGESIZE", 9 ) ) {

			if( 2 != sscanf( ls->line + 10, "%d %d", &w, &h ) ) {
				printf( "Bad IMAGESIZE line '%s'.\n", ls->line );
				exit( 42 );
			}

			w /= scale;
			h /= scale;
		}
		else {
			printf( "Unknown line '%s'.", ls->line );
			exit( 42 );
		}
	}

	delete ls;
	fclose( fp );

//VMStats( stdout );
//exit( 11 );

	if( !images.size() ) {
		printf( "No images in input.\n" );
		exit( 42 );
	}



	GlobalBounds( xmin, xmax, ymin, ymax, images, w, h );

	OffsetCoords( xmin, xmax, ymin, ymax, images,
		x1, y1, x2, y2, x0, y0, xsize, ysize );

















// Now we will need to create one map, and one image, per layer.
if( lspec1 >= 0 ) {  // layer numbers were specified
    lowest = max(lowest, lspec1);
    highest= min(highest, lspec2);
    printf("Layers reset to %d - %d\n", lowest, highest);
    }

// Now find the biggest d we can use such that a grid of points d x d in the image space will hit every
// integer coordinate in global space.  We want d as big as possible to avoid un-needed work,
// but don't want to miss any pixels in global space.
double delta_image_space = 2.0;  // implausibly high value
for(int i=0; i<images.size(); i++) {
    for(int k=1; k<images[i].tf.size(); k++) {
	if( images[i].tf[k].det() == 0.0 )
	    continue;
        Point p[4];
        p[0] = Point(0.0,0.0);  // compose a unit square in image space
        p[1] = Point(1.0,0.0);
        p[2] = Point(0.0,1.0);
        p[3] = Point(1.0,1.0);
        for(int j=0; j<4; j++)
           images[i].tf[k].Transform( p[j] );
        // How big of a square might fit in global space?
        double xmin = min(p[0].x,min(p[1].x,min(p[2].x,p[3].x)));
        double ymin = min(p[0].y,min(p[1].y,min(p[2].y,p[3].y)));
        double xmax = max(p[0].x,max(p[1].x,max(p[2].x,p[3].x)));
        double ymax = max(p[0].y,max(p[1].y,max(p[2].y,p[3].y)));
        double dx = xmax-xmin;
        double dy = ymax-ymin;
        // the allowable scale is 1/(bigger of these)
        double d = 1.0/max(dx,dy);
        printf("d for image %d, patch %d is %f\n", i, k, d);
        // add a tiny safety factor for floating point problems
        delta_image_space = min(delta_image_space, d*0.999999);
	}
    }
printf("Image space delta is %f\n", delta_image_space);

// Since some (many) tiles will not have super-pixel maps,
// create a single blank map they all can point to.
// Likewise for boundary maps.
uint16 *BlankSPMap = (uint16*)malloc(w*h*sizeof(uint16));
memset( BlankSPMap, 0, w * h * sizeof(uint16) );
uint8* BlankBMap = (uint8*)BlankSPMap;

// The following huge loop goes through the layers one at a time, producing 6 outputs
//   (a) a 'before.%d.png', which is the results of simply mapping the input affine
//   (b) an 'after.%d.png', which is the result after edge merging
//   (c) an 'spmap.%d.png', which is the super-pixel map
//   (d) a 'map.%d.png', which describes how each pixel got there
//   (e) a 'mapping.%d.txt', which has inverse transforms for each global point, back
//       to the original image
//   (f) a 'map-sp.txt', which describes how the SP maps of each tile map to the SP of the combined image
//       -1 = SP value was not in the original image
//       -2 = SP was not used in the final output
for(int out_layer = lowest; out_layer <= highest; out_layer++) {  //keep going until we find nothing...

    printf("starting to generate layer %d...\n", out_layer);
    // start by tossing out previous images, if any
    for(int i=0; i<images.size(); i++) {
	if( images[i].raster != NULL )
	    RasterFree( images[i].raster );
	if( images[i].foldmap != NULL )
	    RasterFree( images[i].foldmap );
	if( images[i].spmap != NULL && images[i].spmap != BlankSPMap ) {
	     free(images[i].spmap);
         images[i].spmap = NULL;   // don't keep a pointer to a free'd item.
	}
	if( images[i].bmap != NULL && images[i].bmap != BlankBMap ) {
	    free(images[i].bmap);      // and previous boundary maps
        images[i].bmap = NULL;
    }
	}
    // Find the relevant images, and read their data.
    vector<int> relevant_images;
    // See if we find any boundary maps
    bool AnyBMap = false;
    // Since each superpixel map is numbered independently per image,
    // map all of these to non-overlapping ranges.
    int MaxSPUsed = 0;
    for(int i=0; i<images.size(); i++) {      // for each picture
        if( images[i].layer != out_layer )     // if it's not the layer we're looking for
	    continue;
        if( RuleOutIntersection(w, h, images[i].tf,  0,0, ROUND(xmax), ROUND(ymax)) ) {
            printf("Image %d ruled out\n", i);
	    continue;
            }
        relevant_images.push_back(i);
	uint32 ww, hh;  // if 'scale' option is used, this is the original size
	if( images[i].raster == NULL ) {
	    images[i].raster  = LoadNormImg( images[i].rname, ww, hh );
            ImageResize(images[i].raster, ww, hh, scale);
            }
	if( images[i].foldmap == NULL ) {
	    if( FoldMasks ) {
                // if the fold mask name is not rooted, pre-pend the fold directory name
                string file_name;
                if( images[i].fname[0] == '/' )
		    file_name = string(images[i].fname);
	        else
                    file_name = fold_dir + "/" +string(images[i].fname);  //pre-pend the fold directory name
		images[i].foldmap = Raster8FromAny( file_name.c_str(), ww, hh );
                ImageResize(images[i].foldmap, ww, hh, scale);
		FillInHolesInFoldmap(images[i].foldmap, w, h);
		FoldmapRenumber(images[i]);
		}
	    else {  // no foldmasks wanted; just create a field of 1s.
		images[i].foldmap = (uint8*)malloc(w*h*sizeof(uint8));
		for(int k=0; k<w*h; k++)
		    images[i].foldmap[k] = 1;
		}
	    }
        // Now read the super-pixel maps, if they exist, otherwise point to blanks
	if( images[i].spmap == NULL ) {
            uint16 *test = NULL;
            if( images[i].spname == NULL ) {
		printf("No SPmap defined for image '%s'\n", images[i].rname);
		test = BlankSPMap;
		}
	    else {
		uint32	wdum, hdum;
		test = Raster16FromPng( images[i].spname, wdum, hdum );
	        if( test == NULL ) {
		    printf("Cannot read superpixel map '%s'.\n", images[i].spname);
	            test = BlankSPMap;
		    }
	        else {
                    ImageResize16bits(test, ww, hh, scale);
                    if( RenumberSuperPixels ) {
                        images[i].spbase = MaxSPUsed;
		        RemapSuperPixelsOneImage(test, w, h, MaxSPUsed, images[i].SPmapping);
			}
                    printf("MaxSPUsed is now %d\n", MaxSPUsed);
		    }
		}
	    images[i].spmap = test;
	    }
        // Now read the boundary maps if they exist
	if( images[i].bmap == NULL ) {
            uint8 *test = NULL;
            if( images[i].bname == NULL ) {
		printf("No boundary map defined for image '%s'\n", images[i].rname);
		test = BlankBMap;
		}
	    else { // name is defined; see if file exists
                if( strstr(images[i].bname,".png") != NULL ) {
                    uint32 wj, hj; // not used
		    test = Raster8FromPng(images[i].bname, wj, hj);
                    }
                else
				test = LoadNormImg( images[i].bname, ww, hh );
	        if( test == NULL ) {
		    printf("Cannot read boundary map file '%s'.\n", images[i].bname);
	            test = BlankBMap;
		    }
                else {
		    ImageResize(test, ww, hh, scale);
		    AnyBMap = true;  // at least one bmap was found for this layer
		    }
		}
	    images[i].bmap = test;
	    }
        }





    // First calculate a map, saying for each pixel in global space which image it comes from
    // Likewise, PatchFrom tells which patch of the image this was from
    // Do this by starting at the center of each image, working out, and setting the
    // global pixel if it is not already set.
    // Note that we ignore xmin and ymin, and always start at 0,0;
    uint32 nx = ROUND(xmax)+1;  // number of pixels in X
    uint32 ny = ROUND(ymax)+1;
    printf("Image size is %d by %d, %d images\n", nx, ny, images.size());
    vector<uint16>imap(nx*ny,0);
    if( relevant_images.size() > 0xFFFF ) {  // since we hold this as 16 bits
	printf("Too many images (%d) on layer %d\n", relevant_images.size(), out_layer);
	exit( 42 );
	}
    vector<uint8>PatchFrom(nx*ny,0);

    // Now the choice of fill-in methods
    if( matlab_order ) {
	FillInMatlabOrder(imap, PatchFrom, nx, ny,
	 w, h, delta_image_space, images, relevant_images);
        }
    else {
	FillInFromCenterOut(imap, PatchFrom, nx, ny,
	 w, h, delta_image_space, images, relevant_images);
	}

    // Now compress the map/patch array.  We assume less than 65536 combinations are used, so we can express
    // this as a 16 bit TIF.  The combinations start at 1, and run sequentially.
    map<int,int> temp;
    vector<Triple> Triples(1);  // will fill up as they are found
    vector<int>histo(10);       // small histogram, for testing
    for(uint32 i=0; i<imap.size(); i++) {
        if( imap[i] == 0 )
	    continue;	   // was never set - leave as is
	int id = imap[i] + (PatchFrom[i] << 16);  // create an int.
        map<int,int>::iterator it = temp.find(id);
        if( it == temp.end() ) { // add it
	    Triple t(imap[i], PatchFrom[i], 0);
	    temp.insert(pair<int,int>(id, Triples.size()));
            printf("New combo %5d: %d %d at %d\n", Triples.size(), t.image, t.patch, i);
            Triples.push_back(t);
	    it = temp.find(id);   // now should find it, since we just added it
	    }
        imap[i] = it->second;
        if( it->second < 10 )
	    histo[it->second]++;
	}
    for(int i=0; i<10; i++)
	printf("Entry %d used %d times.\n", i, histo[i]);






    // write the map; only used for debugging for now
    char fname[256];
    if( Debug ) {
        printf("Try writing 16 bit map\n");
        Raster16ToPng16( "test.png", &imap[0], nx, ny );
        sprintf(fname,"map.%d.png", out_layer);
        printf("Writing  8 bit map in png format to %s\n", fname);
        vector<uint8> copy(nx*ny, 0);  // make an 8 bit copy
        for(int i=0; i<nx*ny; i++)
	    copy[i] = imap[i];
        Raster8ToPng8( fname, &copy[0], nx, ny );
        sprintf( fname,"pf.%d.tif", out_layer );
        //Raster8ToTif8( fname, &PatchFrom[0], nx, ny );
        }






    printf("Starting 'before' image\n");
    // Now, create a 'before' picture with seams
    vector<uint8>before(nx*ny,0);
    for(int y=0; y < ny; y++) {
        if( (y & 0x3FF) == 0 ) {
	    printf("."); fflush(stdout);
            }
	for(int x=0; x<nx; x++) {
            uint32 bi= uint32(x) + uint32(y)*nx;  // big index
            uint16 indx = imap[bi];
	    uint16 from = Triples[indx].image;  // what image does this pixel come from?
            uint8 patch = Triples[indx].patch;
	    if( from == 0 )
		continue;  // no image sets this pixel
	    int img = relevant_images[from-1];
	    Point p(x, y);
	    images[img].inv[patch].Transform( p );
	    // of course, this should be in the image, but let's double check
	    if( p.x >= 0.0 && p.x < w-1 && p.y >= 0.0 && p.y < h-1 ) { // then we can interpolate
		int ix = int(p.x);
		int iy = int(p.y);
		double alpha = p.x-ix;
		double beta  = p.y-iy;
		int nn = ix + w*iy;  // index into original image
		uint8* pic = images[img].raster;
		double pix =
		     (1-alpha)*(1-beta) * pic[nn] +
			alpha *(1-beta) * pic[nn+1] +
		     (1-alpha)*   beta  * pic[nn+w] +
			alpha *   beta  * pic[nn+w+1];
		before[bi] = ROUND(pix);
		}
	    } // to aid debugging, draw lines at the image boundaries in the before image.
         }
    printf("Done creating before image; draw lines next\n");
    if( Annotate ) {
	vector<Point> edges;
	for(int x=0; x<w; x++) {
	    edges.push_back(Point(x, 0.0));
	    edges.push_back(Point(x, h-1));
	    }
	for(int y=0; y<h; y++) {
	    edges.push_back(Point(0.0, y));
	    edges.push_back(Point(w-1, y));
	    }
	printf("%d edge pixels\n", edges.size());
	// OK, now transform this edge list by each of the transforms, then color it in.
	for(int i=0; i<images.size(); i++) {
	    if( images[i].layer != out_layer ) // only want edges on current layer
		continue;
	    for(int k=1; k<images[i].tf.size(); k++) {
		if( images[i].tf[k].det() == 0.0 )
		    continue;
		for(int j=0; j<edges.size(); j++) {
		    Point p = edges[j];
		    images[i].tf[k].Transform( p );
		    int ix = ROUND(p.x);
		    int iy = ROUND(p.y);
		    if( 0 <= ix && ix < nx && 0 <= iy && iy < ny )
			before[uint32(ix) + nx*uint32(iy)] = 255;
		    }
		}
        // draw the Tile ID in the center of the tile.
        Point	p = Point( w/2, h/2 );
		images[i].tf[1].Transform( p );
        const char *fn = FileNamePtr( images[i].rname );
        printf( "Draw image %s at %f %f\n", fn, p.x, p.y );
        if( p.x > 0 && p.x < nx-1 && p.y > 0 && p.y < ny-1 )
			DrawText(&before[0], nx, ny, fn, int(p.x), int(p.y));
	    }
	//also to aid in debugging, draw circles at correspondence points
	for(int i=0; i<x1.size(); i++) {
	    if( z1[i] == out_layer || z2[i] == out_layer ) {
		DrawCircle( &before[0], nx, ny, x1[i], y1[i], 10.0 );
		DrawCircle( &before[0], nx, ny, x2[i], y2[i], 10.0 );
		}
	    }
        }
    //sprintf( fname,"before.%05d.tif", out_layer );
    //Raster8ToTif8( fname, &before[0], nx, ny );
    if( make_flat ) {
	if( nx*ny > 0x7FFFFFFFu ) {  // is it *really* big?
	    printf("write half\n");
	    sprintf(fname,"before.%05d.a.png", out_layer);
	    Raster8ToPng8( fname, &before[0], nx, ny/2 );
	    printf("write the other half\n");
	    sprintf(fname,"before.%05d.b.png", out_layer);
	    Raster8ToPng8( fname, &before[nx*ny/2], nx, ny/2 );
	    printf("both written\n");
	    }
	else {   // this is the usual case
	    sprintf(fname,"before.%05d.png", out_layer);
	    Raster8ToPng8( fname, &before[0], nx, ny );
	    }
	}





    // if only simple output is needed, we are done with this layer.  Write map and mapping file
    if( !Warp ) {
        // write the map file, if requested.
        if( make_map ) {
            sprintf(fname,"%s/%s/map.%05d.png", region_dir.c_str(), inv_dir.c_str(), out_layer);
            printf("write Pixel mapping file %s\n", fname);
            Raster16ToPng16( fname, &imap[0], nx, ny );
	    // ------------------------------------------------ write the mapping text file  --------------------
	    // This is the text file that contains the transformations that describe how each pixel got there.
	    sprintf(fname,"%s/%s/mapping.%05d.txt", region_dir.c_str(), inv_dir.c_str(), out_layer);
	    FILE *ftxt = FileOpenOrDie( fname, "w" );
	    //write the image names, and their indexes
	    for(int k=0; k<relevant_images.size(); k++) {      // for each picture
		int i = relevant_images[k];
		fprintf(ftxt,"IMAGE %d '%s'\n", k, images[i].rname);
		}
	    for(int k=1; k<Triples.size(); k++) {
		int from = Triples[k].image;          // index into the relevant pictures array
		int image = relevant_images[from-1];  // actual image number
		int patch = Triples[k].patch;
		//int sector = Triples[k].sector;     // not used since no sectors yet
		TAffine t = images[image].inv[patch];  // transform from global to image
		fprintf(ftxt, "TRANS %d %d ", k, from-1);
		t.TPrint( ftxt );
		}
	    fclose(ftxt);
	    }
        else
            printf("User requested no map file.\n");
	continue;
        }




    // Break each image into triangles
    const int N=16;  // number of triangles per image

    // create the control points.  Triangles numbered CCW from angle 0.  Vertex N is at the center.
    // There is just one list that applies to all images.
    vector<Point> vtx(N+1);
    double ctrx = double(w-1)/2.0;
    double ctry = double(h-1)/2.0;
    vtx[N] = Point(ctrx, ctry);
    double tiny = 1e-5; // move vertices slightly so none exactly on the edge of a region
    for(int i=0; i<N; i++) {
        double theta = 2*PI*(i)/N;
        double x,y;
        if( N/8 <= i && i <= 3*N/8 ) { // top side
            x = ctrx + ctrx*cot(theta);
            y = h-1;
	    }
        else if( 3*N/8 <= i && i <= 5*N/8 ) { // left side
            x = 0;
            y = ctry - ctry*tan(theta);
            }
       else if( 5*N/8 <= i && i <= 7*N/8 ) { // bottom
            x = ctrx - ctrx*cot(theta);
            y = 0;
	    }
       else { // right side
            x = w-1;
            y = ctry + ctry*tan(theta);
            }
       // if to close to the very edge, move outside by a tiny amount
       if( abs(x)       < 0.1 ) x -= tiny;
       if( abs(x-(w-1)) < 0.1 ) x += tiny;
       if( abs(y)       < 0.1 ) y -= tiny;
       if( abs(y-(h-1)) < 0.1 ) y += tiny;
       vtx[i].x = x;
       vtx[i].y = y;
       }
    for(int i=0; i<vtx.size(); i++)
        printf("vtx[%d] = (%f,%f)\n", i, vtx[i].x, vtx[i].y);

    // now, create a map that maps each pixel of an image into the relevant triangle.
    vector<uint8> tmap(w*h,0xFF);
    for(int i=0; i<N; i++) {
        Point va = vtx[N];
        Point vb = vtx[i];
        Point vc = vtx[(i+1)%N];
        // triangle goes va,vb,vc,va in a ccw direction
        double xmin = min(va.x, min(vb.x, vc.x));
        double xmax = max(va.x, max(vb.x, vc.x));
        double ymin = min(va.y, min(vb.y, vc.y));
        double ymax = max(va.y, max(vb.y, vc.y));
        for(int x = int(max(xmin,0.0)); x < w && x <= xmax+1.0; x++) {
	    for(int y = int(max(ymin,0.0)); y < h && y <= ymax+1.0; y++) {
                Point p(x,y);
	        if( LeftSide(va,vb,p) && LeftSide(vb,vc,p) && LeftSide(vc,va,p) )
		    tmap[x+w*y] = i;
		}
	    }
	}
    // Verify all bits are set.
    bool bogus = false;
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
	    if(tmap[x+w*y] == 0xFF) {
		printf("Bogus map! x=%d, y=%d\n", x, y);
		bogus = true;
		}
	    }
	}
    if( bogus )
	exit( 42 );
    //Raster8ToTif8( "tmap.tif", &tmap[0], w, h );

    // Now, for each patch of each relevant image create N triangles using N+1 pts each, all transferred
    // to the global reference frame.
    vector<triangle> tris;
    vector<Point> gvtx;
    for(int m=0; m<relevant_images.size(); m++) {
        int i = relevant_images[m];
	images[i].FirstGlobalPoint = gvtx.size();
        images[i].FirstTriangle = tris.size();
	for(int j=1; j<images[i].tf.size(); j++) {  // patches are from 1 on up
            // copy the points
            int n0 = gvtx.size();
	    for(int k=0; k<N+1; k++) {
	        Point p = vtx[k];
                images[i].tf[j].Transform( p );
		gvtx.push_back(p);
	        }
            for(int k=0; k<N; k++) {
		triangle t;
		t.v[0] = n0+N; // the center
		t.v[1] = n0+k;
                t.v[2] = n0+((k+1)%N);
                tris.push_back(t);
		}
	    }
	}
    // for each triangle, create a matrix that transforms to barycentric coordinates
    for(int k=0; k<tris.size(); k++) {
        Point v0 = gvtx[tris[k].v[0]], v1 = gvtx[tris[k].v[1]], v2 = gvtx[tris[k].v[2]];
        printf("Tri: (%f %f) (%f %f) (%f %f)\n", v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
        double a[3][3];
        a[0][0] = v0.x; a[0][1] = v1.x; a[0][2] = v2.x;
        a[1][0] = v0.y; a[1][1] = v1.y; a[1][2] = v2.y;
        a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
        Invert3x3Matrix( tris[k].a, a );
        }

    // Now write the triangles out for debugging
    if( Debug )
        WriteTriangles("tris", tris, gvtx);

    // Then, look for places where adjacent pixels come from different pictures.
    // For each such spot (or a reasonable subset, find the delta needed make them line up.
    vector<glob_spot> spots;
    const int delta = 40;
    for(int ty=0; ty <ny; ty += delta) {
	printf("ty=%d\n", ty);
        int Ntrans = 0;   // count number of transitions
        int prev_x = -2*delta;  // this line and the next keep track of the last transition
        int prev_i1, prev_i2;   // Needed to avoid problems with nearest-point interleaving
	for(int tx = 0; tx+1 < nx; tx++) {
	    uint32 ind = uint32(ty)*nx + uint32(tx);
            int indx1 = imap[ind];
            int indx2 = imap[ind+1];
            int image1 = Triples[indx1].image;
            int image2 = Triples[indx2].image;
	    if( image1 != image2 && image1 > 0 && image2 > 0 ) { // found a boundary
                if (tx - prev_x < delta && // see if it's a duplicate of previous one
                 (indx1 == prev_i1 && indx2 == prev_i2 || indx1 == prev_i2 && indx2 == prev_i1))
		    printf("Same images, too close, at %d\n", tx);
                else { // it's Ok
		    Ntrans++;
		    int i1 = relevant_images[image1-1];
		    int i2 = relevant_images[image2-1];  // point is on boundary of image i1 and i2
		    //i2 = i1;  // temp, just to make sure we get 0,0
		    Point p1(tx,ty);
		    CommonPoint( spots, images, p1,i1, Triples[indx1].patch, out_layer );
		    prev_x = tx; prev_i1 = indx1; prev_i2 = indx2;
		    }
		}
	    }
        if( Ntrans > 100 ) {
	    printf("Too many transitions (%d) in map!\n", Ntrans);
            if( !Debug ) {
		printf("Writing map as 'test.png'\n");
		Raster16ToPng16( "test.png", &imap[0], nx, ny );
		}
            exit( 42 );
	    }
	}
    for(int tx=0; tx <nx; tx += delta) {
	printf("tx=%d\n", tx);
        int Ntrans = 0;   // count number of transitions
        int prev_y = -2*delta;  // this line and the next keep track of the last transition
        int prev_i1, prev_i2;   // Needed to avoid problems with nearest-point interleaving
	for(int ty = 0; ty+1 < ny; ty++) {
	    uint32 ind = uint32(ty)*nx + uint32(tx);  // note that ind+nx is up by one in y
            int indx1 = imap[ind];
            int indx2 = imap[ind+nx];
            int image1 = Triples[indx1].image;
            int image2 = Triples[indx2].image;
	    if( image1 != image2 && image1 > 0 && image2 > 0 ) { // found a boundary
                if (ty - prev_y < delta && // see if it's a duplicate of previous one
                 (indx1 == prev_i1 && indx2 == prev_i2 || indx1 == prev_i2 && indx2 == prev_i1))
		    printf("Same images, too close, at %d\n", ty);
                else { // it's Ok
		    Ntrans++;
		    int i1 = relevant_images[image1-1];  // -1 since 0 is reserved for 'no image maps here'
		    int i2 = relevant_images[image2-1];  // point is on boundary of image i1 and i2
		    //i2 = i1;  // temp, just to make sure we get 0,0
		    Point p1(tx,ty);
		    CommonPoint( spots, images, p1,i1, Triples[indx1].patch, out_layer );
		    prev_y = ty; prev_i1 = indx1; prev_i2 = indx2;
		    }
		}
	    }
        if( Ntrans > 100 ) {
	    printf("Too many transitions (%d) in map!\n", Ntrans);
            if( !Debug ) {
		printf("Writing map as 'test.png'\n");
		Raster16ToPng16( "test.png", &imap[0], nx, ny );
		}
            exit( 42 );
	    }
	}

    printf("%d global control points; %d constraints\n", gvtx.size(), spots.size() );
    int nvs = gvtx.size();  // number of control points

    // Now generate the constraints.  Declare the (sparse) normal matrix and the RHS
    // This whole section needs to be changed to use patches.
    vector<LHSCol> norm_mat(2*nvs);  // 2 coordinates for each spot
    vector<double> RHS(2*nvs, 0.0); //
    // first, add the constraints that each should be near where they started...
    for(int i=0; i<nvs; i++) {
	int j = 2*i;
        // the central vertex should have a strong constraint.  The others should be more
        // free to move.  The best strength is not clear; stronger will make the results
        // closer to the global alignment; weaker will allow better connections between the
        // tiles.
        double coeff = DontMoveStrength * ((i % (N+1)) == N ? 100 : 1);
        int vx[1] = {j};  // variable index
	AddConstraint(norm_mat, RHS, 1, vx, &coeff, gvtx[i].x*coeff);   // 1.0 * var = global spot
        vx[0] = j+1;
	AddConstraint(norm_mat, RHS, 1, vx, &coeff, gvtx[i].y*coeff);   // then multiply by coeff;
	}
    // Now add the constraints from the common points
    for(int i=0; i<spots.size(); i++) {
	printf("\nSpot %d\n", i);
        // add the spot to spot constraints.  Constraints are pairwise,
        // and up to 4 images may overlap, so create all pairs
        for(int j=0; j<spots[i].which.size()-1; j++) {
	    int im1 = spots[i].which[j];
	    int patch1 = spots[i].patch[j];
	    Point p1 = spots[i].where[j];
            printf("j=%d\n",j);
            int indices[6];   // this will form the constraints.  Variables with these indices
            double fracs[6];  // get these fractions
            FindTriangle(images, im1, patch1, p1, tris, gvtx, N, indices, fracs);
	    for(int k=j+1; k<spots[i].which.size(); k++) {
		// OK, now make a constraint between jth and kth entry. Change weight?
		int im2 =    spots[i].which[k];
                int patch2 = spots[i].patch[k];
                Point p2 =   spots[i].where[k];
                double lfrac[3];  // local fractions
                FindTriangle(images, im2, patch2, p2, tris, gvtx, N, indices+3, lfrac);
                fracs[3] = -lfrac[0]; fracs[4] = -lfrac[1]; fracs[5] = -lfrac[2];
                // we have vertex numbers, but we need variable numbers (2 variables per pt)
                int inds[6];
                for(int m=0; m<6; m++)
		    inds[m] = indices[m]*2;
	        AddConstraint(norm_mat, RHS, 6, inds, fracs, 0.0);  // for X
                for(int m=0; m<6; m++)
		    inds[m] = indices[m]*2 + 1;
	        AddConstraint(norm_mat, RHS, 6, inds, fracs, 0.0);  // for Y
		}
	    }
	}
    // Solve the system
    vector<double> X(2*nvs);
    fflush(stdout);
    WriteSolveRead( X, norm_mat, RHS, "A-mos", 1, !Debug );
    PrintMagnitude( X );
    fflush(stdout);

    vector<Point> NewG(nvs);
    for(int i=0; i<nvs; i++)
	NewG[i] = Point(X[2*i], X[2*i+1]);
    if( Debug )
        WriteTriangles("tris2", tris, NewG);

    // how far did the points move?
    double MoveThresh = 110.0; // should be parameter
    int NTooFar = 0;
    for(int i=0; i<nvs; i++) {
        double dx = gvtx[i].x - NewG[i].x;
        double dy = gvtx[i].y - NewG[i].y;
        double move = sqrt(dx*dx+dy*dy);
	printf("Image %d, pt %d, was (%f,%f), now (%f,%f), moved %f %s\n",
	 i/(N+1), i%(N+1), gvtx[i].x, gvtx[i].y, NewG[i].x, NewG[i].y, move, move>MoveThresh ? "<-----" : "");
        if( move > MoveThresh ) {
	    NewG[i] = gvtx[i];
	    NTooFar++;
	    }
	}
    // Now, for every patch of every image, compute a vector of N transforms, one corresponding
    // to each original triangle.
    for(int m=0; m<relevant_images.size(); m++) {
        int i = relevant_images[m];
        vector<TAffine> empty;
        images[i].sectors.push_back(empty);  // since there is no patch 0
        images[i].sinvs  .push_back(empty);
	for(int j=1; j<images[i].tf.size(); j++) {
	    vector<TAffine> ltfs(N);  // an N element vector of transforms
	    vector<TAffine> invs(N);  // and the inverses to these
            printf("Image %d, patch %d\n", i, j);
            images[i].tf[j].TPrint();
	    for(int k=0; k<N; k++) {
                // original points (in local coords, vtx,  are):
                int i0 = N;
                int i1 = k;
                int i2 = (k+1)%N;
		// Now, find a transformation that maps ORIG into the new cpts
		// first, create a transform that maps a unit right triangle to the original pts
                TAffine o(vtx[i1].x-vtx[i0].x, vtx[i2].x-vtx[i0].x, vtx[i0].x,
                 vtx[i1].y-vtx[i0].y, vtx[i2].y-vtx[i0].y, vtx[i0].y);

		// and the final points, in global space are
		int n0 = images[i].FirstGlobalPoint + (j-1)*(N+1);
                i0 += n0; i1 += n0; i2 += n0;
                // now one that maps the a unit right triangle to the desired final points
                TAffine c(NewG[i1].x-NewG[i0].x, NewG[i2].x-NewG[i0].x, NewG[i0].x,
                 NewG[i1].y-NewG[i0].y, NewG[i2].y-NewG[i0].y, NewG[i0].y);
                 // now, to get from the original to the final, apply o^-1, then c;
                TAffine oi, t;
                oi.InverseOf( o );
                t = c * oi;
                t.TPrint();
                ltfs[k] = t;
                invs[k].InverseOf( t );
		}
            images[i].sectors.push_back( ltfs );
            images[i].sinvs  .push_back( invs );
	    }
	}

    // Find the improvement, if any.
    MeanStd pre, post;
    for(int i=0; i<spots.size(); i++) {
	// should do all pairs, but at least do two.  First, use original affine map
	Point p0 = spots[i].where[0];
        images[spots[i].which[0]].tf[spots[i].patch[0]].Transform( p0 );
        Point p1 = spots[i].where[1];
        images[spots[i].which[1]].tf[spots[i].patch[1]].Transform( p1 );
        printf(" (%f %f) (%f %f) %f\n", p0.x, p0.y, p1.x, p1.y, p0.Dist( p1 ) );
        pre.Element( p0.Dist( p1 ) );

        //Now, re-do with new sector maps
        int ix = int(spots[i].where[0].x);
        int iy = int(spots[i].where[0].y);
        int sect0 = tmap[ix + w*iy];
        ix = int(spots[i].where[1].x);
        iy = int(spots[i].where[1].y);
        int sect1 = tmap[ix + w*iy];
	p0 = spots[i].where[0];
        images[spots[i].which[0]].sectors[spots[i].patch[0]][sect0].Transform( p0 );
        p1 = spots[i].where[1];
        images[spots[i].which[1]].sectors[spots[i].patch[1]][sect1].Transform( p1 );
        printf("   (%f %f) (%f %f) %f\n", p0.x, p0.y, p1.x, p1.y, p0.Dist( p1 ) );
        post.Element( p0.Dist( p1 ) );
        }
    double mean, std, mean2, std2;
    pre. Stats(mean,  std);
    post.Stats(mean2, std2);
    printf("Average error was %f, is now %f; %d points not moved of %d total\n", mean, mean2, NTooFar, nvs);

    // recalculate the map, using sectors.  Add a factor of 0.9 to the image delta since the scale may have
    // changed a little.  Could/Should re-calculate from scratch, but since we only allow minor changes
    // this should be OK.
    printf("Create the new map\n");
    memset( &imap[0], 0, nx * ny * sizeof(uint16) );
    vector<uint8>SectorFrom(nx*ny,0);
    double report = 100;
    // find the center of each picture.  May not be integer
    Point center((w-1)/2.0, (h-1)/2.0);
    for(double r=0.0; r <= max(center.x,center.y)+0.0001; r += delta_image_space*0.9) {
        if( r >= report ) {
            printf("Starting radius %f\n", r);
            report += 100;
            }
	vector<Point> pts;
	PointsInRing(pts, center, r, delta_image_space*0.9, w, h);
        // Now transform this by all relevant images
	for(int k=0; k<relevant_images.size(); k++) {      // for each picture
            int i = relevant_images[k];

	    for(int j=0; j<pts.size(); j++) {
		Point p(pts[j]);
		int ix = int(p.x);  // PointsInRing only returns legal points
		int iy = int(p.y);
                //printf("j=%d ix=%d iy=%d\n", j, ix, iy);
                int patch = images[i].foldmap[ix + w*iy];
                int sector = tmap[ix + w*iy];
                //printf("sector %d %d = %d\n", ix, iy, sector);
                // There can be small patches which were not used, so need to check for legality here
                // Again, should no longer be needed
                if( patch == 0 )  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
		    continue;
		images[i].sectors[patch][sector].Transform( p );  // change to global coordinates
		ix = ROUND(p.x);
		iy = ROUND(p.y);
                if( ix == 8557 && iy == 431 ) {
                    printf("Point %f %f in image %s\n", pts[j].x, pts[j].y, images[i].rname);
		    printf("Maps to pixel %d %d, image %d patch %d sector %d\n", ix, iy, i, patch, sector);
                    images[i].sectors[patch][sector].TPrint( stdout, "Transform is " );
                    images[i].  sinvs[patch][sector].TPrint( stdout, " Inverse  is " );
		    }
                if( ix < 0 || ix >= nx || iy < 0 || iy >= ny )
		    continue;  // outside the image
		//printf("%f i=%d ix,iy=%d %d\n", r, i, ix, iy);
		bool dbg = (2642 <= ix && ix <= 2644 && 3398 <= iy && iy <= 3400);
                if( dbg ) printf("\nDebugging point %d %d, r=%f, p= %f %f, pts[%d]= %f %f\n", ix, iy, r, p.x, p.y, j, pts[j].x, pts[j].y );
		uint32 nn = uint32(ix) + uint32(iy)*nx;   // index into array
		if( imap[nn] == 0 ) {
		    imap[nn] = k+1;                        // this pixel will be set by the kth relevant picture
		    PatchFrom[nn] =  patch;  // patch 'patch', sector 'sector'
                    SectorFrom[nn]= sector;
                    if( dbg ) printf("Not set, setting to %d %d %d %d\n", k, i, patch, sector);
                    }
		else { // already set, but we still might be better (because of finite
		   // grid in source, rounding to int, ordering ).
		    Point pt_us = Point(ix,iy);
		    Point pt_ot = pt_us;
		    // transform back into original frames.
		    int oi = relevant_images[imap[nn]-1];     // other image
		    int op = PatchFrom[nn];                   // other patch
                    int os = SectorFrom[nn];                  // other sector
		    images[oi].sinvs[op][os].Transform( pt_ot );
		    images[i].sinvs[patch][sector].Transform( pt_us );

		    double d_us = max(abs(pt_us.x-center.x), abs(pt_us.y-center.y));  // L-infinity norm
		    double d_ot = max(abs(pt_ot.x-center.x), abs(pt_ot.y-center.y));  // L-infinity norm
			    if( dbg ) printf("Already set to %d %d %d %d, d_us %f, d_ot %f\n", imap[nn]-1, oi, op, os, d_us, d_ot);
		    if( d_us < d_ot ) {
			if( dbg ) printf("It CAN happen #2... d_us= %f, d_ot= %f\n", d_us, d_ot);
			// Set it, even though it was previously set, since this is better
			imap[nn] = k+1;       // this pixel will be set by the ith picture
			PatchFrom[nn] = patch;
			SectorFrom[nn] = sector;
			}
		    }
		}
	    }
	}
    // Now compress the map/patch/sector array.  We assume less than 65536 combinations are used, so we can express
    // this as a 16 bit TIF.  The combinations start at 1, and run sequentially.
    temp.clear();               // clear the mapping file
    Triples.resize(1);          // and the triples; will fill up as they are found; first entry remains (0,0,0)
    for(uint32 i=0; i<imap.size(); i++) {
        if( imap[i] == 0 )
	    continue;       // pixel was never mapped
	int id = imap[i] + (PatchFrom[i] << 16) + (SectorFrom[i] << 24);  // create an int.
        map<int,int>::iterator it = temp.find(id);
        if( it == temp.end() ) { // add it
	    Triple t(imap[i], PatchFrom[i], SectorFrom[i]);
	    temp.insert(pair<int,int>(id, Triples.size()));
            printf("New combo %5d: %3d %3d %3d\n", Triples.size(), t.image, t.patch, t.sector);
            Triples.push_back(t);
	    it = temp.find(id);   // now should find it, since we just added it
	    }
        uint32 iy = i/nx;
        uint32 ix = i - nx*iy;
	if( ix == 8557 && iy == 431 ) {
	    printf("Compress map ix=%d iy=%d\n", ix, iy);
            printf(" imap %d, patch %d, sector %d\n", imap[i], PatchFrom[i], SectorFrom[i]);
            if( imap[i] > 0 )
		printf("Image %d\n", relevant_images[imap[i]-1] );
	    printf(" id = %d = %x\n", id, id);
            printf(" it->second = %d\n", it->second);
	    }
        imap[i] = it->second;
	}
    // write the new map file, if requested
    if( make_map ) {
        sprintf(fname,"%s/%s/map.%05d.png", region_dir.c_str(), inv_dir.c_str(), out_layer);
        printf("write Pixel mapping file %s\n", fname);
        Raster16ToPng16( fname, &imap[0], nx, ny );
        }
    else
        printf("User requested no map file.\n");
    printf("Draw the new image, with warp\n");

    // Now redraw the image, using the new map.  We'll write it into 'before', which is somewhat
    // confusing, but we already have it allocated.
    memset( &before[0], 0, nx * ny * sizeof(uint8) );
    // also, create a super-pixel map.  This needs to be 32 bit, even though that's a humongous
    // array, since there are typically 10K superpixels per image, and 9x9 arrays are typical,
    // so there are way more than 2^16 superpixel IDs.
    vector<uint32> spmap(nx*ny,0);
    for(int y=0; y < ny; y++) {
        if( (y & 0x3FF) == 0 ) {
            printf("."); fflush(stdout);
            }
        for(int x=0; x<nx; x++) {
            uint32 bi = uint32(x) + uint32(y)*nx;  // big index
	    uint16 indx = imap[bi];
	    uint16 from = Triples[indx].image;  // what image does this pixel come from?
            uint8 patch = Triples[indx].patch;  // what patch
            int sector =  Triples[indx].sector; // what sector
	    if( from == 0 )
		continue;  // no image sets this pixel
            int img = relevant_images[from-1];
	    Point p(x, y);
	    images[img].sinvs[patch][sector].Transform( p );
	    // of course, this should be in the image, but let's double check
	    if( p.x >= 0.0 && p.x < w-1 && p.y >= 0.0 && p.y < h-1 ) { // then we can interpolate
		int ix = int(p.x);
		int iy = int(p.y);
		double alpha = p.x-ix;
		double beta  = p.y-iy;
		int nn = ix + w*iy;  // index into original image
		uint8* pic = images[img].raster;
		double pix =
		     (1-alpha)*(1-beta) * pic[nn] +
			alpha *(1-beta) * pic[nn+1] +
		     (1-alpha)*   beta  * pic[nn+w] +
			alpha *   beta  * pic[nn+w+1];
		before[bi] = ROUND(pix);
                // for the super-pixel map we want nearest, not interpolation.
                if( p.x - ix >= 0.5 )
		    ix++;
                if( p.y - iy >= 0.5 )
		    iy++;
                int px = images[img].spmap[ix + w*iy];
                if( px != 0 )  // 0 valued pixels are unassigned, and not translated
                    px += images[img].spbase;
                spmap[bi] = px;
                //if( spmap[x+y*nx] == 65536 ) {  // Just debugging
                    //printf("w=%d h=%d nx=%d ny=%d patch=%d sector=%d\n", w, h, nx, ny, patch, sector);
		    //printf("Set to %d.  x=%d y=%d ix=%d iy=%d img=%d images[img].spbase=%d images[img].spmap[ix + w*iy] = %d\n",
                     //px, x, y, ix, iy, img, images[img].spbase, images[img].spmap[ix + w*iy] );
                    //exit( 42 );
		    //}
                if( Debug && px == 0 )
		    before[bi] = 255;
		}
	    } // to aid debugging, draw lines at the image boundaries in the before image.
         }
    if( Annotate ) {
	vector<Point> edges;
	for(int x=0; x<w; x++) {
	    edges.push_back(Point(x, 0.0));
	    edges.push_back(Point(x, h-1));
	    }
	for(int y=0; y<h; y++) {
	    edges.push_back(Point(0.0, y));
	    edges.push_back(Point(w-1, y));
	    }
	printf("%d edge pixels\n", edges.size());
	// OK, now transform this edge list by each of the transforms, then color it in.
	for(int i=0; i<images.size(); i++) {
	    if( images[i].layer != out_layer ) // only want edges on current layer
		continue;
	    for(int k=1; k<images[i].tf.size(); k++) {
		if( images[i].tf[k].det() == 0.0 )
		    continue;
		for(int j=0; j<edges.size(); j++) {
		    Point p = edges[j];
		    images[i].tf[k].Transform( p );
		    int ix = ROUND(p.x);
		    int iy = ROUND(p.y);
		    if( 0 <= ix && ix < nx && 0 <= iy && iy < ny ) {
                        uint32 bi = uint32(ix) + nx*uint32(iy);
			before[ix + nx*iy] = 255;
                        }
		    }
		}
	    }
        }
    //sprintf( fname, "after.%05d.tif", out_layer );
    //Raster8ToTif8( fname, &before[0], nx, ny );
    sprintf(fname,"%s/%s/after.%05d.png", region_dir.c_str(), gray_dir.c_str(), out_layer);
    if( make_flat )
        Raster8ToPng8( fname, &before[0], nx, ny );
    if( make_tiles )
        WriteImageTiles(rav_dir, out_layer, &before[0], nx, ny);

    // If any boundarymap files were found, write the boundary map
    // Again, use the array 'before'
    memset( &before[0], 0, nx * ny * sizeof(uint8) );
    for(int x=0; x<nx && AnyBMap; x++) {
	for(int y=0; y < ny; y++) {
            uint32 bi = uint32(x) + uint32(y)*nx;  // big index
	    uint16 indx = imap[bi];
	    uint16 from = Triples[indx].image;  // what image does this pixel come from?
            uint8 patch = Triples[indx].patch;  // what patch
            int sector =  Triples[indx].sector; // what sector
	    if( from == 0 )
		continue;  // no image sets this pixel
            int img = relevant_images[from-1];
	    Point p(x, y);
	    images[img].sinvs[patch][sector].Transform( p );
	    // of course, this should be in the image, but let's double check
	    if( p.x >= 0.0 && p.x < w-1 && p.y >= 0.0 && p.y < h-1 ) { // then we can interpolate
		int ix = int(p.x);
		int iy = int(p.y);
		double alpha = p.x-ix;
		double beta  = p.y-iy;
		int nn = ix + w*iy;  // index into original image
		uint8* pic = images[img].bmap;
		double pix =
		     (1-alpha)*(1-beta) * pic[nn] +
			alpha *(1-beta) * pic[nn+1] +
		     (1-alpha)*   beta  * pic[nn+w] +
			alpha *   beta  * pic[nn+w+1];
		before[x+y*nx] = ROUND(pix);
	        }
	    }
	}
    if( AnyBMap ) {
        sprintf(fname,"%s/%s/bmap.%05d.png", region_dir.c_str(), bmap_dir.c_str(), out_layer);
	if( make_flat )
            Raster8ToPng8( fname, &before[0], nx, ny );
	}

    // ------------------------------------------------ write the mapping text file  --------------------
    // This is the text file that contains the transformations that describe how each pixel got there.
    sprintf(fname,"%s/%s/mapping.%05d.txt", region_dir.c_str(), inv_dir.c_str(), out_layer);
    FILE *ftxt = FileOpenOrDie( fname, "w" );
    //write the image names, and their indexes
    for(int k=0; k<relevant_images.size(); k++) {      // for each picture
        int i = relevant_images[k];
	fprintf(ftxt,"IMAGE %d '%s'\n", k, images[i].rname);
        }
    for(int k=1; k<Triples.size(); k++) {
        int from = Triples[k].image;          // index into the relevant pictures array
        int image = relevant_images[from-1];  // actual image number
        int patch = Triples[k].patch;
        int sector = Triples[k].sector;
	TAffine t = images[image].sinvs[patch][sector];  // transform from global to image
        fprintf(ftxt, "TRANS %d %d ", k, from-1);
        t.TPrint( ftxt );
        }
    fclose(ftxt);


    // Write the super pixel map file.  First make sure it has all consecutive numbers, in the interest of
    // efficiency in Raveler.
    int SPmax = 0;
    vector<int> FinalMapping;
    if( RenumberSuperPixels )
        RemapSuperPixelsOneImage( &spmap[0], nx, ny, SPmax, FinalMapping );
    printf("SPmax now %d\n", SPmax);
    sprintf(fname,"%s/%s/sp.%05d.png", region_dir.c_str(), sp_dir.c_str(), out_layer);
    for(uint32 i=0; i<nx*ny; i++)  // set the transparency
	spmap[i] = spmap[i] | 0xFF000000;
    if( make_flat )
        Raster32ToPngRGBA( fname, &spmap[0], nx, ny );



    for(uint32 i=0; i<nx*ny; i++)  // unset the transparency for tiles
	spmap[i] = spmap[i] & 0xFFFFFF;
    if( make_tiles ) {
        WriteSPTiles( rav_dir, out_layer, &spmap[0], nx, ny );
	}
    // Compute the bounds for this layer.  First find the biggest number
    uint32 biggest = 0;
    for(int i=0; i<nx*ny; i++)
	biggest = max(biggest, spmap[i]);
    // Now make the bounds and counts
    {
    vector<int> volume(biggest+1,0);
    vector<int> xmin(biggest+1,nx);
    vector<int> ymin(biggest+1,ny);
    vector<int> xmax(biggest+1,-1);
    vector<int> ymax(biggest+1,-1);
    for(int y=0; y<ny; y++) {
	for(int x=0; x<nx; x++) {
	    uint32 n = uint32(x) + uint32(y)*uint32(nx);
            int p = spmap[n];  // pixel value
            volume[p]++;
            xmin[p] = min(xmin[p], x);
            ymin[p] = min(ymin[p], y);
            xmax[p] = max(xmax[p], x);
            ymax[p] = max(ymax[p], y);
	    }
	}

    sprintf( fname, "bounds.%05d", out_layer );
    FILE *fb = FileOpenOrDie( fname, "w" );

    for(int i=0; i<=biggest; i++) {
        if( volume[i] > 0 )
	    fprintf(fb,"%d %7d %8d %d %d %d %d\n", out_layer, i, xmin[i], ny-1-ymax[i], xmax[i]-xmin[i]+1, ymax[i]-ymin[i]+1, volume[i]);
	}
    fclose(fb);
    }


    bool experiment = false;
    if( experiment ) {
        vector<uint32> copy(nx*ny,0);
	for(int i=0; i<nx*ny; i++)
	    copy[i] = 0xFF000000 + spmap[i];
        printf("try to write the small array as a .png file\n"); fflush(stdout);
        Raster32ToPngRGBA( "try32bita.png", &copy[0], nx, ny );
        int bigx = 45000;
	int bigy = 47000;
        printf("try to allocate a huge array\n"); fflush(stdout);
        vector<uint32> big(bigx*bigy);
        printf("try to fill the huge array\n"); fflush(stdout);
        int p1 = nx*ny;
        int p2 = bigx*bigy;
        for(int i=0; i<p2; i++)
	    big[i] = copy[i % p1];
        printf("try to write the huge array as a .png file\n"); fflush(stdout);
        Raster32ToPngRGBA( "try32bitb.png", &big[0], bigx, bigy );
        printf("done writing the huge array as a .png file\n"); fflush(stdout);
        }

    // Write the file that specifies how the SP IDs in each image are mapped to the SP IDs in the
    // larger output tile.
    // For example, image 1 might have IDs 1-8000.  These are mapped to 1-8000 in the output
    //              image 2 might have IDS 1-7500.  These are mapped to 8001-15500 in the output
    //              image 3 might have IDs 1-9100   These are mapped to 15501-24601 in the output
    //              etc.
    sprintf(fname,"%s/map-sp.%05d.txt", region_dir.c_str(), out_layer);
    FILE *fmap = FileOpenOrDie( fname, "w" );
    // For every image that has a non-blank sp map;
    for(int j=0; j<relevant_images.size(); j++) {
	int i = relevant_images[j];
	if( images[i].spmap != NULL && images[i].spmap != BlankSPMap ) { // one was defined
            fprintf(fmap,"Map '%s'\n", images[i].spname);
            fprintf(fmap,"NUM_MAP %d\n", images[i].SPmapping.size()-1);  // number of entries to follow
	    for(int k=1; k<images[i].SPmapping.size(); k++) {
		int m = images[i].SPmapping[k];
		if( m == 0 )
		    fprintf(fmap,"%d -1\n", k);  // -1 => not in the original image
		else {
		    //printf("%6d maps to %d (global set of maps)\n", k, m);
                    // this could go off the end of the Final Mapping array, if the biggest number
                    // assigned is not used in the final image.  Test for this:
                    int n =  (m < FinalMapping.size()) ? FinalMapping[m] : 0;
		    if( n == 0 )
			fprintf(fmap, "%d -2\n", k); // -2 => not used in final image
		    else
			fprintf(fmap, "%d %d\n", k, n);
		    }
		}
	    }
	}
    fclose(fmap);

    // Now, for every image that was used, look and see if it had annotation.  If so, transform
    // the coordinates and write a file for raveler.
    vector<TBar> tbars;
    vector<Bookmark> bookmarks;
    for(int j=0; j<relevant_images.size(); j++) {
	int i = relevant_images[j];
        // find where the super pixel map is - that is where the annotations will be
        char *p = NULL;
        char *root;
        if( images[i].spname != NULL ) {
	     root = strdup(images[i].spname);
             p = strrchr(root, '/');
             }
        if( p == NULL ) {
	    printf("Oops - no '/' in super-pixel path '%s'.  Try fold mask.\n", images[i].spname == NULL ?"none":images[i].spname);
            root = strdup(images[i].fname);
            p = strrchr(root, '/');
            if( p == NULL ) {
	        printf("Oops - no '/' in fold mask path either '%s'.\n", root);
	        exit( 42 );
                }
	    }
        *(++p) = '\0';  // truncate at trailing '/', leaving directory name
	LookForAndTransform("annotations-body.json",      root, out_layer, i, images, Warp, tmap, tbars); // nop for now
	LookForAndTransform("annotations-bookmarks.json", root, out_layer, i, images, Warp, tmap, bookmarks);
	LookForAndTransform("annotations-synapse.json",   root, out_layer, i, images, Warp, tmap, tbars);
	LookForAndTransform("session-metadata.json",      root, out_layer, i, images, Warp, tmap, tbars); // nop for now.
	}
    // Now, write the resulting synapses out.
    char tbar_name[2048];
    sprintf(tbar_name, "annot-syn-%d.json", out_layer);
    printf("Writing to '%s'\n", tbar_name);
    FILE *fw = FileOpenOrDie( tbar_name, "w" );
    fprintf(fw, "{\"data\": [\n");
    for(int i=0; i<tbars.size(); i++) {
	fprintf(fw, "{\"T-bar\": {\"status\": \"%s\", \"confidence\": %f, \"body ID\": %d, \"location\": [%f, %f, %d]}, \"partners\": [\n",
	 tbars[i].status.c_str(), tbars[i].confidence, tbars[i].body_id, tbars[i].pt.x, ny-tbars[i].pt.y, tbars[i].z);
	for(int j=0; j<tbars[i].partners.size(); j++) {
	    bool last = (j == tbars[i].partners.size()-1);
	    fprintf(fw, " {\"confidence\": %f, \"body ID\": %d, \"location\": [%f, %f, %d]}%c\n",
	     tbars[i].partners[j].confidence, tbars[i].partners[j].body_id, tbars[i].partners[j].pt.x, ny-tbars[i].partners[j].pt.y, tbars[i].partners[j].z,
	     last ? ' ': ',' );  // comma separated if not last
	    }
	bool last = (i == tbars.size()-1);
	fprintf(fw, " ]}%c\n", last ? ' ': ',');  // comma separated if not last
	}
    // ']' ends the list of bookmarks, then the metadata, then '}' ends the 'data' section
    int version = 1;
    fprintf(fw, "], \"metadata\": {\"description\": \"synapse annotations\", \"file version\": %d}}\n", version);
    fclose(fw);
    // Now, write the resulting bookmarks out.
    char book_name[2048];
    sprintf(book_name, "annot-book-%d.json", out_layer);
    printf("Writing to '%s'\n", book_name);
    fw = FileOpenOrDie( book_name, "w" );
    fprintf(fw, "{\"data\": [\n");
    for(int i=0; i<bookmarks.size(); i++) {
	fprintf(fw, "{\"text\": \"%s\", \"body ID\": %d, \"location\": [%f, %f, %d]}",
         bookmarks[i].text.c_str(), bookmarks[i].body_id, bookmarks[i].pt.x, ny-bookmarks[i].pt.y, bookmarks[i].z);
        bool last = (i == bookmarks.size()-1);
        fprintf(fw, "%c\n", last ? ' ': ',');  // comma separated if not last
        }
    // ']' ends the list of bookmarks, then the metadata, then '}' ends the 'data' section
    version = 1;
    fprintf(fw, "], \"metadata\": {\"description\": \"bookmarks\", \"file version\": %d}}\n", version);
    fclose(fw);

    printf("WARNINGS:  %d hit the edge, %d had bad correlation\n", nwarn_edge_interp, nwarn_bad_corr);
    }


VMStats( stdout );

return 0;
}


