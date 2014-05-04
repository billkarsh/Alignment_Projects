
// Make a mosaic


#include	"LinEqu.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"TAffine.h"
#include	"Draw.h"

#include	<string.h>

#include	<map>
#include	<set>
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

vector<glob_spot> spots;

class image {
public:
	uint8*					raster;		// gray scale image
	uint8*					foldmap;	// fold map
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

static bool dp = false;				// debug print
static int nwarn_edge_interp = 0;	// interpolation cannot be done since max is on edge of region
static int nwarn_bad_corr = 0;		// could not find a good correlation
static int nwarn_good_on_edge = 0;	// worst case - good correlation, but on edge






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
/* PointsInRing -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Produce the points at radius R (a square) around the center,
// in steps < 1 in size. Might include duplicates.
//
void PointsInRing(vector<Point> &pts, Point &center, double r, uint32 w, uint32 h)
{
double lx = center.x-r;  // left and right X
double rx = center.x+r;
double by = center.y-r;  // bottom and top Y
double ty = center.y+r;
vector<Point> ppts; // potential points
for(double d=-r; d <= r; d+= 0.9) {
    ppts.push_back(Point(lx, center.y+d)); // left side
    ppts.push_back(Point(rx, center.y+d)); // left side
    ppts.push_back(Point(center.x+d, by)); // bottom
    ppts.push_back(Point(center.x+d, ty)); // top
    }
// Now keep only those within the image (can happen if image is a rectangle, not a square)
for(int i=0; i<ppts.size(); i++)
    if( ppts[i].x >= 0.0 && ppts[i].x <= w-1 && ppts[i].y >= 0.0 && ppts[i].y <= h-1 )
	pts.push_back(ppts[i]);
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
	vector<image>		&images,
	Point				pt,
	int					i1,
	int					patch1,
	int					out_layer )
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
/* PseudoAngle --------------------------------------------------- */
/* --------------------------------------------------------------- */

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

/* --------------------------------------------------------------- */
/* FillInHolesInFoldmap ------------------------------------------ */
/* --------------------------------------------------------------- */

// We consider every 0 in the foldmap, and decide whether to
// fill it in. We will change each 0 to 255 as we examine it
// so we only see it once. After we see each hole, we will
// either fill it in, or add it to the list of 0s to be
// restored at the end.
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
	set<int> boundary;  // all the values we find on the boundary
        vector<int> pixels; // 0 valued pixels in this area
	while( !st.empty() ) {
	    int j = st.top(); st.pop();
            if( map[j] == 0 ) {
	        map[j] = 255;  // so we won't get it again
                pixels.push_back(j);
	        int y = j / w;
	        int x = j - w * y;
	        if( x-1 >= 0 ) st.push(j-1); else boundary.insert(-1);
            if( x+1 <  w ) st.push(j+1); else boundary.insert(-1);
            if( y-1 >= 0 ) st.push(j-w); else boundary.insert(-1);
            if( y+1 <  h ) st.push(j+w); else boundary.insert(-1);
		}
            else {// map value is not 0, add what we ran into to the set of boundary values
		if( map[j] != 255 ) boundary.insert( map[j] );
		}
	    }
        if( boundary.size() == 1 && *(boundary.begin()) != -1 ) { // then we have a hole
	    int pval = *(boundary.begin());
	    for(int k=0; k<pixels.size(); k++)
		map[pixels[k]] = pval;
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
bool Annotate = true;
vector<char *>noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    // process arguments here
    if( argv[i][0] != '-' )
	noa.push_back(argv[i]);
    else if( strcmp(argv[i],"-d") == 0 )
	Debug = true;
    else if( strncmp(argv[i],"-warp", 5) == 0 )
	Warp = true;
    else if( strcmp(argv[i], "-nf") == 0 )
	FoldMasks = false;
    else if( strcmp(argv[i], "-na") == 0 )
	Annotate = false;
    else if( strncmp(argv[i], "-dms=",5) == 0 ) {
	DontMoveStrength = atof(argv[i]+5);
	}
    else
	printf("Ignored option '%s'\n", argv[i]);
    }
if( Warp )
    printf("Don't Move Strength = %f\n", DontMoveStrength);
if( noa.size() < 1 ) {
    printf("Usage: mos <simple file> [region xmin,ymin,dx,dy] [minlayer,maxlayer] [options]\n");
    exit( 42 );
    }

FILE	*fp = FileOpenOrDie( noa[0], "r" );

int ni = 0;  // number of images
vector<char *>dnames;   // directory names
vector<int>   lnums;    // parallel vector of layer numbers
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

// read the file.  Here we assume all images are the same size, so if w and/or h are non-zero, we don't really need to read the
// images, since we already read at least one.
double xmin = BIG, ymin = BIG;
double xmax = -BIG, ymax = -BIG;
vector<double>	x1, y1, x2, y2;
vector<int>		z1, z2;
CLineScan		*ls = new CLineScan;
for(;;){
    if( ls->Get( fp ) <= 0 )
		break;
    if( strncmp(ls->line,"DIR",3) == 0 ) {      // simply store directory string
        char *anum = strtok(ls->line+3," ");  // and layer number in parallel
        char *dname = strtok(NULL," \n");     // vectors
        dnames.push_back(strdup(dname));
        lnums.push_back(atoi(anum));
        }
    else if( strncmp(ls->line,"FOLDMAP",7) == 0 ) {
	char *tname = strtok(ls->line+7," '\n");
	char *mname = strtok(NULL, " '\n");
        if( tname == NULL || mname == NULL ) {
	    printf("Not expecting NULL in FOLDMAP parsing.\n");
	    exit( 42 );
	    }
        map<string,int>::iterator lookup = imap.find(string(tname));
        if( lookup != imap.end() )
	     continue;  // already seen this one
        // at this point, all we need is w and h, so toss image as soon as we read it.
        image ii;
        if( w == 0 ) {  //read a file, then free it, just to set w and h
            ii.raster = Raster8FromAny( tname, w, h );
            RasterFree(ii.raster);
	    }
	ii.raster = NULL;
	ii.foldmap = NULL;
	ii.spname = NULL;      // No super-pixel map is defined for this tile
	ii.spmap = NULL;
	ii.spbase = 0;
	ii.rname = strdup(tname);
        ii.w = w; ii.h = h;
        // find the layer number
        int k;
        for(k=0; k<dnames.size(); k++)
	    if( strstr(tname, dnames[k]) != NULL )
	        break;
        if( k >= dnames.size() ) {
	    printf("Oops - '%s' not in layer directory.  Using 0.\n", tname);
            ii.layer = 0;
            exit( 42 );
	    }
        else // found it
            ii.layer = lnums[k];
        ii.foldmap = NULL;
	ii.fname = strdup(mname);
        //ii.foldmap = Raster8FromAny( mname, w, h );
        // if the foldmap was blahblahblah.tif, also look for blahbladblahd.tif, and use that instead.
        // (this is the fold mask used for drawing, as opposed to that used for correlation.)
        // If layers are specified, do this only for layers that are used, since it requires a file access
        // and hence is slow.
	char *suf = strstr(mname, ".tif");
        if( suf != NULL ) {
            // see if it's a layer we care about.
            if( lspec1 == -1 || (lspec1 <= ii.layer && ii.layer <= lspec2) ) {
                *suf = '\0';
                char temp[2048];
                sprintf(temp, "%sd.tif", mname);
	        *suf = '.';  // put name back, in case we don't find it.
                FILE *fd = fopen(temp, "r");
                if( fd ) {
                    fclose(fd);
		    printf("Swapping to drawing file '%s'\n", temp);
                    free(ii.fname);
                    ii.fname = strdup(temp);
		    }
		}
	    }
        // is this a layer we care about?  If layer range not specified, or specified and in it, save
        if( lspec1 < 0 || (lspec1 <= ii.layer && ii.layer <= lspec2)  ) {
            lowest  = min(lowest,  ii.layer);
            highest = max(highest, ii.layer);
            images.push_back(ii);
            imap[string(tname)] = images.size() - 1;
	    }
        }
    else if( strncmp(ls->line,"TRANSFORM",9) == 0 ) {
        double a,b,c,d,e,f;
        char name[2048];
        if( sscanf(ls->line+9, "%s %lf %lf %lf %lf %lf %lf", name, &a, &b, &c, &d, &e, &f) != 7 ) {
            printf("Not expecting this in TRANSFORM: %s", ls->line);
	    break;
            }
        TAffine tf( a, b, c, d, e, f );
        char *fname = strtok(name," ':");
        //printf("File '%s'\n", fname);
        map<string,int>::iterator imit = imap.find(string(fname));  // imap iterator
        if( imit == imap.end() ) {
	    // this is now normal with single layer image generation
	    //printf("File in TRANSFORM statement has no FOLDMAP - ignored.\n");
	    continue;
	    }
        int k = imit->second;
        printf("Image = %d\n", k);
        char *rest = strtok(NULL," :'\n"); // get the rest of the string, now of the form ::123
        int patch = atoi(rest);
        printf("rest = '%s', patch = %d\n", rest, patch);
        // Make sure vector is big enough
        if( images[k].tf.size() <= patch ) {
             images[k].tf.resize(patch+1, TAffine(0,0,0,0,0,0));  // initialize to an illegal transform
             images[k].inv.resize(patch+1);
             }
        images[k].tf[patch] = tf;
	}
    else if( strncmp(ls->line,"SPMAP",5) == 0 ) {
        char name[2048], where[2048];
        if( sscanf(ls->line+5, "%s %s", name, where) != 2 ) {
            printf("Not expecting this in SPMAP: %s", ls->line);
	    break;
            }
        char *fname = strtok(name," ':");
        //printf("File '%s'\n", fname);
        map<string,int>::iterator imit = imap.find(string(fname));  // imap iterator
        if( imit == imap.end() ) {
	    // this is now normal with single layer image generation
	    printf("File in SPMAP statement has no FOLDMAP - ignored.\n");
	    continue;
	    }
        int k = imit->second;
        printf("Image = %d\n", k);
        images[k].spname = strdup(where);
        }
    else if( strncmp(ls->line,"MPOINTS",7) == 0 ) {
        double a,b,c,d;
        int za, zc;
        if( sscanf(ls->line+7, "%d %lf %lf %d %lf %lf", &za,  &a, &b, &zc, &c, &d) != 6 ) {
            printf("Not expecting this: %s", ls->line);
	    break;
            }
        if( lspec1 < 0 || (lspec1-1 <= za && za <= lspec2+1) || (lspec1-1 <= zc && zc <= lspec2) ) {
            z1.push_back(za); x1.push_back(a); y1.push_back(b);
            z2.push_back(zc); x2.push_back(c); y2.push_back(d);
	    }
	}
    else if( strncmp(ls->line,"BBOX",4) == 0 ) {
        if( sscanf(ls->line+4, "%lf %lf %lf %lf", &xmin, &ymin, &xmax, &ymax) != 4 ) {
	     printf("Bad BBOX statement %s\n", ls->line);
	     return 42;
             }
	}
    else if( strncmp(ls->line,"IMAGESIZE",9) == 0 ) {
	// ignore for now - so we do not read to compute BBs
	}
    else {
	printf("Unknown line '%s'", ls->line);
	return 42;
        }
}

delete ls;

if( images.size() == 0 ) {
    printf("No images in input\n");
    return 42;
    }
// Now find the bounding box in global space, if not already specified
if( xmin > BIG/2 ) {
    for(int i=0; i<images.size(); i++) {
	for(int k=1; k<images[i].tf.size(); k++) {
	    double det = images[i].tf[k].det();  // no legal transform will have a 0 determinant
	    if( det == 0.0 )
		continue;
	    vector<Point> pts;
	    pts.push_back(Point(0.0,0.0));
	    pts.push_back(Point(w-1,0.0));
	    pts.push_back(Point(w-1, h-1));
	    pts.push_back(Point(0.0, h-1));
	    for(int j=0; j<4; j++) {
		images[i].tf[k].Transform( pts[j] );
		xmin = fmin(xmin, pts[j].x);
		ymin = fmin(ymin, pts[j].y);
		xmax = fmax(xmax, pts[j].x);
		ymax = fmax(ymax, pts[j].y);
		}
	    }
	}
    }
printf("Bounds of global image are x=[%f %f] y=[%f %f]\n", xmin, xmax, ymin, ymax);

// Modify all transformations so the min is (0,0).  , then apply x0,y0.
// Always change by an integer amount, to make later layer-to-layer alignment easier
int xfl = int(floor(xmin));
int yfl = int(floor(ymin));
xmin -= xfl; xmax -= xfl;              // change the bounding box
ymin -= yfl; ymax -= yfl;
for(int i=0; i<images.size(); i++) {
    for(int k=1; k<images[i].tf.size(); k++) {
	if( images[i].tf[k].det() == 0.0 )
	    continue;
        images[i].tf[k].t[2] -= (xfl+x0);          // and each transform
        images[i].tf[k].t[5] -= (yfl+y0);
        images[i].inv[k].InverseOf( images[i].tf[k] );
	}
    }
printf("Bounds of global image are x=[%f %f] y=[%f %f]\n", xmin, xmax, ymin, ymax);
if( xsize != -1 ) {
    xmax = xsize;
    ymax = ysize;
    }

// Print transforms, post offset
//for(int i=0; i<images.size(); i++) {
    //for(int k=1; k<images[i].tf.size(); k++) {
        //printf("Post offset: image %d, layer %d tf[%d] = %f %f %f   %f %f %f\n", i, images[i].layer, k,
	 //images[i].tf[k].t[0], images[i].tf[k].t[1], images[i].tf[k].t[2],
	 //images[i].tf[k].t[3], images[i].tf[k].t[4], images[i].tf[k].t[5]);
	//}
    //}
// fix the corresponding points
for(int i=0; i<x1.size(); i++) {
    x1[i] -= (xfl+x0); x2[i] -= (xfl+x0);
    y1[i] -= (yfl+y0); y2[i] -= (yfl+y0);
    }

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

// Since some (many) tiles will not have super-pixel maps, create a single blank map they all
// can point to.
uint16 *BlankSPMap = (uint16*)malloc(w*h*sizeof(uint16));
memset( BlankSPMap, 0, w * h * sizeof(uint16) );

// The following huge loop goes through the layers one at a time, producing 3 outputs
//   (a) a 'before.%d.png', which is the results of simply mapping the input affine
//   (b) an 'after.%d.png', which is the result after edge merging
//   (c) an 'spmap.%d.png', which is the super-pixel map
for(int out_layer = lowest; out_layer <= highest; out_layer++) {  //keep going until we find nothing...

    printf("starting to generate layer %d...\n", out_layer);
    // start by tossing out previous images, if any
    for(int i=0; i<images.size(); i++) {
	if( images[i].raster != NULL )
	    RasterFree(images[i].raster);
	if( images[i].foldmap != NULL )
	    RasterFree(images[i].foldmap);
	if( images[i].spmap != NULL && images[i].spmap != BlankSPMap )
	     free(images[i].spmap);
	}
    // Find the relevant images, and read their data.
    vector<int> relevant_images;
    // Since each superpixel map is numbered independently per image,
    // map all of these to non-overlapping ranges.
    int MaxSPUsed = 0;
    for(int i=0; i<images.size(); i++) {      // for each picture
        if( images[i].layer != out_layer )     // in the correct layer
	    continue;
        relevant_images.push_back(i);
	uint32 ww, hh;
	if( images[i].raster == NULL )
	    images[i].raster  = Raster8FromAny( images[i].rname, ww, hh );	// Normalize param=true
	if( images[i].foldmap == NULL ) {
	    if( FoldMasks ) {
		images[i].foldmap = Raster8FromAny( images[i].fname, ww, hh );
		FillInHolesInFoldmap(images[i].foldmap, ww, hh);
		FoldmapRenumber(images[i]);
		}
	    else {  // no foldmasks wanted; just create a field of 1s.
		images[i].foldmap = (uint8*)malloc(w*h*sizeof(uint8));
		for(int k=0; k<w*h; k++)
		    images[i].foldmap[k] = 1;
		}
	    }
	if( images[i].spmap == NULL ) {
            uint16 *test = NULL;
            if( images[i].spname == NULL ) {
		printf("No SPmap defined for image '%s'\n", images[i].rname);
		test = BlankSPMap;
		}
	    else {
			uint32	xw, xh;
			test = Raster16FromPng( images[i].spname, xw, xh, stdout );
	        if( test == NULL ) {
		    printf("No luck reading '%s'.\n", images[i].spname);
	            test = BlankSPMap;
		    }
	        else {
                    images[i].spbase = MaxSPUsed;
		    RemapSuperPixelsOneImage(test, w, h, MaxSPUsed, images[i].SPmapping);
                    printf("MaxSPUsed is now %d\n", MaxSPUsed);
		    }
		}
	    images[i].spmap = test;
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
    // find the center of each picture.  May not be integer
    Point center((w-1)/2.0, (h-1)/2.0);
    for(double r=0.0; r <= max(center.x,center.y)+0.0001; r += delta_image_space) {
        //compute a ring of radius R
	vector<Point> pts;
	PointsInRing(pts, center, r, w, h);
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
                if( patch == 0)  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
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
		    int oi = imap[nn]-1;     // other image
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
    // Now compress the map/patch array.  We assume less than 65536 combinations are used, so we can express
    // this as a 16 bit TIF.  The combinations start at 1, and run sequentially.
    map<int,int> temp;
    vector<Triple> Triples(1);  // will fill up as they are found
    for(int i=0; i<imap.size(); i++) {
        if( imap[i] == 0 )
	    continue;	   // was never set - leave as is
	int id = imap[i] + (PatchFrom[i] << 16);  // create an int.
        map<int,int>::iterator it = temp.find(id);
        if( it == temp.end() ) { // add it
	    Triple t(imap[i], PatchFrom[i], 0);
	    temp.insert(pair<int,int>(id, Triples.size()));
            printf("New combo %5d: %d %d\n", Triples.size(), t.image, t.patch);
            Triples.push_back(t);
	    it = temp.find(id);   // now should find it, since we just added it
	    }
        imap[i] = it->second;
	}

    // write the map; only used for debugging for now
    char fname[256];
    if( Debug ) {
        printf("Try writing 16 bit map\n");
        Raster16ToPng16("test.png", &imap[0], nx, ny);
        sprintf(fname,"map.%d.png", out_layer);
        printf("Writing  8 bit map in png format to %s\n", fname);
        vector<uint8> copy(nx*ny, 0);  // make an 8 bit copy
        for(int i=0; i<nx*ny; i++)
	    copy[i] = imap[i];
        Raster8ToPng8(fname, &copy[0], nx, ny);
        //Raster8ToTif8( fname, &(copy[0]), nx, ny );
        sprintf(fname,"pf.%d.tif", out_layer);
        //Raster8ToTif8( fname, &(PatchFrom[0]), nx, ny );
        }

    printf("Writing map.tif done - starting 'before' image\n");
    // Now, create a 'before' picture with seams
    vector<uint8>before(nx*ny,0);
    for(int x=0; x<nx; x++) {
	for(int y=0; y < ny; y++) {
            uint16 indx = imap[x + nx*y];
	    uint16 from = Triples[indx].image;  // what image does this pixel come from?
            uint8 patch = Triples[indx].patch;
	    if( from == 0 )
		continue;  // no image sets this pixel
	    int img = relevant_images[from-1];
	    Point p(x, y);
	    images[img].inv[patch].Transform( p );
	    // of course, this should be in the image, but let's double check
	    if( p.x >= 0.0 && p.x < w-1 && p.y >= 0.0 && p.y < h-1 ) { // then we can interpolate
		double pix = InterpolatePixel( p.x, p.y, images[img].raster, w );
		before[x + nx*y] = ROUND(pix);
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
		    if( 0 <= ix && ix < nx && 0 <= iy && iy < ny )
			before[ix + nx*iy] = 255;
		    }
		}
	    }
	//also to aid in debugging, draw circles at correspondence points
	for(int i=0; i<x1.size(); i++) {
	    if( z1[i] == out_layer || z2[i] == out_layer ) {
		DrawCircle(&before[0], nx, ny, x1[i], y1[i], 10.0);
		DrawCircle(&before[0], nx, ny, x2[i], y2[i], 10.0);
		}
	    }
        }
    sprintf(fname,"before.%d.tif", out_layer);
    //Raster8ToTif8( fname, &(before[0]), nx, ny );
    sprintf(fname,"before.%d.png", out_layer);
    Raster8ToPng8(fname, &before[0], nx, ny);

    // if only simple output is needed, we are done with this layer
    if( !Warp )
	continue;

    // Break each image into triangles
    const int N=16;  // number of triangles per image

    // create the control points.  Triangles numbered CCW from angle 0.  Vertex N is at the center.
    // There is just one list that applies to all images.
    vector<Point> vtx(N+1);
    double ctrx = double(w-1)/2.0;
    double ctry = double(h-1)/2.0;
    vtx[N] = Point(ctrx, ctry);
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
       vtx[i].x = x;
       vtx[i].y = y;
       }
    for(int i=0; i<vtx.size(); i++)
        printf("vtx[%d] = (%f,%f)\n", i, vtx[i].x, vtx[i].y);

    // now, create a map that maps each pixel of an image into the relevent triangle.
    vector<uint8> tmap(w*h,0);
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
    //Raster8ToTif8( "tmap.tif", &(tmap[0]), w, h );

    // Now, for each patch of each image create N triangles using N+1 pts each, all transferred
    // to the global reference frame.
    vector<triangle> tris;
    vector<Point> gvtx;
    for(int i=0; i<images.size(); i++) {
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
        Invert3x3Matrix(tris[k].a, a);
        }

    // Now write the triangles out for debugging
    if( Debug )
        WriteTriangles("tris", tris, gvtx);

    // Then, look for places where adjacent pixels come from different pictures.
    // For each such spot (or a reasonable subset, find the delta needed make them line up.
    const int delta = 40;
    for(int ty=0; ty <ny; ty += delta) {
	printf("ty=%d\n", ty);
        int Ntrans = 0;   // count number of transitions
	for(int tx = 0; tx+1 < nx; tx++) {
	    int ind = ty*nx + tx;
            int indx1 = imap[ind];
            int indx2 = imap[ind+1];
            int image1 = Triples[indx1].image;
            int image2 = Triples[indx2].image;
	    if( image1 != image2 && image1 > 0 && image2 > 0 ) { // found a boundary
		Ntrans++;
		int i1 = relevant_images[image1-1];
		int i2 = relevant_images[image2-1];  // point is on boundary of image i1 and i2
		//i2 = i1;  // temp, just to make sure we get 0,0
		Point p1(tx,ty);
		CommonPoint(images, p1,i1, Triples[indx1].patch, out_layer);
		}
	    }
        if( Ntrans > 100 ) {
	    printf("Too many transitions (%d) in map!\n", Ntrans);
            if( !Debug ) {
		printf("Writing map as 'test.png'\n");
        	Raster16ToPng16("test.png", &imap[0], nx, ny);
		}
            return 42;
	    }
	}
    for(int tx=0; tx <nx; tx += delta) {
	printf("tx=%d\n", tx);
        int Ntrans = 0;   // count number of transitions
	for(int ty = 0; ty+1 < ny; ty++) {
	    int ind = ty*nx + tx;  // note that ind+nx is up by one in y
            int indx1 = imap[ind];
            int indx2 = imap[ind+nx];
            int image1 = Triples[indx1].image;
            int image2 = Triples[indx2].image;
	    if( image1 != image2 && image1 > 0 && image2 > 0 ) { // found a boundary
		Ntrans++;
		int i1 = relevant_images[image1-1];  // -1 since 0 is reserved for 'no image maps here'
		int i2 = relevant_images[image2-1];  // point is on boundary of image i1 and i2
		//i2 = i1;  // temp, just to make sure we get 0,0
		Point p1(tx,ty);
		CommonPoint(images, p1,i1, Triples[indx1].patch, out_layer);
		}
	    }
        if( Ntrans > 100 ) {
	    printf("Too many transitions (%d) in map!\n", Ntrans);
            if( !Debug ) {
		printf("Writing map as 'test.png'\n");
        	Raster16ToPng16("test.png", &imap[0], nx, ny);
		}
            return 42;
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
        double move = NewG[i].Dist( gvtx[i] );
	printf("Image %d, pt %d, was (%f,%f), now (%f,%f), moved %f %s\n",
	 i/(N+1), i%(N+1), gvtx[i].x, gvtx[i].y, NewG[i].x, NewG[i].y, move, move>MoveThresh ? "<-----" : "");
        if( move > MoveThresh ) {
	    NewG[i] = gvtx[i];
	    NTooFar++;
	    }
	}
    // Now, for every patch of every image, compute a vector of N transforms, one corresponding
    // to each original triangle.
    for(int i=0; i<images.size(); i++) {
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
    for(double r=0.0; r <= max(center.x,center.y)+0.0001; r += delta_image_space*0.9) {
	vector<Point> pts;
	PointsInRing(pts, center, r, w, h);
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
                if( patch == 0)  // || patch >= images[i].tf.size() || images[i].tf[patch].det() == 0.0 )
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
		uint32 nn = ix + uint32(iy)*uint32(nx);   // index into array
		if( imap[nn] == 0 ) {
                    if( ix == 8557 && iy == 431 )
		        printf("Setting pixel %d %d\n", ix, iy);
		    imap[nn] = k+1;                        // this pixel will be set by the kth relevant picture
		    PatchFrom[nn] =  patch;  // patch 'patch', sector 'sector'
                    SectorFrom[nn]= sector;
                    }
		}
	    }
	}
    // Now compress the map/patch/sector array.  We assume less than 65536 combinations are used, so we can express
    // this as a 16 bit TIF.  The combinations start at 1, and run sequentially.
    temp.clear();               // clear the mapping file
    Triples.resize(1);          // and the triples; will fill up as they are found; first entry remains (0,0,0)
    for(int i=0; i<imap.size(); i++) {
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
        imap[i] = it->second;
	}
    // write the new map file
    sprintf(fname,"map.%d.png", out_layer);
    printf("write Pixel mapping file %s\n", fname);
    Raster16ToPng16(fname, &imap[0], nx, ny);
    printf("Draw the new image, with warp\n");
    // Now redraw the image, using the new map.  We'll write it into 'before', which is somewhat
    // confusing.
    memset( &before[0], 0, nx * ny * sizeof(uint8) );

    // also, create a super-pixel map.  This needs to be 32 bit, even though that's a humongous
    // array, since there are typically 10K superpixels per image, and 9x9 arrays are typical,
    // so there are way more than 2^16 superpixel IDs.
    vector<uint32> spmap(nx*ny,0);
    for(int x=0; x<nx; x++) {
	for(int y=0; y < ny; y++) {
	    uint16 indx = imap[x + nx*y];
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

		double pix = InterpolatePixel( p.x, p.y, images[img].raster, w );
		int ix = (int)p.x;
		int iy = (int)p.y;
		before[x + nx*y] = ROUND(pix);
                // for the super-pixel map we want nearest, not interoplation.
                if( p.x - ix >= 0.5 )
		    ix++;
                if( p.y - iy >= 0.5 )
		    iy++;
                int px = images[img].spmap[ix + w*iy];
                if( px != 0 )  // 0 valued pixels are unassigned, and not translated
                    px += images[img].spbase;
                spmap[x + nx*y] = px;
                if( spmap[x + nx*y] == 65536 ) {  // Just debugging
                    printf("w=%d h=%d nx=%d ny=%d patch=%d sector=%d\n", w, h, nx, ny, patch, sector);
		    printf("Set to %d.  x=%d y=%d ix=%d iy=%d img=%d images[img].spbase=%d images[img].spmap[ix + w*iy] = %d\n",
                     px, x, y, ix, iy, img, images[img].spbase, images[img].spmap[ix + w*iy] );
                    return 42;
		    }
                if( Debug && px == 0 )
		    before[x + nx*y] = 255;
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
		    if( 0 <= ix && ix < nx && 0 <= iy && iy < ny )
			before[ix + nx*iy] = 255;
		    }
		}
	    }
        }
    sprintf(fname,"after.%d.tif", out_layer);
    //Raster8ToTif8( fname, &(before[0]), nx, ny );
    sprintf(fname,"after.%d.png", out_layer);
    Raster8ToPng8(fname, &before[0], nx, ny);

    // ------------------------------------------------ write the mapping text file  --------------------
    sprintf( fname, "mapping.%d.txt", out_layer );
    FILE	*ftxt = FileOpenOrDie( fname, "w" );

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
	TAffine t = images[image].sinvs[patch][sector];  // transform from image to global
        fprintf( ftxt, "TRANS %d %d: ", k, from-1 );
        t.TPrint( ftxt );
        }
    fclose(ftxt);


    // Write the super pixel map file.  First make sure it has all consecutive numbers, in the interest of
    // efficiency in Raveler.
    int SPmax = 0;
    vector<int> FinalMapping;
    RemapSuperPixelsOneImage(&(spmap[0]), nx, ny, SPmax, FinalMapping);
    printf("SPmax now %d\n", SPmax);
    sprintf(fname,"sp.%d.png", out_layer);
    for(int i=0; i<nx*ny; i++)  // set the transparecy
	spmap[i] = spmap[i] | 0xFF000000;
    Raster32ToPngRGBA(fname, &spmap[0], nx, ny);
    bool experiment = false;
    if( experiment ) {
        vector<uint32> copy(nx*ny,0);
	for(int i=0; i<nx*ny; i++)
	    copy[i] = 0xFF000000 + spmap[i];
        printf("try to write the small array as a .png file\n"); fflush(stdout);
        Raster32ToPngRGBA("try32bita.png", &copy[0], nx, ny);
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
        Raster32ToPngRGBA("try32bitb.png", &big[0], bigx, bigy);
        printf("done writing the huge array as a .png file\n"); fflush(stdout);
        }

    // write out the mapping files.  For every image that has a non-blank sp map;
    sprintf( fname, "map-sp.%d.txt", out_layer );

    FILE	*fmap = FileOpenOrDie( fname, "w" );

    for(int j=0; j<relevant_images.size(); j++) {
	int i = relevant_images[j];
	if( images[i].spmap != NULL && images[i].spmap != BlankSPMap ) { // one was defined
            fprintf(fmap,"Map '%s'\n", images[i].spname);
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

    printf("WARNINGS:  %d hit the edge, %d had bad correlation\n", nwarn_edge_interp, nwarn_bad_corr);
    }
return 0;
}


