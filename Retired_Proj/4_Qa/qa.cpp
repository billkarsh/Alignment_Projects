

#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"

#include	<stdlib.h>






#if 0	// not used here, but interesting ideas on fold masking

// Here we convert an input image to a map where 0 = on fold, 1=region 1, 2=region 2, etc.
// Normally the pipeline does this, but we do it here for standalone testing.
void ImageToFoldMap(Picture &pic, uint8* FoldMask)  //pipeline does this in normal operation
{

// Here, define the connected region(s) on the above layer.
int w = pic.w;
int h = pic.h;
uint8 *raster = pic.raster;
int npixels = w * h;

// First, find the mean and standard deviation of the 'real' non-saturated pixels.
// Ignore any that are too close to saturated light or dark.
int SAT=3;  // how close to the boundary do you need to be to be considered 'saturated'
MeanStd m;
for(int i=0; i<npixels; i++) {
    uint8 pix = raster[i];
    if( pix >= SAT && pix <= 255-SAT )
	m.Element(raster[i]);
    }
printf("%d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);

if(m.HowMany()/double(npixels) < 0.9) {
    printf("Saturated image!  Retrying with SAT=1\n");
    SAT = 1;
    m.Reset();
    for(int i=0; i<npixels; i++) {
	uint8 pix = raster[i];
	if( pix >= SAT && pix <= 255-SAT )
	    m.Element(raster[i]);
	}
    printf("%d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);
    }

double mean, std;
m.Stats(mean, std);  // find existing statistics
printf("Of the above image points, mean= %f and std dev = %f\n", mean, std);

// If pure white is outside the range of practical values, then it too
// should be ignored, as it is not useful info.  But if 255 is within the
// 'practical' range, we don't want to cut these pixels out
if( mean + 2.5*std < 255 ) {
    printf("Removing white pixels\n");
    for(int i=0; i<npixels; i++) {
	uint8 pix = raster[i];
	if (pix > 255-SAT)           // too bright is just as bad as too dark - not useful info
	    raster[i] = 0;
	}
    }

// Now find the connected regions.  We want to ignore very black pixels as folds, and those
// near them (distance = D) as unreliable, and might connect regions which should be disconnected.
// However, if the image as a whole is too black, this will result in lots of disconnected
// regions, as there will be many black pixels and many near them.  So reset the parameters in
// these cases.
int D=10;
int nbig = 0;
vector<ConnRegion> cr;
double thresh = 4.0;          // we would like 4 std dev, if possible
if( mean - thresh*std < 1 ) {  // but if not, pick a value that is feasible
    thresh = mean/std*0.95;   // set to 95% of the way to 0
    printf("Forced to reduce threshold to %f std dev.\n", thresh);
    if( thresh < 2.0 ) {       // we'll get too many black, and fragment the area
	thresh = (mean - 0.5)/std;  // set threshold to a pixel value of 0.5 (on scale of 255)
	printf("Desperate measures.  Changing threshold to %f, so only v=0 pixels are out.\n", thresh);
        printf("Also disabling black pixel expansion\n");
	D = 0;
	}
    }
const char *p = getenv("FoldMaskThreshold");
if( p != NULL ) {
   thresh = atof(p);
   printf("Environment variable over-ride of threshold to %f\n", thresh);
   }

// pixels to vector
vector<double> v(npixels);
for(int i=0; i<npixels; i++) {
    int y = i / w;
    int x = i - w * y;   // pixels next to background pixels should also be black
    int pix = raster[i];
    v[i] = (pix-mean)/std;
    }

// If there are chunks with no contrast, get rid of them.
double dx = (w-128)/64.0;  // idea is to tile with overlapping 128x128 squares
double dy = (h-128)/64.0;
vector<int> zap_me;        // save list of pixels to be zapped, since cannot zap on the fly without
                           // munging the overlapping squares
for(double y=0; y< (h-127); y += dy) {
    int ymin = min(int(y), h-128);  // just in case
    int ymax = ymin+127;
    for(double x=0; x< (w-127); x += dx) {
	int xmin = min(int(x), w-128);
	int xmax = xmin+127;
        vector<double>local(128*128);
        printf("xmin, ymin = (%d %d)\n", xmin, ymin);
        for(int ix=xmin; ix <= xmax; ix++)
	    for(int iy=ymin; iy <= ymax; iy++)
		local[(ix-xmin) + 128*(iy-ymin)] = v[ix + w*iy];
        if( IsLowContrast(local, std) ) {
	    for(int ix=xmin; ix <= xmax; ix++)
		for(int iy=ymin; iy <= ymax; iy++)
		    zap_me.push_back(ix + w*iy);
	    }
	}
    }
for(int i=0; i<zap_me.size(); i++)
    v[zap_me[i]] = -thresh - 100.0;  // sure to be bad

// for those that should not be included, remove all pixels within distance D as well.  This will
// remove tiny filaments connected big regions together.
vector<int> remove;  // indices of points to be removed.
for(int i=0; i<npixels; i++)
    if( v[i] < -thresh )
	remove.push_back(i);
for(int ii=0; ii<remove.size(); ii++) {
    int i = remove[ii];
    int y = i / w;
    int x = i - w * y;
    int x0 = max(x-D,0);
    int x1 = min(x+D, w-1);
    int y0 = max(y-D,0);
    int y1 = min(y+D,h-1);
    for(int xx = x0; xx <= x1; xx++) {
	for(int yy = y0; yy <= y1; yy++) {
	    v[xx + w*yy] = -thresh - 1.0;   // set to a value that is more than enough black
            }
        }
    }

// Now find the connected regions

	for( int i = 0; i < npixels; ++i ) {

		if( v[i] > -thresh ) {

			ConnRegion	c;
			int			npts;

			npts = Propagate( c.pts, v, w, h, i,
					-thresh, -thresh - 1.0 );

			printf(
			"ImageToFoldMap: ConnRegion with %d pixels.\n", npts );

			if( npts > 90000 ) {	// want 100k, but we shrank it
				cr.push_back( c );
				++nbig;
			}
		}
	}

// Final accounting

   SetBoundsAndColors( cr, FoldMask, w, h );
}

#endif


// is it a legal region?
bool lr(int sx, int sy, void *arg)
{
	//printf( "lr: sx, sy %d %d\n", sx, sy );
	return sx == 21 && sy == 21;
}


// is the pixel count legal?
// We want at least half the pixels
// in each sample to be non-zero
bool lc(int c1, int c2, void *arg)
{
	//printf( "lc: c1,c2 %d %d\n", c1, c2 );
	return c1 > 220 && c2 > 220;
	return true;
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main(int argc, char **argv)
{
if( argc < 4 ) {
    printf("Usage: qa <tif-file> <tif-file> <out-file> [option]\n");
    exit( 42 );
    }

uint32 w, h;
uint8* a = Raster8FromTif( argv[1], w, h );
uint32 ww,hh;
uint8* b = Raster8FromTif( argv[2], ww, hh );
if( w != ww || h != hh ) {
    printf("Images are different sizes w,h=%d,%d w,h=%d %d\n", w, h, ww, hh);
    return 42;
    }
const int P=105;  // patch size is PxP
const int L=21;   // size of region to look for
int ow = w/P;
int oh = h/P;
int oxi=0, oyi;
vector<uint32>	Result( ow * oh, 0xFF000000 );	// set alpha chan

for(double fx = 0.0; fx < w-P/2; fx += double(w)/ow) {
    int x = int(fx);
    printf("x=%d\n", x);
    oyi = 0;
    for(double fy = 0.0; fy < h-P/2; fy += double(h)/oh) {
        int y = int(fy);
	// now correlate a PxP square at this location.  Normalize both
        MeanStd sa, sb;  // stats for both
        vector<Point> pa;
	vector<double> aa;
	int n0 = x + w*y;
	for(int dx=0; dx < P; dx++) {
            int n1 = n0+dx;
	    for(int dy=0; dy < P; dy++) {
		aa.push_back(a[n1]);
                pa.push_back(Point(dx,dy));
                sa.Element(a[n1]);
		n1 += w;
		}
	    }
        // now find a smaller version from the center of the patch on B.
        vector<Point>pb;
        vector<double> bb;
        int offset = (P-1)/2 - (L-1)/2;
        n0 = x+offset + w*(y+offset);
	for(int dx=0; dx < L; dx++) {
            int n1 = n0+dx;
	    for(int dy=0; dy < L; dy++) {
		bb.push_back(b[n1]);
                pb.push_back(Point(dx,dy));
                sb.Element(b[n1]);
		n1 += w;
		}
	    }
        // get the mean and std deviation of each.
        double meana, meanb, stda, stdb;
        sa.Stats(meana, stda);
        sb.Stats(meanb, stdb);
        //printf("m,s=%f %f, m,s= %f %f\n", meana, stda, meanb, stdb);
        double corr = 0.0;
        double dx = 0.0,dy= 0.0;
        if( stda > 0.001 && stdb > 0.001 ) { // then measure it
            vector<CD> ftc;  // fourier transform cache
            corr = CorrPatches(
					stdout, false, dx, dy,
					pb, bb, pa, aa, 0, 0, 4000,
					lr, NULL, lc, NULL, ftc );
            dx -= offset;
            dy -= offset;
            //printf("at x,y=%d,%d ---> dx, dy= %f, %f, corr %f\n\n", x, y, dx, dy, corr);
	    }
        //printf("x=%d, y=%d, corr=%f\n", x, y, corr);
        double inten =  255.0*max(0.0, min(corr, 1.0));
        //now distribute between the green and red channels
        double alpha = fmin(sqrt(dx*dx + dy*dy)/40.0, 1.0);
        int green = int((1-alpha)*inten);
        int   red = int((  alpha)*inten);
        if( 0.25 < alpha && alpha < 0.75 )
	    printf("alpha %f inten %f red %d green %d\n", alpha, inten, red, green);
        Result[oxi + oyi*ow] |= ((green << 8) | red);
        oyi++;
	}
    oxi++;
    }

// Write it out
//Raster8ToTif8( argv[3], FoldMask, ow, oh );
Raster32ToTifRGBA( argv[3], &Result[0], ow, oh );

return 0;
}


