

#include	"FoldMask.h"

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"CTemplate.h"
#include	"CCorrImages.h"
#include	"CCorrCand.h"
#include	"Cffmap.h"

#include	"ls_svd.h"
#include	"tinyxml.h"

#include	<string.h>

#include	<queue>
using namespace std;


bool TrOnly = false;      // Use translation only?

int OVERLAP = 75;  //expect about this much overlap, in pixels
int RADIUS = 256;  // correction for alignment should be found within this distance
double THRESHOLD = 0.25;  // lowest correlation considered a match


class Picture : public PicBase {

public:
	double	dist_from_center;	// distance from center if exists

public:
	Picture()	{dist_from_center = 0.0;};

	bool operator < (const Picture &rhs) const
		{return dist_from_center < rhs.dist_from_center;};
};






// Improve the correlation by tweaking the location of the control points.
// returns the best correlation obtained.
double ImproveControlPts(
 vector<Point> &cpts,              // the control points
 vector<vector<double> > &lambda,  // pixels coords expressed as linear combo of control points
 vector<double> &spv,              // the values at these pixels
 vector<double> &image2,           // the image we are mapping to
 int w,                            // and its width
 int h,							   // and its length
 FILE *flog)                       // file for logging
{
printf("Contains %d pixels\n", spv.size() );
Normalize(spv);

double best_so_far = 0.0;

vector<double> Tpixels;          // Computed target pixel values go here
vector<Point> cderivs;           // and control point derivatives go here
GradDescStep( Tpixels, cderivs, cpts, lambda, spv, image2, w, h );
int nnz;
double corr = CorrVectors(NULL, spv, Tpixels, &nnz);
const double cthresh = 0.20;
//now do an initial check for plausibility
if( corr < cthresh ) {  // if unexpected, dump out data passed in for debugging
    printf("Correlation should not be less than %f at the start of iterations\n", cthresh);
    printf("Control points are:");
    for(int i=0; i<cpts.size(); i++)
        printf("(%.3f %.3f) ", cpts[i].x, cpts[i].y);
    printf("corr=%f\n", corr);
    for(int i=0; i<lambda.size(); i++) {
	printf("---i=%d\n",i);
	for(int j=0; j<lambda[i].size(); j++)
	    printf("%.4f ", lambda[i][j]);
        printf("\n spv=%f Tpixels=%f\n", spv[i], Tpixels[i]);
        }
    exit( 42 );
    }

// Now try to tweak the control points for a good match
for(double step=10; step > 0.05; ) {
    printf("\n");
    printf("Control points are:");
    for(int i=0; i<cpts.size(); i++)
        printf("(%f %f) ", cpts[i].x, cpts[i].y);
    printf("corr=%f, step is %f\n", corr, step);
    // compute a downhill vector of length 'step'.
    double sum2 = 0.0;  // sum of squares
    for(int i=0; i<cderivs.size(); i++) {
        Point p = cderivs[i];
	sum2 += p.x*p.x + p.y*p.y;
        }
    double mag = sqrt(sum2);  // length of the vector
    // find the new (and hopefully improved) control points
    vector<Point> newpts = cpts;
    for(int i=0; i<cpts.size(); i++) {
       newpts[i].x += step*cderivs[i].x/mag;
       newpts[i].y += step*cderivs[i].y/mag;
       }
    // is the new spot better?
    vector<Point>  newderivs;
    GradDescStep( Tpixels, newderivs, newpts, lambda, spv, image2, w, h );
    double c = CorrVectors(NULL, spv, Tpixels, &nnz);
    if( c > corr ) { // new point is better
	corr = c;
        cpts = newpts;
        cderivs = newderivs;
        }
    else
	step = step / 2;
    }
return corr;
}

void TryNewOptimizer(vector<Point> &plist, vector<double> &spv, vector<double> &image2, TForm &t, FILE *flog)
{
// compute the bounding box
printf("\n---------- Try new optimizer on %d points----------------\n", plist.size());
DBox	B;
BBoxFromPoints( B, plist );
printf("region size is [%f %f] in x, [%f %f] in y\n", B.L, B.R, B.B, B.T);

// create 3 control points
vector<Point> cpts;
cpts.push_back(Point(B.L, B.B));
cpts.push_back(Point(B.R, B.B));
cpts.push_back(Point((B.L+B.R)/2, B.T));
printf("Control points are (%f %f) (%f %f) (%f %f)\n", cpts[0].x, cpts[0].y, cpts[1].x, cpts[1].y, cpts[2].x, cpts[2].y);

// find each point as a linear combination of control points
double a[3][3];
a[0][0] = cpts[0].x; a[0][1] = cpts[1].x; a[0][2] = cpts[2].x;
a[1][0] = cpts[0].y; a[1][1] = cpts[1].y; a[1][2] = cpts[2].y;
a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
double inv[3][3];
Print3x3Matrix( stdout, a );
Invert3x3Matrix( inv, a );
Print3x3Matrix( stdout, inv );
vector<vector<double> > lambda;
for(int j=0; j<plist.size(); j++) {
    //printf(" Point is (%f %f)\n", plist[j].x, plist[j].y);
    vector<double> lam;
    lam.push_back(inv[0][0]*plist[j].x + inv[0][1]*plist[j].y + inv[0][2]*1.0);
    lam.push_back(inv[1][0]*plist[j].x + inv[1][1]*plist[j].y + inv[1][2]*1.0);
    lam.push_back(inv[2][0]*plist[j].x + inv[2][1]*plist[j].y + inv[2][2]*1.0);
    //printf(" lambdas are %f %f %f\n", lam[0], lam[1], lam[2]);
    lambda.push_back(lam);
    }

// Transform the control points to the target frame
vector<Point> orig = cpts;
t.Transform( cpts );

// call the optimizer
ImproveControlPts(cpts, lambda, spv, image2, 4096, 4096, flog);

// Now, find a transformation that maps ORIG into the new cpts
// first, create a transform that maps a unit right triangle to the original pts
TForm o(orig[1].x-orig[0].x, orig[2].x-orig[0].x, orig[0].x,
        orig[1].y-orig[0].y, orig[2].y-orig[0].y, orig[0].y);
// now one that maps the final optimized control points to the unit right triangle
TForm c(cpts[1].x-cpts[0].x, cpts[2].x-cpts[0].x, cpts[0].x,
        cpts[1].y-cpts[0].y, cpts[2].y-cpts[0].y, cpts[0].y);
// now, to get from the original to the final, apply o^-1, then c;
TForm oi;
InvertTrans(oi, o);
//TForm temp;
MultiplyTrans(t, c, oi);
t.PrintTransform();
}

// Improve the correlation, if possible, by tweaking the transform.  pts are the points
// in the original, and pts are their values.  image is a 4Kx4K matrix of doubles, already
// normalized.  dx and dy are the initial estimates of how to map the original points into
// the array image2.
// returns the best correlation obtained.
double ImproveCorrelation(vector<Point> &Plist, vector<double> &spv, vector<double>image2,
 double dx, double dy, TForm &t, FILE *flog)
{
printf("Contains %d pixels\n", Plist.size() );
Normalize(spv);
if (dx != BIG)  // if dx == BIG, start improving from the transform we have
    t = TForm(1.0, 0.0, dx, 0.0, 1.0, dy);  // otherwise, create a transform with just dx, dy

double best_so_far = 0.0;

TryNewOptimizer(Plist, spv, image2, t, flog);

/*
// Now try to find a transform with a good match.
for(double step=1; step > 0.05; ) {
    printf("\n");
    double dvx, dvy;
    vector<Point> Tpoints = Plist;   //   Start with locs in source
    t.Transform( Tpoints );            // Transform to locations in target
    vector<double> Tpixels;          // Get the pixel values
    ValuesFromImageAndPoints( Tpixels, dvx, dvy, image2, 4096, Tpoints, spv );
    if( fabs(dvx) < 1.0E-6 && fabs(dvy) < 1.0E-6 ) {
	printf("No derivatives!\n");
	fprintf(flog,"No derivatives!\n");
        exit( 42 );
	}
    //Normalize(Tpixels);
    Point cog;  // Center of gravity
    cog = FindCOG(Plist, Tpixels);
    //printf("original COG is %f %f\n", cog.x, cog.y);
    int nnz;
    double start = CorrVectors(NULL, spv, Tpixels, &nnz);
    printf(" %f %f %f %f %f %f: correlation %f\n",
     t.t[0], t.t[1], t.t[2], t.t[3], t.t[4], t.t[5], start);
    best_so_far = start;
    TForm tbest;  // best transform found
    int bdir = -1;  // flag to tell what the best direction is

    for(int rot = -1; rot<2; rot += 2) {  // try two rotations
	TForm t2 = t;
	Point cogt(cog.x, cog.y);
	t2.Transform( cogt ); // transform the center of gravity
	t2.RotateAround( cog, cogt, rot*step*PI/180 ); // +- step degrees
	vector<Point> Tpoints = Plist;   //   Start with locs in source
	t2.Transform( Tpoints );            // Transform to locations in target
        double junk1, junk2;
	ValuesFromImageAndPoints( Tpixels, junk1, junk2, image2, 4096, Tpoints, spv );
	//Normalize(Tpixels);
	double nc = CorrVectors(NULL, spv, Tpixels);
	printf(" rotated case %d: correlation is %f\n", rot, nc);
	if (nc > best_so_far){
	    best_so_far = nc;
	    bdir = 10;
            tbest = t2;
	    }
	}

    // make steps of this size until no more improvement
    double best_rot = best_so_far;  // best result from rotations
    double prev;
    double nc = start;     // new correlation
    TForm t2 = t;
    do {
        prev = nc;
	double norm = sqrt(dvx*dvx + dvy*dvy);
	double sx = step*dvx/norm;  // direction sin()
	double sy = step*dvy/norm;
	printf("step is %f %f\n", sx, sy);
	t2.AddXY( sx, sy );
	Tpoints = Plist;   //   Start with locs in source
	t2.Transform( Tpoints );            // Transform to locations in target
        double pdvx=dvx, pdvy=dvy; // just for comparing real with prediction
	ValuesFromImageAndPoints( Tpixels, dvx, dvy, image2, 4096, Tpoints, spv );
	//Normalize(Tpixels);
	nc = CorrVectors(NULL, spv, Tpixels, &nnz);
	printf(" predict %f, real correlation is %f\n", prev+sx*pdvx/nnz+sy*pdvy/nnz, nc);
	if (nc > best_so_far){
	    best_so_far = nc;
            tbest = t2;
	    bdir = 1;
	    }
        }
    while(nc > prev && nc > best_rot);

    // Now tried two rotations and a direction; pick the best if any were better
    if( bdir >= 0 ) { //we found a better transform
         t = tbest;
	 }
    else // nothing was better; reduce step size
	step = step/2;
    }
*/
// Now t is the best transform we can find.
printf(
"Best transform is %9.4f %9.4f %10.2f\n"
"                  %9.4f %9.4f %10.2f\n",
t.t[0], t.t[1], t.t[2], t.t[3], t.t[4], t.t[5] );

return best_so_far;
}


// Align the jth entry to the ith entry.
void AlignAPair(vector<Picture> &vp, int i, int j, FILE *flog, int pass)
{
// Compute the mean and std deviation of the first using only the 'real' pixels -
// those that are >0, and hence not part of the background.
printf("Aligning images %d and %d\n", i, j);
vector<double> v2;
int npixels = vp[i].w*vp[i].h;
for(int k=0; k<npixels; k++)
    if( vp[i].raster[k] > 0 )
	v2.push_back(vp[i].raster[k]);
printf("Image i: %d real pixels, %f percent\n", v2.size(), v2.size()*100.0/npixels);
double mean2, std2;
Stats(v2, mean2, std2);  // find existing statistics
printf("Of the target image,  mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy.  Background pixels map to zero
vector<double> image2(4096*4096, 0.0);
for(int k=0; k<npixels; k++) {
    int y = k / vp[i].w;
    int x = k - vp[i].w * y;
    double pix = vp[i].raster[k];
    if (pix == 0)    // background pixels
	pix = mean2;  // will be set to the mean value (0)
    image2[x + 4096*y] = (pix-mean2)/std2;
    }

// Now, for the other image, do this
// Find the points in image that map onto the first picture
vector<Point> pts;
vector<double> vals;
for(int k=0; k<vp[j].w*vp[j].h; k++) {
    int y = k / vp[j].w;
    int x = k - vp[j].w * y;
    Point pt(x,y);            // Coordinates in picture j
    vp[j].tr.Transform( pt );  // now in global space
    vp[i].Inverse.Transform( pt ); // Into the space of picture i
    if( pt.x >= 0.0 && pt.x < vp[i].w && pt.y > 0.0 && pt.y < vp[i].h ) {
	//printf("pt: x,y=%f %f\n", pt.x, pt.y);
	pts.push_back(pt);
	vals.push_back(vp[j].raster[k]);
	}
    }
printf("There are %d overlap pixels\n", pts.size() );
if( pts.size() < 10000 ) {
    printf("Not enough pixels.\n");
    fprintf(flog,"Not enough pixels.\n");
    exit( 42 );
    }
Normalize(vals);
// assume radius no worse than 400(?)
double	dx, dy;
double	c = CorrPatchToImage( dx, dy, pts, vals, image2, 0, 0, RADIUS, true );
if( pass == 2 && c < 0.25 ) {
    printf("Dubious alignment improvement: %d pts, corr=%f\n", pts.size(), c);
    printf("file 1: %s\n", vp[i].fname.c_str() );
    printf("file 2: %s\n", vp[j].fname.c_str() );
    fprintf(flog,"Dubious alignment improvement: %d pts, corr=%f\n", pts.size(), c);
    exit( 42 );
    }
printf(" deltas are %f %f\n", dx, dy);
Point zero(0.0, 0.0);
Point delta(dx, dy);      // what we want in space of image i
vp[i].tr.Transform( delta );        vp[i].tr.Transform( zero );
// make the fix in global space, since this is where the tr[2] and [5] are applied
vp[j].tr.t[2] += (delta.x - zero.x);
vp[j].tr.t[5] += (delta.y - zero.y);
printf("deltas in global image space %f %f\n", delta.x-zero.x, delta.y-zero.y);
InvertTrans(vp[j].Inverse, vp[j].tr);  // Recompute j's inverse
}

// are two images neighbors?  .  If pass 2, just tell if they overlap by at last half in
// one dimension.  But in pass 1, if they are 'close' and 'should' be neighbors, 'fixes' the transform
bool neighbor(vector<Picture> &vp, int i, int j, int pass)
{
    // First do crude computation Transform opposing corners if j into i's space
    Point pll(0.0,0.0);
    vp[j].tr.Transform( pll );			// now in global space
    vp[i].Inverse.Transform( pll );		// Into the space of picture i
    Point pur( vp[j].w-1, vp[j].h-1 );	// upper right
    vp[j].tr.Transform( pur );			// now in global space
    vp[i].Inverse.Transform( pur );		// Into the space of picture i
    double xmin = max(0.0, pll.x);
    double xmax = min(double(vp[i].w-1), pur.x);
    double ymin = max(0.0, pll.y);
    double ymax = min(double(vp[i].h-1), pur.y);
    printf("Compare pass %d, %d to %d: [%f %f] [%f %f]\n", pass, i, j, xmin, xmax, ymin, ymax);
    double half = vp[i].w/2;
    if( pass == 2 )
	return xmax > xmin && ymax > ymin && (xmax - xmin > half || ymax - ymin > half);
    bool status = false; // assume no overlap
    double xshift=0, yshift=0;   // possible shifts to correct
    double area = (xmax-xmin)*(ymax-ymin);
    if( xmax > xmin && ymax > ymin && area > vp[i].w*0.666*OVERLAP ) {  // at least 2/3 the expected overlap
	status = true;  // they do overlap in any case
        if (area < vp[i].w*1.3333*OVERLAP)  // if more than this, trim it back
	    return true;
        // otherwise maybe too big; could need to shift
        if( ymax - ymin > half && xmin == 0 && xmax > OVERLAP )
            xshift = OVERLAP - xmax;
        if( ymax - ymin > half && xmax == vp[i].w-1 && xmin < vp[i].w-OVERLAP )
	    xshift = vp[i].w-OVERLAP - xmin;
        if( xmax - xmin > half && ymin == 0 && ymax > OVERLAP )
            yshift = OVERLAP - ymax;
        if( xmax - xmin > half && ymax == vp[i].h-1 && ymin < vp[i].h-OVERLAP )
	    yshift = vp[i].h-OVERLAP - ymin;
	}
    // Another case is that they should overlap, but don't really.  Tell this if there is
    // at least a half picture overlap in one direction, but a less than half picture gap
    // in the other.
    if( ymax-ymin > half && xmin > vp[i].w-OVERLAP && xmin < vp[i].w + half ) { // to the right
        status = true;
	xshift = vp[i].w - OVERLAP - xmin;
        }
    if( ymax-ymin > half && xmax < OVERLAP && xmax > -half ) { // to the left
        status = true;
	xshift = OVERLAP - xmax;
        }
    if( xmax-xmin > half && ymin > vp[i].h-OVERLAP && ymin < vp[i].h + half ) { // to the top
        status = true;
	yshift = vp[i].h - OVERLAP - ymin;
        }
    if( xmax-xmin > half && ymax < OVERLAP && ymax > -half ) { // to the bottom
        status = true;
	yshift = OVERLAP - ymax;
        }
    if( status ) {
        printf("Crude shift: dx=%f dy=%f\n", xshift, yshift);
        Point zero(0.0,0.0);
        Point vec(xshift, yshift);
        vp[i].tr.Transform( zero ); vp[i].tr.Transform( vec );
        vec.x -= zero.x; vec.y -= zero.y;  // vector we want in global space
        vp[j].tr.AddXY( vec.x, vec.y );
        }
     return status;
}
// if there is more than one picture that overlaps a tile, improve the alignment
// by considering them pairwise.  Start with the first (biggest) one, do the
// north-east-south-west neighbors (diagonal does not have enough overlap). Then
// do the neighbors of those already aligned, etc.
void ImproveAlignment(vector<Picture> &vp, FILE *flog, int pass)
{
vector<bool> aligned(vp.size(), false);
aligned[0] = true;  // First one, the biggest one, is aligned by definition
int Nnot = vp.size()-1;
for(; Nnot > 0; Nnot--) {
    // Find an unaligned entry with an aligned neighbor (there should be at least one)
    int i,j;
    bool found = false;
    for(i=0; i<vp.size() && !found; i++) {
        for(j=0; j<vp.size() && !found; j++)
	    found = aligned[i] && (!aligned[j]) && neighbor(vp,i,j, pass);
        }
    if( !found ) {
        printf("No aligned-unaligned neighboring pair\n");
        fprintf(flog,"No aligned-unaligned neighboring pair\n");
	exit( 42 );
	}
    i--; j--;  //loops incremented once more after finding
    AlignAPair(vp, i, j, flog, pass);
    aligned[j] = true;
    }
}



// Find the correlation between a patch of j and the whole picture i.
// Picture i has already been fft'd.
Point SubCorrelate(vector<Picture> &vp, int i, int j, int xmin, int ymin, int xmax, int ymax)
{
int N = 4096;
int M = N*(N/2+1);  // number of complex numbers in FFT of 2D real

int w = vp[j].w;
int h = vp[j].h;
int np = w * h; // number of pixels
// create a vector of all the pixels in the piece, so we can normalize their values.
vector<double> px;
for( int k = 0; k < np; ++k ) {
    int y = k/w;
    int x = k-w*y;
    if( x >= xmin && x <= xmax && y >= ymin && y <= ymax )
	px.push_back(vp[j].raster[k]);
    }
double avg, std;
Stats(px, avg, std);  // we now have the mean and standard deviation
printf("sub-frame %d, mean %f, std %f\n", i, avg, std);

// Now make the frame of normalized values.
// Make it in a nx by ny buffer so we can FFT it.

vector<double>	fr( N * N, 0.0 );
vector<CD>		temp( M );
vector<CD>		ff;

for(int k=0; k<np; k++) {
    int y = k/w;
    int x = k-w*y;
    if( x >= xmin && x <= xmax && y >= ymin && y <= ymax )
	fr[x+N*y] = (vp[j].raster[k] - avg)/std;
    }

// Now fft it
FFT_2D( ff, fr, N, N );

// Now multiply image i by the conjugate of the sub-image j
for(int k=0; k<M; k++)
    temp[k] = vp[i].fft_of_frame[k] * conj(ff[k]);

// reverse the FFT to find the lags
IFT_2D( fr, temp, N, N );

// Now find the maximum value
double biggest = -1.0E30;
int bigx = -1;
int bigy = 0;
double norm = (xmax-xmin+1)*(ymax-ymin+1)*double(N)*double(N);
for(int k=0; k<N*N; k++) {
    int y = k / N;
    int x = k - N * y;
    if (x >= N/2) x = x-N;
    if (y >= N/2) y = y-N;
    if( fr[k]/norm > biggest ) {
	biggest = fr[k]/norm;
	bigx = x; bigy = y;
	}
    }

PrintCorLandscape( biggest, bigx, bigy, 0, 0, 0,
	6, 2, &fr[0], N, N, norm, stdout );

printf( "Sub-correlate: Maximum %f at (%d, %d).\n",
	biggest, bigx, bigy );

Point	pt( bigx, bigy );
ParabPeakFFT( pt.x, pt.y, 1, &fr[0], N, N );

printf( "Sub-correlate: Final at (%f, %f).\n",
	pt.x, pt.y );

return pt;
}


// The vector i1-i2 in picture i should match the vector j1->j2 in j.  Make this happen
void MakeMatch(vector<Picture> &vp, int i, Point i1, Point i2, int j, Point j1, Point j2)
{
printf("Make match (%f %f) to (%f %f) in %d, and (%f %f) to (%f %f) in %d\n",
 i1.x, i1.y, i2.x, i2.y, i,   j1.x, j1.y, j2.x, j2.y, j);

// We will adjust the transform of picture j. Transform the first vector from is space
// to global space.
vp[i].tr.Transform( i1 );
vp[i].tr.Transform( i2 );
printf("Now match (%f %f) to (%f %f) in global space, and (%f %f) to (%f %f) in %d\n",
 i1.x, i1.y, i2.x, i2.y,  j1.x, j1.y, j2.x, j2.y, j);
//first, what's the scale change:
double scale = i1.Dist( i2 ) / j1.Dist( j2 );
// and the angle
double ai = atan2(i2.y-i1.y, i2.x-i1.x);
double aj = atan2(j2.y-j1.y, j2.x-j1.x);
printf("Scale %f, angle %f radians, distance %f\n", scale, ai-aj, i1.Dist( i2 ) );
TForm tf(
	scale*cos(ai-aj),
	scale*(-sin(ai-aj)),
	0.0,
	-tf.t[1],
	tf.t[0],
	0.0 );
// now transform the first point and make sure it ends up in the right place
Point p(j1.x,j1.y);
tf.Transform( p );
tf.t[2] += i1.x - p.x;
tf.t[5] += i1.y - p.y;
// OK, test the new transform
Point p1(j1.x,j1.y);
Point p2(j2.x,j2.y);
tf.Transform( p1 );
tf.Transform( p2 );
printf("Point (%f %f), distance %f\n", p1.x, p1.y, p1.Dist( p2 ) );
printf(" Finally get (%f %f) to (%f %f) in global space, and (%f %f) to (%f %f) in %d\n",
 i1.x, i1.y, i2.x, i2.y,  p1.x, p1.y, p2.x, p2.y, j);
p1.x = p1.y = 0.0;
p2.x = vp[j].w-1; p2.y = vp[j].h-1;
tf.Transform( p1 );
tf.Transform( p2 );
printf(" Closure data: corners map to (%f %f) (%f %f)\n", p1.x, p1.y, p2.x, p2.y);
vp[j].tr = tf;
InvertTrans(vp[j].Inverse, tf);
}

// Point Pi in image i should align with point Pj in image j.  i < j.  Direction 'unconstrained'
// is not tighly constrained by the data, and needs an extra constraint to make sure it does
// not contract to 0.
void WriteEqns6(FILE *feq, vector<Picture> &vp, int i, Point Pi, int j, Point Pj,
 vector<vector<double> > &eqns, vector<double> &rhs)
{
int jj = 6*j-5;
printf("Point (%f %f) in pic %d = point (%f %f) in pic %d\n", Pi.x, Pi.y, i, Pj.x, Pj.y, j);
if( i == 0 ) { // tranformation for 0 is fixed.  Shoot for a fixed target
    vp[i].tr.Transform( Pi );  // the point in global space
    // first write the eqn for x.
    fprintf(feq,"%f X%d  %f X%d  1.0 X%d = %f\n", Pj.x, jj, Pj.y, jj+1, jj+2, Pi.x);
    vector<double> eq1(jj+2, 0.0);
    eq1[jj-1] = Pj.x; eq1[jj] = Pj.y; eq1[jj+1] = 1.0;
    eqns.push_back(eq1); rhs.push_back(Pi.x);
    fprintf(feq,"%f X%d  %f X%d  1.0 X%d = %f\n", Pj.x, jj+3, Pj.y, jj+4, jj+5, Pi.y);
    vector<double> eq2(jj+5, 0.0);
    eq2[jj+2] = Pj.x; eq2[jj+3] = Pj.y; eq2[jj+4] = 1.0;
    eqns.push_back(eq2); rhs.push_back(Pi.y);
    }
else { // the two should map to the same point, or Tr(i) - Tr(j) = 0
    int ii = i*6-5;
    fprintf(feq,"%f X%d  %f X%d  1.0 X%d  %f X%d  %f X%d  -1 X%d = 0\n",
     Pi.x, ii, Pi.y, ii+1, ii+2,  -Pj.x, jj, -Pj.y, jj+1, jj+2);
    vector<double> eq1(jj+2, 0.0);
    eq1[ii-1] =  Pi.x; eq1[ii] =  Pi.y; eq1[ii+1] =  1.0;
    eq1[jj-1] = -Pj.x; eq1[jj] = -Pj.y; eq1[jj+1] = -1.0;
    eqns.push_back(eq1); rhs.push_back(0.0);
    fprintf(feq,"%f X%d  %f X%d  1.0 X%d  %f X%d  %f X%d  -1 X%d = 0\n",
     Pi.x, ii+3, Pi.y, ii+4, ii+5,  -Pj.x, jj+3, -Pj.y, jj+4, jj+5);
    vector<double> eq2(jj+5, 0.0);
    eq2[ii+2] =  Pi.x; eq2[ii+3] =  Pi.y; eq2[ii+4] =  1.0;
    eq2[jj+2] = -Pj.x; eq2[jj+3] = -Pj.y; eq2[jj+4] = -1.0;
    eqns.push_back(eq2); rhs.push_back(0.0);
    }
}
// call this one if we are only doing translation.  Just two variables per picture,
// and all transforms are scale and rotation free
void WriteEqns2(FILE *feq, vector<Picture> &vp, int i, Point Pi, int j, Point Pj,
 vector<vector<double> > &eqns, vector<double> &rhs)
{
int jj = 2*(j-1); // equation number
printf("Point (%f %f) in pic %d = point (%f %f) in pic %d\n", Pi.x, Pi.y, i, Pj.x, Pj.y, j);
if( i == 0 ) { // tranformation for 0 is fixed.  Shoot for a fixed target
    vp[i].tr.Transform( Pi );  // the point in global space
    vector<double> eq1(jj+1, 0.0);
    eq1[jj] = 1.0;
    eqns.push_back(eq1); rhs.push_back(Pi.x - Pj.x);
    vector<double> eq2(jj+2, 0.0);
    eq2[jj+1] = 1.0;
    eqns.push_back(eq2); rhs.push_back(Pi.y - Pj.y);
    }
else { // the two should map to the same point, or Tr(i) - Tr(j) = 0
    int ii = (i-1)*2;
    vector<double> eq1(jj+1, 0.0);
    eq1[ii] = 1.0; eq1[jj] = -1.0;
    eqns.push_back(eq1); rhs.push_back(Pj.x-Pi.x);
    vector<double> eq2(jj+2, 0.0);
    eq2[ii+1] = 1.0; eq2[jj+1] = -1.0;
    eqns.push_back(eq2); rhs.push_back(Pj.y-Pi.y);
    }
}












void FindCorrPoints(vector<Picture> &vp, int i, int j, vector<Point> &cpi, vector<Point> &cpj)
{
// clear the point lists, so if we return early they are corrrectly zero-sized.
cpi.clear(); cpj.clear();
// real to complex generates an array of size N *(n/2+1)
int N = 4096;
int M = N*(N/2+1);  // number of complex numbers in FFT of 2D real

vp[i].MakeFFTExist(i);  // make sure the FFTs exist
vp[j].MakeFFTExist(j);
// Now multiply image i by the conjugate of image j

vector<double>	fr;
vector<CD>		temp( M );

for(int k=0; k<M; k++)
    temp[k] = vp[i].fft_of_frame[k] * conj(vp[j].fft_of_frame[k]);

// reverse the FFT to find the lags
IFT_2D( fr, temp, N, N );

// Now find the maximum value
CorrCand best[4];
double biggest = -1.0E30;
int bigx = -1;
int bigy = 0;
for(int k=0; k<N*N; k++) {
    int y = k / N;
    int x = k - N * y;
    if (x >= N/2) x = x-N;
    if (y >= N/2) y = y-N;
    int xspan = vp[i].w - iabs(x);
    int yspan = vp[i].h - iabs(y);
    double norm;
    if (xspan > 0 && yspan > 0 && xspan*yspan > 40*vp[i].w)  // need at least 40 pixels
	norm = double(xspan)*yspan;  // avoid overflow
    else
	norm = 1.0E30; // really big, so this can't be a winner
    double nfeat = norm / 2500;  // feature is typically 50x50, at least
    norm = norm *(1.0+ 2/sqrt(nfeat)); // add 2 standard deviations
    if( fr[k]/norm > biggest ) {
	biggest = fr[k]/norm;
	bigx = x; bigy = y;
	}
    int xp = x-y;  // rotate coordinate frame by 45 deg ccw (plus scaling, which
    int yp = x+y;  // we don't care about
    int indx = (yp>0)*2 + (xp>0); //0 = west, 1 = south, 2 = north, 3 = east
    if( fr[k]/norm > best[indx].val )
	best[indx] = CorrCand(x,y,fr[k]/norm);
    }
// Look at different directions.  Not used for now.
printf("west %d %d %f, south %d %d %f, north %d %d %f, east %d %d %f\n",
 best[0].x, best[0].y, best[0].val/(N*N), best[1].x, best[1].y, best[1].val/(N*N),
 best[2].x, best[2].y, best[2].val/(N*N), best[3].x, best[3].y, best[3].val/(N*N));

// For interest, write the landscape around the max.
// Also, for a true max, we'd expect sort of a quadratic behavior around the peak.
// Also, compute how much peak sticks above surrounding terrain
double ls = 0.0; // local sum
int ln = 0;
double peak;
vector<double>ring;  // ring of non-weighted correlations bigger (expect none for a real peak)
		     // will contain all those RAD away, plus the center point
const int RAD = 3;
for(int iy=-4; iy <=4; iy += 1) {
    int ay = bigy + iy;
    //printf("biggest %f, y=%d, tx, ty, radius= %d %d %d\n", biggest, ay, tx, ty, radius);
    if( ay < 0 )
	ay += N;
    for(int ix=-4; ix<=4; ix += 1) {
	int ax = bigx + ix;
	if( ax < 0 )
	    ax += N;
	double val = fr[ax + N*ay]/(N*N);
	printf("%8.1f ",val);
	if( (iabs(iy) == RAD && iabs(ix) <= RAD) || (iabs(iy)<=RAD && iabs(ix)==RAD) || (ix == 0 && iy == 0) )
	    ring.push_back(val);
        ls += val; ln++;  // for average
        if (ix == 0 && iy == 0) peak = val;
	}
    printf("\n");
    }
int non = 0;
for(int k=0; k<ring.size(); k++)
    non += (ring[k] > ring[ring.size()/2]);  // none should be bigger than middle element

// Biggest correlation we could expect is w (or l) * overlap
// N^2 term is since the FFTs are not normalized
printf("Maximum correlation of %f at [%d,%d], < %d of %d in ring, peak/avg = %f\n",
 biggest/(N*N), bigx, bigy, non, ring.size(), peak/(ls/ln));
// find the region of overlap in frame j's coordinates
int xmin = max(0, -bigx);
int xmax = min(vp[j].w-1, vp[i].w-1-bigx);
int ymin = max(0, -bigy);
int ymax = min(vp[j].h-1, vp[i].h-1-bigy);
printf("overlap region is x=[%d %d] y=[%d %d]\n", xmin, xmax, ymin, ymax);
Point Pcenter1j, off1, Pcenter1i, Pcenter2j, off2, Pcenter2i;
// create the template from the overlap portion in frame i
Template t(vp[i], i, xmin+bigx, xmax+bigx, ymin+bigy, ymax+bigy);
if (xmax-xmin < vp[i].w/2 && ymax-ymin < vp[i].h/2)  // we are looking for full edges
    return;
if (biggest/(N*N) < 0.04)                             // we need a certain amount of quality
    return;
//if (biggest/(N*N) < 0.25 && non > 0)  // not a great correlation
    //return;
//if (biggest/(N*N) < 0.30 && non > 1)  // not a great correlation
    //return;
int nbad = 0;
if( xmax-xmin > ymax-ymin ) {  // cut with a vertical line
    int xmid = (xmin+xmax)/2;
    for(int q=0; q<4; q++) {
	int x1 = xmin + int(double(q)/4*(xmax-xmin));
	int x2 = xmin + int(double(q+1)/4*(xmax-xmin));
	printf("\n");
	printf("---Small area test--- x=[%d %d] y=[%d %d]\n", x1, x2, ymin, ymax);
	Point out = t.Match(vp[j], j, x1+1, ymin+1, x2-1, ymax-1);
	double d = sqrt(pow(out.x-bigx,2.0) + pow(out.y-bigy, 2.0));
	printf("  consistent within %f pixels\n", d);
	if( d < ((q==1 || q==2)? 15.0 : 30.0) ) {
	    Point Pj(double(x1+x2)/2.0, double(ymin+ymax)/2.0); // middle of patch, in j's frame
	    Point Pi(out.x+Pj.x, out.y+Pj.y);
	    cpi.push_back(Pi); cpj.push_back(Pj);
	    }
	else
	    nbad++;
	}
    }
else { // cut with a horizontal line
    int ymid = (ymin+ymax)/2;
    for(int q=0; q<4; q++) {
	int y1 = ymin + int(double(q)/4*(ymax-ymin));
	int y2 = ymin + int(double(q+1)/4*(ymax-ymin));
	printf("\n");
	printf("---Small area test--- x=[%d %d] y=[%d %d]\n", xmin, xmax, y1, y2);
	//Point out = SubCorrelate(vp,i,j,xmin+1,y1+1,xmax-1,y2-1);
	Point out=t.Match(vp[j], j, xmin+1, y1+1, xmax-1, y2-1);
	double d = sqrt(pow(out.x-bigx,2.0) + pow(out.y-bigy, 2.0));
	printf("  consistent within %f pixels\n", d);
	if( d < ((q==1 || q==2)? 15.0 : 30.0) ) {  // limit of 10 for center matches, 20 otherwise
	    Point Pj(double(xmin+xmax)/2.0, double(y1+y2)/2.0); // middle of patch, in j's frame
	    Point Pi(out.x+Pj.x, out.y+Pj.y);
	    cpi.push_back(Pi); cpj.push_back(Pj);
	    }
	else
	    nbad++;
	}
    }
if( nbad > 1 ) {  // get rid of any eqns that were generated
    printf("%d bad points - no correlations generated\n", nbad);
    cpi.clear();
    cpj.clear();
    }
}


// New correlation code
void NewAlign(vector<Picture> &vp, FILE *flog)
{
CCorrImages* CI = CCorrImages::Read( "Corr.xml", flog );

// real to complex generates an array of size N *(n/2+1)
int N = 4096;
int M = N*(N/2+1);  // number of complex numbers in FFT of 2D real

// Step 1 - for each picture, create a frame (just pixels within OVERLAP of edge).
// Then FFT the frame, and keep the FFT

printf("------------------ Experimental correlation -------------------------------\n");
FILE	*feq = FileOpenOrDie( "eqns", "w", flog );

// Here's where we do this internally..
vector<vector<double> > eqns;
vector<double> rhs;
// Now find the correlations among the pairs
for(int i=0; i<vp.size(); i++) {
    for(int j=i+1; j<vp.size(); j++) {
	printf("\nChecking %d to %d\n", i, j);
        printf(" Initial transform %d:", i); vp[i].tr.PrintTransform();
        printf(" Initial transform %d:", j); vp[j].tr.PrintTransform();
        // transform jth image's bounding box to ith image coords, just for fun
        Point bb1(0.0,0.0); Point bb2(vp[j].w, vp[j].h);
	vp[j].tr.Transform( bb1 ); vp[j].tr.Transform( bb2 );
	vp[i].Inverse.Transform( bb1 ); vp[i].Inverse.Transform( bb2 );
	printf("In frame of image %d, image %d is [%f %f] to [%f %f]\n",
	 i, j, bb1.x, bb1.y, bb2.x, bb2.y);

        // first look and see if we have any cached correspondences:
        vector<Point>	cpi, cpj;
        int nc = CI->Find( vp[i].fname, vp[j].fname, cpi, cpj );
        printf(" Got %d saved points\n", nc);
        //for(int k=0; k<cpi.size(); k++) {
	    //Mangle(cpi[k], vp[i].w, vp[i].h);
	    //Mangle(cpj[k], vp[j].w, vp[j].h);
            //}

        //if( nc == 0 ) {
            //FindCorrPoints(vp, i, j, cpi, cpj);
            // are there any correspondence points?  If so add to list
            //if( cpi.size() > 0 )
	        //CI->Add( vp[i].fname, vp[j].fname, cpi, cpj );
            //}
        // generate equations from points, if any
        for(int k=0; k<cpi.size(); k++) {
            if( TrOnly )
	        WriteEqns2(feq, vp, i, cpi[k], j, cpj[k], eqns, rhs);
            else
	        WriteEqns6(feq, vp, i, cpi[k], j, cpj[k], eqns, rhs);
            }
	}

    // for all except the first, add some preference for a square soln (rotation and scaling only)
    if( i != 0 && !TrOnly ) {
        int id = (i-1)*6;
        double strength = 5000.0;  // how strong is this preference?  pixel error vs matrix asymmetry
                                  // these are both true for pure rotation + scaling
        fprintf(feq, "%f X%d  %f X%d = 0\n", strength, id+1, -strength, id+5);  // a11 = a22
        fprintf(feq, "%f X%d  %f X%d = 0\n", strength, id+2,  strength, id+4);  // a12 = -a21
        vector<double>eq1(id+5,0.0);
        eq1[id] = strength; eq1[id+4] = -strength;
        eqns.push_back(eq1);  rhs.push_back(0.0);
        vector<double>eq2(id+4,0.0);
        eq2[id+1] = strength; eq2[id+3] = strength;
        eqns.push_back(eq2);  rhs.push_back(0.0);
        double ss = 10.0;       // a strength to constraint the scale.
        fprintf(feq, "%f X%d = %f\n", ss, id+1, ss);  // a11 = 1.0
        vector<double>eq3(id+1,0.0);
        eq3[id] = ss;
        eqns.push_back(eq3);  rhs.push_back(ss);
        }
    }
fclose( feq );
//CI->Write( "Corr.xml" );
delete CI;
//system( "../svdfit eqns -print" );
//feq = FileOpenOrDie( "coeff", "r", flog );

vector<double>x;
vector<double>rslt;
int ns = LeastSquaresBySVD(eqns, x, rhs, rslt);
if( ns > 0 ) {
    printf("Singular values in least-squares fit!\n");
    fprintf(flog,"Singular values in least-squares fit!\n");
    //exit( 42 );
    }
for(int j=1; j<vp.size(); j++) {
    double scale = sqrt(pow(vp[j].tr.t[0], 2.0) + pow(vp[j].tr.t[1],2.0) );
    double angle = RadiansFromAffine( vp[j].tr ) * 180.0/PI;
    printf("Was %6.4f %7.3f: ", scale, angle); vp[j].tr.PrintTransform();
    if( TrOnly ) {
        vp[j].tr.t[2] = x[(j-1)*2];
        vp[j].tr.t[5] = x[(j-1)*2+1];
        }
    else {
        vp[j].tr.CopyIn( &x[(j-1)*6] );
	}
    scale = sqrt(pow(vp[j].tr.t[0], 2.0) + pow(vp[j].tr.t[1],2.0) );
    angle = RadiansFromAffine( vp[j].tr ) * 180.0/PI;
    printf(" is %7.4f %6.3f: ", scale, angle); vp[j].tr.PrintTransform();
    InvertTrans(vp[j].Inverse, vp[j].tr);
    }
//fclose(feq);
//
//Let's see how good the fit was

double errx = 0.0;
double erry = 0.0;
double err2 = 0.0;
double avg = 0.0;
double max_err = -1.0E30;
int npts = 0;
for(int i=0; i<eqns.size(); i++) {
    // does it contain an X term?  If so, this and the next equation form a pair
    bool xterm = false;
    for(int j=2; j<eqns[i].size(); j += 6)
	xterm |= (eqns[i][j] != 0.0);
    if( xterm ) {
        double xerr = rhs[i] - rslt[i];
        double yerr = rhs[i+1] - rslt[i+1];
        double emag = sqrt(xerr*xerr + yerr*yerr);
        errx += xerr; erry += yerr;
        avg += emag;
        err2 += emag*emag;
        max_err = fmax(max_err, emag);
	npts++;
        i++;
        }
    }
double rms = sqrt(err2/npts);
printf("%d points, x error %f, yerror %f, avg mag %f, RMS %f, max err %f\n",
 npts, errx/npts, erry/npts, avg/npts, rms, max_err);

if( rms < 0.01 ) {
    printf("Fit too good to be true! (%f)\n", rms);
    printf("Fit too good to be true! (%f)\n", rms);
    exit( 42 );
    }

fprintf(flog,"[%4.2f]", rms);
}

// routine to create b&w and color composite images.   Takes a target area, a list
// of tiles, and makes two images.  Returns number of valid pixels.
// 'col' specifies color for monochrome image, or 'W' for white
int CreateCompositeImage( int xmin, int ymin, int w2, int h2, vector<Picture> &vp,
 uint32 *bw, char col, uint32 *color, uint8 *map, vector<TForm> &tfs)
{
int npixels2 = w2*h2;
int nvpix = 0;
int sc = vp[0].scale;  // all scales must be the same, so this is OK
for(int j=0; j<npixels2; j++) {
    //bw[j] = 0;
    color[j] = 0xFF000000;
    int y = j/w2;                // coords within the image
    int x = j-w2*y;
    Point pt(x+xmin,y+ymin);     // coords in global space
    //printf("in global space %f %f, %d candidates\n", pt.x, pt.y, vp.size() );
    bool first = true; // for composite image, will take first match, since they
                       // are sorted by size.  For color image, combine all
    for(int k=0; k<vp.size(); k++) {
	Point p2(pt.x, pt.y);
	vp[k].Inverse.Transform( p2 );
        //printf("Image %d: %f %f\n", p2.x, p2.y);
	// is the point within the image?
	if( p2.x >= 0 && p2.x < vp[k].w-1 && p2.y >= 0 && p2.y < vp[k].h-1 ) {
	    int xll = int(p2.x);
	    int yll = int(p2.y);
            if( map != NULL ) { // the image is sub-divided by a map.
                // the map is in full scale coordinates
                int mx = xll*sc;
                int my = yll*sc;
		int id = map[mx+(vp[k].w*sc)*my];
                if( id >= 10 ) { // we have a better mapping
                    Point p3(pt); // global point
		    tfs[id-10].Transform( p3 );
		    if( p3.x >= 0 && p3.x < vp[k].w-1 && p3.y >= 0 && p3.y < vp[k].h-1 ) {
			p2 = p3;
			xll = int(p2.x);
			yll = int(p2.y);
			}
		    }
		}

	    double val = InterpolatePixel( p2.x, p2.y, vp[k].raster, vp[k].w );

        // create a negative image of higher contrast
	    int pix = int(127 + 2*(127-val) +0.5);  // rounding
		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;
	    //printf("set [%d,%d] to %f\n", x, y, val);
	    if( first ) {  // if the first of the image group to set this pixel
                if (col == 'W') // set R, G, and B
	            bw[x+w2*y] = 0xFF000000 | (pix<<16) | (pix << 8) | pix;
                else {
                    int shift = (col == 'R' ? 0 : (col == 'G' ? 8 : 16));
                    //printf("shift, pix %d %d\n", shift, pix);
                    //printf("was 0x%x\n", bw[x+w2*y]);
		    bw[x+w2*y] |= (pix << shift);
                    //printf(" is 0x%x\n", bw[x+w2*y]);
                    }
                first = false;
	        nvpix++;
                }
            int dcolor = k % 3; // drawing color
            // invert the pixels for the color image
            pix = 255 - pix;
	    color[x+w2*y] |= (pix << (8*dcolor));
            }
	}
    }
return nvpix;
}

// Find a bunch of points, spread out as far from center, that still fall within both pictures
// a = above, b = below, looking for points in a's coordinate system

int dxs[8] = {1, 1, 0, -1, -1, -1,  0,  1};
int dys[8] = {0, 1, 1,  1,  0, -1, -1, -1};
void FindMoreSpots(vector<Picture> &b, int j, vector<Picture> &a, int k, vector<Point> &centers, double xc, double yc)
{
double rad = RADIUS + 10; // 10 pixels for margin
for(int dir=0; dir<8; dir++) { // try 8 directions
    double x = xc;
    double y = yc;
    Point save(x,y);
    printf("------------------------------------------------------\n");
    for(int d=0; d<1000; d++) {
        x += dxs[dir];
        y += dys[dir];
        if( x > rad && x < a[k].w-rad && y > rad && y < a[k].h-rad ) { // inside a
	    Point p(x,y), s(x,y);
            printf(" In frame a, %f %f\n", p.x, p.y);
            a[k].tr.Transform( p ); // into global space
	    b[j].Inverse.Transform( p ); // back into b's space
            printf(" In frame b, %f %f\n", p.x, p.y);
            if (p.x > rad && p.x < b[j].w-rad && p.y > rad && y < b[j].h-rad)  // inside b
                save = s; // OK, it's good.  Save it
            }
	}
    centers.push_back(save);
    }
}


// read an existing local mapping file

int ReadMappingFile(vector<ffmap> &mapv, char *name, FILE *flog)
{
mapv.clear();  // remove any existing stuff from mapv
TiXmlDocument doc(name);
bool loadOK = doc.LoadFile();
printf("XML load gives %d\n", loadOK);
if (!loadOK)  {
    printf("Could not open XML file '%s'\n", name);
    fprintf(flog,"Could not open XML file '%s'\n", name);
    exit( 42 );
    }
TiXmlHandle hDoc(&doc);
TiXmlElement* child;
TiXmlHandle hRoot(0);

// block: should be <trakem2>
TiXmlNode *node=0;
node =doc.FirstChild();
child=hDoc.FirstChild("local_maps").FirstChild("entry").ToElement();
// should always have a valid root but handle gracefully if it does not
if( !child ) {
    printf(      "Either <local_maps> or <entry> is missing\n");
    fprintf(flog,"Either <local_maps> or <entry> is missing\n");
    exit( 42 );
    }
printf("child element value %s\n", child->Value() );
for( child; child; child=child->NextSiblingElement() ) {
    ffmap im;
    printf("Got a <entry>\n");
    im.fname = child->Attribute("overlap");
    im.mname = child->Attribute("map");
    //printf("z= %s\n", attr);
    // now look through all the <entry> elements in each
    TiXmlElement *c2;
    c2 = child->FirstChildElement("map");
    for( c2; c2; c2=c2->NextSiblingElement() ) {
        const char *tf  = c2->Attribute("transform");
        double a,b,c,d,e,f;
        sscanf(tf,"%lf %lf %lf %lf %lf %lf", &a, &b, &d, &e, &c, &f);
        printf("%7.4f %8.4f %12.2f\n%7.4f %8.4f %12.2f\n", a,b,c,d,e,f);
        //printf("%7.4f %8.4f %12.2f\n%7.4f %8.4f %12.2f\n", a,b,c,d,e,f);
        im.transforms.push_back( TForm(a, b, c, d, e, f) );
	}
    mapv.push_back(im);
    }
return 0;
}






// creates vectors of coordinates and values from within 'radius' of 'point'.

void CircleFromRaster(uint8 *raster, int w, int h, Point center, int radius, vector<Point> &pts, vector<double> &vals)
{
pts.clear();
vals.clear();
for(int i=0; i<w*h; i++) {
    int y = i / w;
    int x = i - w * y;
    Point pt(x,y);            // Coordinates in picture j
    Point d(pt.x-center.x, pt.y-center.y);
    double dist = sqrt(d.x*d.x + d.y*d.y);
    if( dist <= radius ) {
	//printf("pt: x,y=%f %f\n", pt.x, pt.y);
	pts.push_back(pt);
	vals.push_back(raster[i]);
	}
    //else  // so we can look at it...
	//above[k].raster[i] = 0;
    }
}

bool sort_by_x(const Point &a, const Point &b) { return a.x < b.x; }
bool sort_by_y(const Point &a, const Point &b) { return a.y < b.y; }

void MakeSubregions(vector<Point> &pts, double xmin, double xmax, double ymin, double ymax, vector<ConnRegion> &out)
{
// count the points, and find the BB of all within the region
//
vector<Point> vp;
int xll =  BIG, yll =  BIG;
int xur = -BIG, yur = -BIG;
for(int i=0; i<pts.size(); i++) {
    if( pts[i].x > xmin && pts[i].x <xmax && pts[i].y > ymin && pts[i].y < ymax ) {
	vp.push_back(pts[i]);
        xll = min(xll, int(pts[i].x));
        yll = min(yll, int(pts[i].y));
        xur = max(xur, int(pts[i].x));
        yur = max(yur, int(pts[i].y));
        }
    }
printf(" MakeSub: %d pts, [%d %d], [%d %d]\n", vp.size(), xll, xur, yll, yur);
if( vp.size() > 300000 ) { // too big; should be programmable
    // slice in X or Y?
    if( xur-xll > yur - yll ) { // wider than tall.  Slice at X
        sort(vp.begin(), vp.end(), sort_by_x);
	double split = vp[vp.size()/2].x + 0.5;  // make 2 equal pt-count groups
        printf("splitting in x at %f\n", split);
	MakeSubregions(vp, xll-0.5, split, yll-0.5, yur+0.5, out);
	MakeSubregions(vp, split, xmax+0.5,yll-0.5, yur+0.5, out);
        }
    else {
        sort(vp.begin(), vp.end(), sort_by_y);
	double split = vp[vp.size()/2].y + 0.5;
        printf("splitting in y at %f\n", split);
	MakeSubregions(vp, xll-0.5, xur+0.5, yll-0.5, split, out);
	MakeSubregions(vp, xll-0.5, xur+0.5, split, yur+0.5, out);
	}
    }
else { // OK size, just push the region onto the output stack
    ConnRegion cr;
    cr.pts = vp;  // pixels within the region
    cr.B.L = xll;
    cr.B.R = xur;
    cr.B.B = yll;
    cr.B.T = yur; // bounding box in original image
    out.push_back(cr);
    }
}

// we want edges between minl and 2*minl;
const int minl = 400;
const int minl2 = minl*minl;
bool EdgeIsShort(vector<vertex> &small, int i)
{
const int M = small.size();
int next = (i+1)%M;
return small[i].DistSqr( small[next] ) < minl2;
}



// try to draw a line from entry s to entry e of vertices.  See how big the error is...
double DistFrom(int s, int e, vector<vertex> &v)
{
const int N = v.size();
double dmax = 0.0;
for(int i=(s+1)%N; i != e; i = (i+1)%N)
    dmax = fmax(dmax, SegPointDist(v[s], v[e], v[i]) );
return dmax;
}





// classes for creating a queue
class qe {
  public:
    qe(int t, double c){to = t; cost = c;}
    bool operator<(const qe &a) const {return cost > a.cost;};  // priority is less if cost is higher
    int to;  // we know how to get to the node 'to' for total cost 'cost'
    double cost;
    };


class gr{
  public:
    gr(int xx, int yy, int b, double c){x = xx; y = yy; back = b; cost = c;};
    int x, y, back;
    double cost;
    };

double SegPointDist(gr &v0, gr &v1, gr &v2)
{return SegPointDist(v0.x, v0.y,  v1.x, v1.y,  v2.x, v2.y);}

// try to draw a line from entry s to entry e of vertices.  The cost of this
// approximation is the sum of all the errors at intermediate vertices.
double CostOfApprox(int s, int e, vector<gr> &v)
{
double sum = 0.0;
for(int i=s+1; i != e; i = i+1)
    sum += SegPointDist(v[s], v[e], v[i]);
return sum;
}
// given a vector of integer-valued points, create a bounding polygon that approximates it.
// We can afford a few pixel of errors (perhaps 5 pixels) and do not want too fine of a
// fragmentation of edges.
// Returns a vector of control points and triangles.
void CreateOutline(vector<Point> &pts, vector<vertex> &ControlPoints, vector<triangle> &tris)
{
// compute min and max.
IBox	B;
BBoxFromPoints( B, pts );
printf("region size is [%d %d] in x, [%d %d] in y\n", B.L, B.R, B.B, B.T);
// create a bit-map image
int w = B.R-B.L+3;  // a border of blanks all around
int h = B.T-B.B+3;
vector<uint8> map(w*h,0);  // create a map and fill it with zeros
for(int j=0; j<pts.size(); j++) {
    int x = int(pts[j].x) - B.L + 1;
    int y = int(pts[j].y) - B.B + 1;
    map[x + w*y] = 1;
    }
// start from the middle left middle, progress till we hit something
// We work in map coordinates since that is easier
int y0 = h/2;
int x0 = 0;
for(x0=0; x0<w && map[w*y0+x0] == 0; x0++)
    ;
if( x0 == w ) {
    printf("Can't find starting point??\n");
    exit( 42 );
    }
// Now work our way around CCW, till we get back the the beginning
vector<vertex> vertices;
int x = x0, y = y0, dir = 0;
for(; !(x == x0 && y == y0 && dir != 0); ) {
    //printf("-- %6d %6d %3d\n", x, y, dir);
    int op = (dir+4)%8;  // opposite to incoming direction
    for(int k=1; k<=8; k++) {  // one of these directions MUST be non-zero
        int j = (op+k)%8;
        int dx = dxs[j];
	int dy = dys[j];
	if( map[x+dx + w*(y+dy)] ) {
            vertex vtx(x,y,j);
            vertices.push_back(vtx);
	    x = x+dx;
	    y = y+dy;
	    dir = j;
	    break;
	    }
	}
    }

// Now compress all the identical edges.  Only save vertices whose direction is different than previous one
vector<vertex> small;
int N = vertices.size();
for(int i=0; i<N; i++) {
    int prev = (i+N-1)%N;  // since using mod N, adding N-1 is same as subtracting 1
    //printf("[%d].dir = %d, prev [%d].dir = %d\n", i, vertices[i].dir, prev, vertices[prev].dir);
    if( vertices[i].dir != vertices[prev].dir ) {
	vertex v = vertices[i];
        v.orig = i;       // was originally the ith vertex
        small.push_back(v);
	}
    }
printf("Had %d original vertices, now %d\n", N, small.size() );

// Now, if any edges are too long, split them
for(int i=0; i<small.size(); i++) {
    int next = (i+1)%small.size();
    int dx = small[next].x - small[i].x;
    int dy = small[next].y - small[i].y;
    double len = sqrt(dx*dx + dy*dy);
    //printf("len=%f\n", len);
    if( len > minl*2.0 ) {
        int ins = int(len/(minl*1.5)); // how many to insert
        for(int j=1; j<=ins; j++) {
	    double frac = double(j)/(ins+1);
            vertex nv(int(small[i].x+dx*frac), int(small[i].y+dy*frac), 0);
            printf("Inserting at (%d %d)\n", nv.x, nv.y);
            vector<vertex>::iterator it = small.begin();
	    small.insert(it+i+j, nv);
	    }
	}
    }

for(int i=0; i<small.size(); i++)
    printf("--After split edges: vertex (%d %d)\n", small[i].x, small[i].y);


// now find the lowest cost path from the first to last vertex.
//
// set all back pointers to -1, and all costs to infinity
vector<gr> graph;
for(int i=0; i<small.size(); i++)
    graph.push_back(gr(small[i].x, small[i].y, -1, double(BIG)));
// now duplicate the first entry after the last.  Since the first step will be a big one, would
// like the first point to be a junction between two long edges.  Could be assured by a rotation
// but will usually happen by our construction
graph.push_back(graph[0]);
// Push the [0] entry on the queue
priority_queue<qe> q;
q.push(qe(0, 0.0));
graph[0].cost = 0.0;

//
while (!q.empty()) {
    qe tp = q.top(); // highest priority == lowest cost
    q.pop();
    printf("can get from node 0 to %d for a total cost of %f\n", tp.to, tp.cost);
    if (tp.cost > graph[tp.to].cost) // we already knew how to get here for cheaper
	continue;
    if( tp.to == graph.size()-1 )
	break;
    // we need to examine all the fanouts of this node
    int s = tp.to;  // the source node
    for(int j=s+1; j<graph.size(); j++) {
        int dx = graph[j].x - graph[s].x;
	int dy = graph[j].y - graph[s].y;
	double len = sqrt(dx*dx + dy*dy);
	//printf("len=%f\n", len);
	if( minl <= len && len <= 2*minl ) {
	    double cost =  graph[s].cost + CostOfApprox(s, j, graph);
            if( cost < graph[j].cost ) {// we've found a better way to get there
		graph[j].cost = cost;
		graph[j].back = s;
		q.push(qe(j, cost));
	 	}
	    }
	}
    }

if( graph[graph.size()-1].cost == BIG ) {
    printf("No path to final vertex??\n");
    exit( 42 );
    }
//backtrack from the final vertex, adding output edges as we go
vector<lineseg> edges;
for(int i=graph.size()-1; i != 0; i = graph[i].back) {
    vertex v(graph[i].x, graph[i].y, 0);
    RemoveFromMap( map, w, h, v, minl );   // remove all points near existing points
    gr prev = graph[graph[i].back];
    double d=sqrt(pow(double(prev.x-graph[i].x),2.0) + pow(double(prev.y-graph[i].y),2.0) );
    printf("edge from (%d %d) to (%d %d),length %7.2f\n",prev.x, prev.y, graph[i].x, graph[i].y, d);
    edges.push_back(lineseg(prev.x, prev.y, graph[i].x, graph[i].y));
    }

// Now add any internal vertices
vector<vertex> internal;
for(int i=0; i<w*h; i++) {
     if(map[i]) {
        int y = i / w;
        int x = i - w * y;
        printf("Add vertex at (%d %d)\n", x, y);
	vertex newv(x,y,0);
        RemoveFromMap( map, w, h, newv, minl );
        internal.push_back(newv);
        }
    }

// now divide into triangles.  First, compile a list of (ordered) exterior edges.

// Now, take each edge one at a time.  Find the best triangle for that edge, then remove
// triangle from the figure.  Repeat until no edges left.  Result will be a list of
// control points, and triangles that refer to them.
ControlPoints.clear();
tris.clear();
for(; edges.size() > 0; ) {
    // print all edges in matlab format
    printf("\n--------- Start matlab format ------------------\n");
    for(int j=0; j<edges.size(); j++) {
	printf("x=[%d %d]; y=[%d %d]; plot(x,y); hold on;", edges[j].v[0].x, edges[j].v[1].x, edges[j].v[0].y, edges[j].v[1].y);
        if( j > 0 && (j%2) == 0 )
	    printf("\n");
	}
    printf("\n");
    // we'll use the last edge, since it's the easiest one to delete.
    vertex va(edges[edges.size()-1].v[0]);
    vertex vb(edges[edges.size()-1].v[1]);
    edges.pop_back();
    vertex vmid((va.x+vb.x)/2, (va.y+vb.y)/2, 0); // midpoint
    printf("working on (%d %d) -> (%d %d)\n", va.x, va.y, vb.x, vb.y);
    // find the closest vertex that's on the right side.
    int id = -2;   // will be -1 if it's an internal vertex, 0 for v[0] of an edge, 1 if v[1] of an edge
    int ii = -1;   // the index into either edges or internal
    double bestd = BIG;
    for(int j=0; j<edges.size(); j++) {
        vertex vc(edges[j].v[0]);
        double d = vmid.DistSqr( vc );
        //printf("#1 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if( LeftSide(va, vb, vc) && d < bestd && (!AnyCrossing(edges, va, vc)) && (!AnyCrossing(edges, vb, vc))  ) {
	    ii = j; id = 0; bestd = d;
	    }
        vc = edges[j].v[1];
        d = vmid.DistSqr( vc );
        //printf("#2 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if( LeftSide(va, vb, vc) && d < bestd && (!AnyCrossing(edges, va, vc)) && (!AnyCrossing(edges, vb, vc))  ) {
	    ii = j; id = 1; bestd = d;
	    }
        }
    for(int j=0; j<internal.size(); j++) {
        vertex vc(internal[j]);
        double d = vmid.DistSqr( vc );
        //printf("#3 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if( LeftSide(va, vb, vc) && d < bestd && (!AnyCrossing(edges, va, vc)) && (!AnyCrossing(edges, vb, vc))  ) {
	    ii = j; id = -1; bestd = d;
	    }
        }
    vertex vc = id < 0 ? internal[ii] : edges[ii].v[id];
    // now the triangle goes va->vb->vc->va
    printf("Triangle (%d %d) (%d %d) (%d %d)\n", va.x, va.y, vb.x, vb.y, vc.x, vc.y);
    triangle t;
    vector<vertex> temp; temp.push_back(va); temp.push_back(vb); temp.push_back(vc);
    for(int j=0; j<3; j++) {
	int k;
	for(k=0; k<ControlPoints.size(); k++)
	    if( ControlPoints[k] == temp[j] ) {
		t.v[j] = k;
		break;
	    }
        if( k == ControlPoints.size() ) { // did not find it; add it
            ControlPoints.push_back(temp[j]);
            t.v[j] = k;
	    }
        }
    tris.push_back(t);
    // for the edges vb->vc and vc->va, if they are in the list of edges, remove them
    // if they are not in the list, then add them.
    lineseg newone(vb,vc);
    int start = edges.size();
    for(int j=0; j<start; j++) {
	if( edges[j] == newone ) {
	    edges.erase(edges.begin()+j);
	    break;
	    }
        }
    if(edges.size() == start)
	edges.push_back(lineseg(vc,vb));
    lineseg newtwo(vc,va);
    start = edges.size();
    for(int j=0; j<start; j++) {
	if( edges[j] == newtwo ) {
	    edges.erase(edges.begin()+j);
	    break;
	    }
        }
    if(edges.size() == start)
	edges.push_back(lineseg(va,vc));
    // if the edge came from the internal list, it's internal no longer, so remove it.
    if( id == -1 )
	internal.erase(internal.begin() + ii);
    }
printf("Got %d triangles using %d control points\n", tris.size(), ControlPoints.size() );

// Convert back to input coordinates
for(int j=0; j<ControlPoints.size(); j++) {
    ControlPoints[j].x += (B.L - 1 );
    ControlPoints[j].y += (B.B - 1 );
    }
}



// Tries to find a mapping from one connected region to another.
// a = above, b = below, we look for a mapping from a->b
// The picture has the whole image, the connected region tells what part we are concerned with
void RegionToRegionMap(Picture &pa, ConnRegion &a, Picture &pb, ConnRegion &b, uint16 *ids, ffmap &map1,
FILE *flog)
{

// now we need to find how each of these regions map to the picture below
// if indeed they do at all.
int w = pa.w;		// must be the same for a and b, verified earlier
int h = pa.h;
int sc = pa.scale;
int fullw = pa.w*sc;
int fullh = pa.h*sc;
int npixels = b.pts.size();
printf(" In RegionToRegion - scale = %d\n", sc);
MeanStd mm;
for(int i=0; i<npixels; i++) {
    int ix = int(b.pts[i].x);
    int iy = int(b.pts[i].y);
    mm.Element(pb.raster[ix + w*iy]);
    }
printf("below: %d real pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/(w*h));
double mean2, std2;
mm.Stats(mean2, std2);  // find existing statistics
printf("On the target image,  mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy, with image in lower left 2Kx2K, so FFT will work.
// Background pixels map to zero.
// Also, create a map that just tells whether each point is in the region, or not
vector<double> image2(4096*4096, 0.0);
vector<char> IsIn(w*h,0);
for(int i=0; i<npixels; i++) {
    int x = int(b.pts[i].x);
    int y = int(b.pts[i].y);
    double pix = (pb.raster[x + w*y] - mean2)/std2;
    if (abs(pix) > 2)    // outlying pixel value
	pix = 0.0;  // will be set to the mean value (0)
    image2[x + 4096*y] = pix;
    IsIn[x + w*y] = 1;
    }

// if needed, make a bigger copy that can be used for full resolution
//vector<double>image3;
//if( sc > 1 ) { // then we have a copy with more than 2Kx2K resolution
//    int np3 = npixels*sc*sc;
//    MeanStd m3;
//    for(int i=0; i<np3; i++) {
//	if( below[0].orig[i] > 0 )
//	    m3.Element(below[0].orig[i]);
//	}
//    double mean3, std3;
//    m3.Stats(mean3, std3);  // find existing statistics
//    printf("Of the full resolution target image,  mean= %f and std dev = %f\n", mean3, std3);
//    image3.push_back(0.0);
//    image3.insert(image3.begin(), 4096*4096-1, 0.0); // is there a better way?
//    int bigw = below[0].w*sc;
//    for(int i=0; i<np3; i++) {
//	int y = i / bigw;
//	int x = i - bigw * y;
//	double pix = below[0].orig[i];
//	if (pix == 0)    // background pixels
//	    pix = mean3;  // will be set to the mean value (0)
//	image3[x + 4096*y] = (pix-mean3)/std3;
//	}
//    }

// create a map image.  For every pixel in A, tell the closest control point.  0-9  reserved for 'no data',
// for cases such as folds, off the target image, piece too small, etc. .
// Otherwise if the pixel has value N, the closest point is N-10
//

// split the connected region into sub-regions of limited size.
printf("\n-----Starting region ----\n");
vector<vertex> ControlPoints;
vector<triangle> tris;
CreateOutline(a.pts, ControlPoints, tris);

vector<ConnRegion> subs(tris.size());

// assign each point to one of the sub-regions
for(int k=0; k<a.pts.size(); k++) {
    Point pt = a.pts[k];
    int reg = BestTriangle(tris, ControlPoints, pt);
    subs[reg].pts.push_back(pt);
    subs[reg].B.L = min(subs[reg].B.L, int(pt.x));
    subs[reg].B.B = min(subs[reg].B.B, int(pt.y));
    subs[reg].B.R = max(subs[reg].B.R, int(pt.x));
    subs[reg].B.T = max(subs[reg].B.T, int(pt.y));
    }
// Now, try to match each constituent triangle by itself.  Save the results in vector tcorr and dx,dy of vector 'subs'
vector<double> tcor;
for(int i=0; i<subs.size(); i++) {
    printf("\n %d pts, x=[%d %d], y = [%d %d]\n", subs[i].pts.size(), subs[i].B.L, subs[i].B.R, subs[i].B.B, subs[i].B.T);
    vector<double> vals;
    for(int k=0; k<subs[i].pts.size(); k++) {
	int ix = int(subs[i].pts[k].x);
	int iy = int(subs[i].pts[k].y);
	vals.push_back(pa.raster[ix + w*iy]);
	}
    Normalize(vals);
    //WriteDebugImage(subs[i].pts, vals, 100 + 100*j + i);
    double	c = CorrPatchToImage( subs[i].dx, subs[i].dy, subs[i].pts, vals, image2, 0, 0, 4000, true );
    tcor.push_back(c);
    printf(" c = %f, deltas are %f %f\n", c, subs[i].dx, subs[i].dy);
    }
// Now, find a 'concensus' transform - one that many triangles agree on
int best = -1;
int cbest = 0;  // best one found
double corr_of_best = 0.0;   // and correlation of that one
for(int i=0; i<subs.size(); i++) {
    int close = 0; // count how many are close to this offset (use 15 pixels)
    for(int k=0; k<subs.size(); k++) {
	if( fabs(subs[i].dx-subs[k].dx) < 15.0 && fabs(subs[i].dy - subs[k].dy) < 15.0 )
	    close++;
	}
    if( close > cbest || (close == cbest && tcor[i] > corr_of_best) ) {  // more close matches than any other point?
	best = i;
	cbest = close;
	corr_of_best = tcor[i];
	}
    }
printf("OK, best correlation is sub-region %d, %d are close, correlation %f\n", best, cbest, tcor[best]);
TForm tr_guess; // our best guess at the transform relating
if( tcor[best] > 0.25 ) {  // just guessing here
    vector<double> vals;                             // should make a routine
    for(int k=0; k<subs[best].pts.size(); k++) {
	int ix = int(subs[best].pts[k].x);
	int iy = int(subs[best].pts[k].y);
	vals.push_back(pa.raster[ix + w*iy]);
	}
    Normalize(vals);
    double c=ImproveCorrelation(subs[best].pts, vals, image2, subs[best].dx, subs[best].dy, tr_guess, flog);
    printf("Best guess transform: "); tr_guess.PrintTransform();

    // now using the guess, keep only those points that map to the below image.  (this could be off
    // by 10 pixels or so, so we'll need to fix it later, but it's a good guess
    vector<Point> newpts;
    for(int i=0; i<a.pts.size(); i++) {
	Point p(a.pts[i]);
	tr_guess.Transform( p );
        int ix = int(p.x);   // approximate target pixel
        int iy = int(p.y);
	if( ix >= 0 && ix < w && iy >= 0 && iy < h && IsIn[ix + w*iy] != 0 )
	    newpts.push_back(a.pts[i]);
	}
    // some quality checks should go here
    printf("Roughly %d pixels map to the image below\n", newpts.size());
    // Now divide again, could be very different
    CreateOutline(newpts, ControlPoints, tris);
    // for each triangle, create a matrix that transforms to barycentric coordinates
    for(int k=0; k<tris.size(); k++) {
	vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
	printf("Tri: (%d %d) (%d %d) (%d %d)\n", v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
	double a[3][3];
	a[0][0] = v0.x; a[0][1] = v1.x; a[0][2] = v2.x;
	a[1][0] = v0.y; a[1][1] = v1.y; a[1][2] = v2.y;
	a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
	Invert3x3Matrix(tris[k].a, a);
	}
    for(int k=0; k<tris.size(); k++) {  // write in matlab form
	vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
	printf("x=[%d %d %d %d]; y=[%d %d %d %d]; plot(x,y); hold on;\n", v0.x, v1.x, v2.x, v0.x,
	 v0.y, v1.y, v2.y, v0.y);
	}
    subs.clear(); subs.insert(subs.begin(), tris.size(), ConnRegion());

    // assign each point to one of the sub-regions
    int next_id = map1.transforms.size() + 10;
    vector<vector<double> >lambda;
    for(int k=0; k<newpts.size(); k++) {
	Point pt = newpts[k];
	int reg = BestTriangle( tris, ControlPoints, pt );
	subs[reg].pts.push_back(pt);
	subs[reg].B.L = min(subs[reg].B.L, int(pt.x));
	subs[reg].B.B = min(subs[reg].B.B, int(pt.y));
	subs[reg].B.R = max(subs[reg].B.R, int(pt.x));
	subs[reg].B.T = max(subs[reg].B.T, int(pt.y));
	// Map the point into barycentric coordinates for its triangle
	//vertex v0 = ControlPoints[tris[reg].v[0]], v1 = ControlPoints[tris[reg].v[1]], v2 = ControlPoints[tris[reg].v[2]];
	//printf(" pt (%f %f) in tri: (%d %d) (%d %d) (%d %d)\n", pt.x, pt.y, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
	double l[3];
	l[0] = tris[reg].a[0][0]*pt.x + tris[reg].a[0][1]*pt.y + tris[reg].a[0][2]*1.0;
	l[1] = tris[reg].a[1][0]*pt.x + tris[reg].a[1][1]*pt.y + tris[reg].a[1][2]*1.0;
	l[2] = tris[reg].a[2][0]*pt.x + tris[reg].a[2][1]*pt.y + tris[reg].a[2][2]*1.0;
	//printf("prop %f %f %f\n", l[0], l[1], l[2]);
	vector<double>lam(ControlPoints.size(), 0.0);
	for(int m=0; m<3; m++)
	    lam[tris[reg].v[m]] = l[m];
	//for(int m=0; m<lam.size(); m++)
	    //printf("%.4f ", lam[m]);
	//printf("\n");
	lambda.push_back(lam);
        for(int x=0; x<sc; x++)  // array ids is at original scale
	    for(int y=0; y<sc; y++)
	        ids[sc*int(pt.x)+x + fullw*(sc*int(pt.y)+y)] = next_id + reg;
	}
    // optimize the mesh.
    vector<double> spv;  // create a normalized array of source pixel values
    for(int k=0; k<newpts.size(); k++) {
	int ix = int(newpts[k].x);
	int iy = int(newpts[k].y);
	spv.push_back(pa.raster[ix + w*iy]);
	}
    Normalize(spv);
    // next, transform the control points to the target frame
    vector<Point> cpts; // make a real-valued copy
    for(int k=0; k<ControlPoints.size(); k++)
	cpts.push_back(Point(ControlPoints[k].x, ControlPoints[k].y));
    vector<Point> orig = cpts;
    tr_guess.Transform( cpts );
    // call the optimizer to find a better spot for the control points
    ImproveControlPts(cpts, lambda, spv, image2, 4096, 4096, flog);
    // compute a transform, and a center, for each triangle
    for(int k=0; k<tris.size(); k++) {
	int i0 = tris[k].v[0];
	int i1 = tris[k].v[1];
	int i2 = tris[k].v[2];
	// Now, find a transformation that maps ORIG into the new cpts
	// first, create a transform that maps a unit right triangle to the original pts
	TForm o(orig[i1].x-orig[i0].x, orig[i2].x-orig[i0].x, orig[i0].x,
	 orig[i1].y-orig[i0].y, orig[i2].y-orig[i0].y, orig[i0].y);
	// now one that maps the final optimized control points to the unit right triangle
	TForm c(cpts[i1].x-cpts[i0].x, cpts[i2].x-cpts[i0].x, cpts[i0].x,
	 cpts[i1].y-cpts[i0].y, cpts[i2].y-cpts[i0].y, cpts[i0].y);
	// now, to get from the original to the final, apply o^-1, then c;
	TForm oi,t;
	InvertTrans(oi, o);
	//TForm temp;
	MultiplyTrans(t, c, oi);
	t.PrintTransform();
	map1.transforms.push_back(t);
	double sumx = orig[i0].x + orig[i1].x + orig[i2].x;
	double sumy = orig[i0].y + orig[i1].y + orig[i2].y;
	map1.centers.push_back(Point(sumx/3.0, sumy/3.0));
	}
    }

for(int i=0; i<map1.centers.size(); i++) {
    printf(" center %f %f :", map1.centers[i].x, map1.centers[i].y);
    map1.transforms[i].PrintTransform();
    }
}

// The main routine the pipeline should call.
//
// Ntrans				- returned tform count
// array_of_transforms	- array of these values
// map_mask				- <10=no map, else tform[pix-10]
// AboveRaster			- the higher layer
// fold_mask_above		- 0=fold, 1,2,3...=region #
// BelowRaster			- the lower layer
// fold_mask_below		- 0=fold, 1,2,3...=region #
// w, h					- size of all images
// fout					- most output goes here
// flog					- some details written here
//
void PipelineDeformableMap(
	int				&Ntrans,
	double*			&array_of_transforms,
	uint16*			map_mask,
	const uint8*	AboveRaster,
	const uint8*	fold_mask_above,
	const uint8*	BelowRaster,
	const uint8*	fold_mask_below,
	int				w,
	int				h,
	FILE*			fout,
	FILE*			flog )
{
Picture above, below; // convert to this for sake of compatibility
	above.SetExternal( AboveRaster, w, h );
	below.SetExternal( BelowRaster, w, h );

//create smaller (at most 2Kx2K pixels, if original was bigger
above.DownsampleIfNeeded( flog );
below.DownsampleIfNeeded( flog );

// Create the connected region lists.  Note that the connected region lists are always at the reduced resolution, if this is used.
	vector<ConnRegion>	Acr, Bcr;

	ConnRgnsFromFoldMask( Acr, fold_mask_above,
		w, h, above.scale, BIG, flog );

	ConnRgnsFromFoldMask( Bcr, fold_mask_below,
		w, h, below.scale, BIG, flog );

memset( map_mask, 0, w * h * sizeof(uint16) );

ffmap map1;  // data structure holding feature to feature mapping
// now find all the mappings from each connected region to each connected region
for(int i = 0; i <Acr.size(); i++) {
    for(int j=0; j < Bcr.size(); j++) {
        printf("Looking for mapping from region %d of 'above' to region '%d' of below\n", i, j);
        RegionToRegionMap(above, Acr[i], below, Bcr[j], map_mask, map1, flog);
	}
    }
// Now package the transforms in the way the MatLab code would like to see them.
printf("Modifying transforms: above.scale=%d\n", above.scale);
Ntrans = map1.transforms.size();
array_of_transforms = (double *) malloc(Ntrans*6*sizeof(double));
for(int i=0; i<Ntrans; i++) {
    map1.transforms[i].MulXY( above.scale );
    map1.transforms[i].ToMatlab(); // re-order the components
    map1.transforms[i].CopyOut( array_of_transforms + i*6 );
    }
}


// Here we convert an input image to a map where,
// 0 = on fold, 1=region 1, 2=region 2, etc.
// Normally the pipeline does this, but we do it here
// for standalone testing.
//
// Some issues here:
//	- cr is not returned, so is inefficient.
//	- Assumes fold has value exactly zero.
//	- Unclear if -thresh is high enough to exclude fold value
//		-mean/std for arbitrary stats.
//	- Fold may be counted as a region cr??
//	- foldMask coloring (0 = fold) only works if folds not in cr.
//	- Experiment section appears to have no purpose.
//
void ImageToFoldMap( Picture &pic, uint8* foldMask )
{
	vector<ConnRegion>	cr;
	vector<double>		v;
	uint8				*raster	= pic.raster;
	double				mean, std;
	double				thresh	= 4.0;
	int					w		= pic.w;
	int					h		= pic.h;
	int					npixels	= w * h;
	int					nbig	= 0;
	int					i;

// Experiment: overlap with multiple regions

	bool within_section = false;
	if( within_section ) {
		SetWithinSectionBorders( foldMask, w, h );
		return;
	}

// Find mean and std dev of non-zero pixels

	StatsRasterNonZeros( raster, npixels, mean, std );

// Make a normalized copy of the raster.
// Do not copy any pixels next to black (val = 0) pixels,
// since they themselves are dubious.

	v.reserve( npixels );

	for( i = 0; i < npixels; ++i ) {

		int	y	= i / w;
		int	x	= i - w * y;
		int	pix;

		if( !(pix	= raster[i]) )
			;
		else if( x-1 >= 0 && raster[i-1] == 0 ) pix = 0;
		else if( x+1 <  w && raster[i+1] == 0 ) pix = 0;
		else if( y-1 >= 0 && raster[i-w] == 0 ) pix = 0;
		else if( y+1 <  h && raster[i+w] == 0 ) pix = 0;

		v.push_back( (pix-mean)/std );
	}

// Set threshold for connected pixels

	if( mean - thresh * std < 1 ) {

		thresh = mean / std * 0.95;

		printf(
		"Foldmap: Forced to reduce threshold to %f.\n", thresh );
	}

	const char *p = getenv( "FoldMaskThreshold" );

	if( p ) {

		thresh = atof( p );

		printf(
		"Foldmap: Environment variable over-ride of"
		" threshold to %f.\n", thresh );
	}

// Now find the connected regions

	for( i = 0; i < npixels; ++i ) {

		if( v[i] > -thresh ) {

			ConnRegion	c;
			int			npts;

			npts = Propagate( c.pts, v, w, h, i,
					-thresh, -thresh - 1.0 );

			printf(
			"ImageToFoldMap: ConnRegion with %d pixels.\n", npts );

			if( npts > 100000 ) {
				cr.push_back( c );
				++nbig;
			}
		}
	}

// Final accounting

	SetBoundsAndColors( cr, foldMask, w, h );
}


// reads in two files and an approximate transform.  Creates a map of more exact
// transforms, and updates (or writes) the XML file with these transforms
// arguments:
//    Arg[1] = below.tif
//    arg[2] = above.tif
//    arg[3] = transform from "above" to "below".  Should be within a few pixels
//    arg[4] = name of xml file
int main(int argc, char* argv[])
{
FILE	*flog = FileOpenOrDie( "two.log", "a" );

if( argc < 5 ) {
    printf("Usage: two <below-image> <above-image> <transform above->below> <xml file>\n");
    exit( 42 );
    }
bool WriteBack = false;   // overwrite the original?

time_t t0 = time(NULL);
char atime[32];
strcpy(atime, ctime(&t0));
atime[24] = '\0';  // remove the newline
fprintf(flog,"two: %s ", atime);
for(int i=1; i<argc; i++) {
    // process arguments here
    fprintf(flog,"%s ", argv[i]);
    if( strcmp(argv[i],"-o") == 0 && i+1 < argc )
	OVERLAP = atoi(argv[i+1]);
    if( strcmp(argv[i],"-t") == 0 && i+1 < argc )
	THRESHOLD = atof(argv[i+1]);
    if( strcmp(argv[i],"-r") == 0 && i+1 < argc )
	RADIUS = atoi(argv[i+1]);
    if( strcmp(argv[i],"-f") == 0 )
	WriteBack = true;
    if( strcmp(argv[i],"-tr") == 0 )
	TrOnly = true;
    }
fflush(flog);
printf("Overlap %d, radius %d, threshold %f\n", OVERLAP, RADIUS, THRESHOLD);

// for historical reasons, we think of an array of 'above' images and an array of 'below'
// images.  Here we will have only one of each.  Set the 'below' tranform to 1:1
// and the 'above' transform to that provided.
//
vector<Picture> below;
vector<Picture> above;

Picture p;
p.dist_from_center = 0.0;  // prevent reference to undefined fields

p.fname = argv[1];
p.tr = TForm();  // make explicit 1<->1 transform; already true from initialization
below.push_back(p);

p.fname = argv[2];
double a,b,c,d,e,f;
sscanf(argv[3],"%lf %lf %lf %lf %lf %lf", &a, &b, &d, &e, &c, &f);
p.tr = TForm(a, b, c, d, e, f);
p.tr.PrintTransform();
above.push_back(p);

// sort the higher layer vector by dist from center.  Read each of them in
// Here we do not use arrays, so this is over-elaborate, but it should work.
sort(above.begin(), above.end());
for(int j=0; j<above.size(); j++) {
    //printf("Above: z=%d, distance=%f\n", above[j].z, above[j].dist_from_center);
	above[j].LoadOriginal( argv[2], flog, false );
    InvertTrans(above[j].Inverse, above[j].tr);
    }
// same for the layer below
sort(below.begin(), below.end());
for(int j=0; j<below.size(); j++) {
    //printf("Below: z=%d, distance=%f\n", below[j].z, below[j].dist_from_center);
    below[j].LoadOriginal( argv[1], flog, false );
    InvertTrans(below[j].Inverse, below[j].tr);
    }
if( above[0].w != below[0].w || above[0].h != below[0].h ) {
    printf("Different scales for input pictures\n");
    fprintf(flog,"Different scales for input pictures\n");
    exit( 42 );
    }

int w = above[0].w;
int h = above[0].h;
int npixels = w * h;
uint8 *fold_mask_above = (uint8*)malloc(npixels*sizeof(uint8));
uint8 *fold_mask_below = (uint8*)malloc(npixels*sizeof(uint8));
uint16 *rmap = (uint16*)malloc(npixels*sizeof(uint16));

ImageToFoldMap(below[0], fold_mask_below);  //pipeline does this in normal operation
Raster8ToTif8( "bmap.tif", fold_mask_below, w, h );
ImageToFoldMap(above[0], fold_mask_above);  //but we will do it in standalone mode
Raster8ToTif8( "amap.tif", fold_mask_above, w, h );

//Now call the routine just as
int		Ntrans = 0;
double	*array_of_transforms;

PipelineDeformableMap(
	Ntrans, array_of_transforms, rmap,
	above[0].raster, fold_mask_above,
	below[0].raster, fold_mask_below,
	w, h, flog, flog );

uint8 *fake = (uint8*)malloc(npixels*sizeof(uint8)); // convert to 8 bit for debugging, since gimp can't do 16.
for(int i=0; i<npixels; i++)
    fake[i] = rmap[i];
Raster8ToTif8( "map.tif", fake, w, h );
printf("Got %d mapping regions\n", Ntrans);

// convert the double array back to an array of tforms.
TForm* tfs = new TForm[Ntrans];
TForm* ifs = new TForm[Ntrans];  // inverse of these transforms

for(int i=0; i<Ntrans; i++) {
    printf("Transform %3d:", i);
    tfs[i].CopyIn( array_of_transforms + i*6 );
    tfs[i].FromMatlab();
    tfs[i].PrintTransform();
    InvertTrans(ifs[i], tfs[i]);
    printf("\n");
    }






/* --------- */
/* ABOverlay */
/* --------- */

// Here we write out a diagnostic image.  This requires several steps:
int xmin = -1000, xmax = w+1000, ymin = -1000, ymax = h+1000;
uint32 w2 = uint32(xmax-xmin+1), h2=uint32(ymax-ymin+1);   // These describe the comparison image
uint32 *raster2;
size_t npixels2 = w2 * h2;
printf("  w %d, h %d, %d pixels in comparison image; w,h,npixels %d %d %d\n", w2, h2, npixels2, w, h, npixels);

raster2 = (uint32*)RasterAlloc( npixels2 * sizeof(uint32 ));

float *rt  = (float*)RasterAlloc( npixels2 * sizeof(float) );
									 // since we will write partial pixel
									 // values here, and so should not truncate
// Find a center for each region.   We have a map for image 'above'.  Make one for image 'below'
vector<uint16>bmap(npixels,0);
printf("Going to allocate the counts\n");
vector<int>  pcount(Ntrans, 0);
printf("allocated the counts\n");
vector<Point> sums(Ntrans, Point(0.0,0.0));
printf("allocated the sums\n");
fflush(stdout);
for(int x=0; x<w; x++) {
    for(int y=0; y<h; y++) {
        int n = rmap[x + w*y];
        if( n-10 >= Ntrans ) {
	    printf("odd - n=%d, Ntrans=%d\n", n, Ntrans);
            exit(-1);
            }
	if( n >= 10 ) {
	    sums[n-10].x += x;
            sums[n-10].y += y;
            pcount[n-10]++;
            Point p(x,y);
	    tfs[n-10].Transform( p );
            int ix = int(p.x);
            int iy = int(p.y);
            if( ix >= 0 && ix < w && iy >= 0 && iy < h )
		bmap[ix + w*iy] = n;
            }
        }
    }
for(int i=0; i<Ntrans; i++) {
    sums[i].x = sums[i].x/pcount[i];
    sums[i].y = sums[i].y/pcount[i];
    printf("region %d, center (%f %f) in 'above', npix=%d\n", i, sums[i].x, sums[i].y, pcount[i]);
    tfs[i].Transform( sums[i] );
    printf("   maps to (%f %f) in image 'below'.\n", sums[i].x, sums[i].y);
    }

for(int k=0; k<npixels2; k++) {
    raster2[k] = 0xFF000000;       // set the alpha channel
    rt[k] = 0.0;
    }

// for the image below.  For each pixel, look at the map. If the map is defined, write the data into that spot
// of the image.  If it is not defined, then find the closest center and use that transform.  If using that transform
// you fall within the below image, then do nothing (we wish to see through).  Otherwise copy it (we wish to see the
// surroundings).
//
for(int i=0; i<npixels; i++) {
    int iy = i / w;
    int ix = i - w * iy;
    int mv = bmap[i];  // get the map value (map is at full size
    Point p(ix, iy);
    if (mv >= 10)            // transform is defined, see where it goes
	ifs[mv-10].Transform( p );
    else { // not strictly defined...
	double dbest = 1.0E30;
        int best = -1;
        for(int q=0; q < Ntrans; q++) {
	    double dx = ix-sums[q].x;
            double dy = iy-sums[q].y;
	    double d = dx*dx + dy*dy;   // no need for sqrt, it's monotonic
	    if( d < dbest ) {
		dbest = d;
	        best = q;
	        }
	    }
        // now see, if using the best transform, if it's in the target image
        ifs[best].Transform( p );
        if( p.x >= 0.0 && p.x < w && p.y >= 0.0 && p.y < h )
	    continue; // do not want to set it.
	}
    // so now p is point, transformed into image 'above'.  Subtract (xmin, ymin) to get coordinate in image
    p.x -= xmin;
    p.y -= ymin;
    DistributePixel( p.x, p.y, below[0].raster[ix+w*iy], rt, w2, h2 );
}
// transfer temp to the RED channel.  Invert and stretch contrast
for(int i=0; i<npixels2; i++) {
    if( rt[i] > 0.0 ) {
        int pix = int(127 + (127-rt[i])*2);
		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;
        raster2[i] |= pix;
        }
    }
// now just copy the image 'above' to the GREEN channel.
for(int i=0; i<npixels; i++) {
    int iy = i / w;
    int ix = i - w * iy;
    int pix = 127 + (127 - above[0].raster[i]) * 2;  // invert and contrast stretch
	if( pix < 0 )
		pix = 0;
	else if( pix > 255 )
		pix = 255;
    ix -= xmin;
    iy -= ymin;
    if( ix >= 0 && ix < w2 && iy >= 0 && iy < h2 )  // This test should not be needed
	raster2[ix + w2*iy] |= (pix << 8);
    }

// Write it out
Raster32ToTifRGBA( "comp.tif", raster2, w2, h2 );

exit(-1);  // temporary






exit(0);

// allocate two fold arrays, one output array

// start the processing of images here
for(int i=0; i<above.size(); i++)
    above[i].DownsampleIfNeeded( flog );
for(int i=0; i<below.size(); i++)
    below[i].DownsampleIfNeeded( flog );

// Here, define the connected region(s) on the above layer.
// First, find the mean and standard deviation of the 'real' (non-zero) pixels
w = above[0].w;
h = above[0].h;
uint8 *raster = above[0].raster;
npixels = w * h;

// Find mean and std dev of non-zero pixels
	double	mean, std;
	StatsRasterNonZeros( raster, npixels, mean, std );

// Make a normalized copy of the raster.  Do not copy any pixels next to black pixels
// since they themselves are dubious.
// Make 2 copies since the first (v) will be destroyed during connected region evaluation.
vector<double> v;
for(int i=0; i<npixels; i++) {
    int y = i / w;
    int x = i - w * y;   // pixels next to background pixels should also be black
    int pix = raster[i];
    if (x-1 >= 0 && raster[i-1] == 0 ) pix = 0;
    if (x+1 <  w && raster[i+1] == 0 ) pix = 0;
    if (y-1 >= 0 && raster[i-w] == 0 ) pix = 0;
    if (y+1 <  h && raster[i+w] == 0 ) pix = 0;
    v.push_back((pix-mean)/std);
    }

// Now find the connected regions
int nbig = 0;
vector<ConnRegion> cr;
double thresh = 4.0;
if( mean - thresh*std < 1 ) {
    thresh = mean/std*0.95;
    printf("Forced to reduce threshold to %f\n", thresh);
    }
int start = 0;
// find all connected regions of 'reasonable' values.  The background should be about
// -4 or -5 on this scale.
for(int k=0;;){
    int i;
    for(i=start; i<npixels; i++)  // find first good pixel
	 if( v[i] > -thresh )
	    break;
    if (i >= npixels)             // if there are not any, then all done
	break;
    // found at least one pixel.  Find all connected ones.
    //printf("found v[%d] = %f\n", i, v[i]);
    start = i+1;  // next time, start here...
    ConnRegion c;
    stack<int> st;
    //printf("push %d\n", i);
    st.push(i);
    while(!st.empty()) {
	int j = st.top();
	//printf("pop %d, val %f\n",j, v[j]);
	st.pop();
	//if( fabs(v[j]) < thresh ) {
	if( v[j] > -thresh ) {
	    int y = j / w;
	    int x = j - w * y;
	    //printf("j, x, y = %d %d %d, v=%f\n", j, x, y, v[j]);
	    Point p(x,y);
	    c.pts.push_back(p);
	    v[j] = -5.0;
	    if (x-1 >= 0) st.push(j-1);
	    if (x+1 < w)  st.push(j+1);
	    if (y-1 >= 0) st.push(j-w);
	    if (y+1 < h)  st.push(j+w);
	    }
	}
    printf("Connected region of %d pixels\n", c.pts.size());
    if( c.pts.size() > 200000 ) {
	nbig++;
	cr.push_back(c);
	BBoxFromPoints( cr[k].B, cr[k].pts );
	printf("region size is [%d %d] in x, [%d %d] in y\n",
	 cr[k].B.L, cr[k].B.R, cr[k].B.B, cr[k].B.T);
	k++;
        }
    }

// Now find the biggest chunk.
int bc;  // the biggest chuck
int bc_size = -1; // and its size
for(int k=0; k<cr.size(); k++) {
    //printf(" looking %d, %d\n", cr[k].pts.size(), bc_size );
    if( int(cr[k].pts.size()) > bc_size ) {
	bc_size = cr[k].pts.size();
	bc = k;
	}
    }
if( bc_size <= 0 ) {
    printf("This cannot happen - bc_size=%d\n", bc_size);
    fprintf(flog,"This cannot happen - bc_size=%d\n", bc_size);
    exit( 42 );
    }
printf("connected region %d is the biggest\n", bc);

// now we need to find how each of these regions map to the picture below
// if indeed they do at all.
int sc = above[0].scale;
int fullw = above[0].w*sc;
int fullh = above[0].h*sc;
npixels = below[0].w*below[0].h;
MeanStd mm;
for(int i=0; i<npixels; i++) {
    if( below[0].raster[i] > 0 )
	mm.Element(below[0].raster[i]);
    }
printf("below: %d real pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/npixels);
double mean2, std2;
mm.Stats(mean2, std2);  // find existing statistics
printf("Of the target image,  mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy, with image in lower left 2Kx2K, so FFT will work.
// Background pixels map to zero
vector<double> image2(4096*4096, 0.0);
for(int i=0; i<npixels; i++) {
    int y = i / below[0].w;
    int x = i - below[0].w * y;
    double pix = (below[0].raster[i] - mean2)/std2;
    if (abs(pix) > 2)    // outlying pixel value
	pix = 0.0;  // will be set to the mean value (0)
    image2[x + 4096*y] = pix;
    }

// if needed, make a bigger copy that can be used for full resolution
vector<double>image3;
if( sc > 1 ) { // then we have a copy with more than 2Kx2K resolution
    int np3 = npixels*sc*sc;
    MeanStd m3;
    for(int i=0; i<np3; i++) {
	if( below[0].original[i] > 0 )
	    m3.Element(below[0].original[i]);
	}
    double mean3, std3;
    m3.Stats(mean3, std3);  // find existing statistics
    printf("Of the full resolution target image,  mean= %f and std dev = %f\n", mean3, std3);
    image3.push_back(0.0);
    image3.insert(image3.begin(), 4096*4096-1, 0.0); // is there a better way?
    int bigw = below[0].w*sc;
    for(int i=0; i<np3; i++) {
	int y = i / bigw;
	int x = i - bigw * y;
	double pix = below[0].original[i];
	if (pix == 0)    // background pixels
	    pix = mean3;  // will be set to the mean value (0)
	image3[x + 4096*y] = (pix-mean3)/std3;
	}
    }

// create a map image.  For every pixel in A, tell the closest control point.  0-9  reserved for 'no data',
// for cases such as folds, off the target image, piece too small, etc. .
// Otherwise if the pixel has value N, the closest point is N-10
//
ffmap map1;
uint8 *ids = (uint8*)malloc(w*h*sizeof(uint8));
memset( ids, 0, w * h );

int next_id = 10;  //ids of regions will start with 10

for(int j=0; j<cr.size(); j++) {
    // split the connected regions into sub-regions of limited size.
    printf("\n-----Starting region %d-----\n", j);
    vector<vertex> ControlPoints;
    vector<triangle> tris;
    CreateOutline(cr[j].pts, ControlPoints, tris);

    vector<ConnRegion> subs(tris.size());

    // assign each point to one of the sub-regions
    for(int k=0; k<cr[j].pts.size(); k++) {
        Point pt = cr[j].pts[k];
	int reg = BestTriangle( tris, ControlPoints, pt );
        subs[reg].pts.push_back(pt);
	subs[reg].B.L = min(subs[reg].B.L, int(pt.x));
	subs[reg].B.B = min(subs[reg].B.B, int(pt.y));
	subs[reg].B.R = max(subs[reg].B.R, int(pt.x));
	subs[reg].B.T = max(subs[reg].B.T, int(pt.y));
        }
    // Now, try to match each constituent triangle by itself.  Save the results in vector tcorr and dx,dy of vector 'subs'
    vector<double> tcor;
    for(int i=0; i<subs.size(); i++) {
	printf("\n %d pts, x=[%d %d], y = [%d %d]\n", subs[i].pts.size(), subs[i].B.L, subs[i].B.R, subs[i].B.B, subs[i].B.T);
        vector<double> vals;
        for(int k=0; k<subs[i].pts.size(); k++) {
	    int ix = int(subs[i].pts[k].x);
	    int iy = int(subs[i].pts[k].y);
            vals.push_back(raster[ix + w*iy]);
	    }
	Normalize(vals);
        //WriteDebugImage(subs[i].pts, vals, 100 + 100*j + i);
	double	c = CorrPatchToImage( subs[i].dx, subs[i].dy, subs[i].pts, vals, image2, 0, 0, 4000, true );
        tcor.push_back(c);
	printf(" c = %f, deltas are %f %f\n", c, subs[i].dx, subs[i].dy);
	}
    // Now, find a 'concensus' transform - one that many triangles agree on
    int best = -1;
    int cbest = 0;  // Number of close matchs to the best one found
    double corr_of_best = 0.0;   // and correlation of that one
    for(int i=0; i<subs.size(); i++) {
	int close = 0; // count how many are close to this offset (use 15 pixels)
        for(int k=0; k<subs.size(); k++) {
	    if( fabs(subs[i].dx-subs[k].dx) < 15.0 && fabs(subs[i].dy - subs[k].dy) < 15.0 )
		close++;
	    }
        if( close > cbest || (close == cbest && tcor[i] > corr_of_best) ) {  // more close matches than any other point?
            best = i;
            cbest = close;
            corr_of_best = tcor[i];
            }
        }
    printf("OK, best correlation is sub-region %d, %d are close, correlation %f\n", best, cbest, tcor[best]);
    TForm tr_guess; // our best guess at the transform relating
    if( cbest >= 2 && tcor[best] > 0.25 ) {  // just guessing here.  Need two triangle match, better corr > 0.25
        vector<double> vals;                             // should make a routine
        for(int k=0; k<subs[best].pts.size(); k++) {
	    int ix = int(subs[best].pts[k].x);
	    int iy = int(subs[best].pts[k].y);
            vals.push_back(raster[ix + w*iy]);
	    }
	Normalize(vals);
	c=ImproveCorrelation(subs[best].pts, vals, image2, subs[best].dx, subs[best].dy, tr_guess, flog);
        printf("Best guess transform: "); tr_guess.PrintTransform();

	// now using the guess, keep only those points that map to the below image.  (this could be off
	// by 10 pixels or so, so we'll need to fix it later, but it's a good guess
	vector<Point> newpts;
        for(int i=0; i<cr[j].pts.size(); i++) {
	    Point p(cr[j].pts[i]);
	    tr_guess.Transform( p );
	    if( p.x >= 0 && p.x < below[0].w && p.y >= 0 && p.y < below[0].h )
		newpts.push_back(cr[j].pts[i]);
	    }
        // some quality checks should go here
        printf("About %d pixels map to the image below\n", newpts.size());
	//MakeSubregions(newpts, -cb, cb, -cb, cb, subs);  // divide again, could be very different
	CreateOutline(newpts, ControlPoints, tris);
        // for each triangle, create a matrix that transforms to barycentric coordinates
        for(int k=0; k<tris.size(); k++) {
            vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
            printf("Tri: (%d %d) (%d %d) (%d %d)\n", v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
	    double a[3][3];
	    a[0][0] = v0.x; a[0][1] = v1.x; a[0][2] = v2.x;
	    a[1][0] = v0.y; a[1][1] = v1.y; a[1][2] = v2.y;
	    a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
	    Invert3x3Matrix(tris[k].a, a);
	    }
        for(int k=0; k<tris.size(); k++) {  // write in matlab form
            vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
            printf("x=[%d %d %d %d]; y=[%d %d %d %d]; plot(x,y); hold on;\n", v0.x, v1.x, v2.x, v0.x,
	     v0.y, v1.y, v2.y, v0.y);
	    }
        subs.clear(); subs.insert(subs.begin(), tris.size(), ConnRegion());

	// assign each point to one of the sub-regions
	vector<vector<double> >lambda;
	for(int k=0; k<newpts.size(); k++) {
	    Point pt = newpts[k];
	    int reg = BestTriangle( tris, ControlPoints, pt );
	    subs[reg].pts.push_back(pt);
	    subs[reg].B.L = min(subs[reg].B.L, int(pt.x));
	    subs[reg].B.B = min(subs[reg].B.B, int(pt.y));
	    subs[reg].B.R = max(subs[reg].B.R, int(pt.x));
	    subs[reg].B.T = max(subs[reg].B.T, int(pt.y));
	    // Map the point into barycentric coordinates for its triangle
            //vertex v0 = ControlPoints[tris[reg].v[0]], v1 = ControlPoints[tris[reg].v[1]], v2 = ControlPoints[tris[reg].v[2]];
            //printf(" pt (%f %f) in tri: (%d %d) (%d %d) (%d %d)\n", pt.x, pt.y, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
	    double l[3];
	    l[0] = tris[reg].a[0][0]*pt.x + tris[reg].a[0][1]*pt.y + tris[reg].a[0][2]*1.0;
	    l[1] = tris[reg].a[1][0]*pt.x + tris[reg].a[1][1]*pt.y + tris[reg].a[1][2]*1.0;
	    l[2] = tris[reg].a[2][0]*pt.x + tris[reg].a[2][1]*pt.y + tris[reg].a[2][2]*1.0;
            //printf("prop %f %f %f\n", l[0], l[1], l[2]);
            vector<double>lam(ControlPoints.size(), 0.0);
            for(int m=0; m<3; m++)
		lam[tris[reg].v[m]] = l[m];
            //for(int m=0; m<lam.size(); m++)
		//printf("%.4f ", lam[m]);
            //printf("\n");
            lambda.push_back(lam);
            ids[int(pt.x) + w*int(pt.y)] = next_id + reg;
	    }
        next_id += tris.size();
        // optimize the mesh.
        vector<double> spv;  // create a normalized array of source pixel values
        for(int k=0; k<newpts.size(); k++) {
	    int ix = int(newpts[k].x);
	    int iy = int(newpts[k].y);
            spv.push_back(raster[ix + w*iy]);
	    }
	Normalize(spv);
        // next, transform the control points to the target frame
        vector<Point> cpts; // make a real-valued copy
	for(int k=0; k<ControlPoints.size(); k++)
	    cpts.push_back(Point(ControlPoints[k].x, ControlPoints[k].y));
        vector<Point> orig = cpts;
        tr_guess.Transform( cpts );
	// call the optimizer to find a better spot for the control points
	ImproveControlPts(cpts, lambda, spv, image2, 4096, 4096, flog);
        // compute a transform, and a center, for each triangle
        for(int k=0; k<tris.size(); k++) {
            int i0 = tris[k].v[0];
            int i1 = tris[k].v[1];
            int i2 = tris[k].v[2];
	    // Now, find a transformation that maps ORIG into the new cpts
	    // first, create a transform that maps a unit right triangle to the original pts
	    TForm o(orig[i1].x-orig[i0].x, orig[i2].x-orig[i0].x, orig[i0].x,
	     orig[i1].y-orig[i0].y, orig[i2].y-orig[i0].y, orig[i0].y);
	    // now one that maps the final optimized control points to the unit right triangle
	    TForm c(cpts[i1].x-cpts[i0].x, cpts[i2].x-cpts[i0].x, cpts[i0].x,
	     cpts[i1].y-cpts[i0].y, cpts[i2].y-cpts[i0].y, cpts[i0].y);
	    // now, to get from the original to the final, apply o^-1, then c;
	    TForm oi,t;
	    InvertTrans(oi, o);
	    //TForm temp;
	    MultiplyTrans(t, c, oi);
	    t.PrintTransform();
            map1.transforms.push_back(t);
            double sumx = orig[i0].x + orig[i1].x + orig[i2].x;
            double sumy = orig[i0].y + orig[i1].y + orig[i2].y;
            map1.centers.push_back(Point(sumx/3.0, sumy/3.0));
	    }
	}
    }

for(int i=0; i<map1.centers.size(); i++) {
    printf(" center %f %f :", map1.centers[i].x, map1.centers[i].y);
    map1.transforms[i].PrintTransform();
    }

// Read any existing mapping files
vector<ffmap> mapv;
ReadMappingFile(mapv, argv[4], flog);

// If we already have data for this image, replace it.
map1.fname = argv[1];
map1.mname = "map.tif";
int k;
for(k=0; k<mapv.size(); k++) {
    if( mapv[k].fname == argv[1] ) {  // found the name; use old map file name
        map1.mname = mapv[k].mname;
	mapv[k] = map1;
        break;
	}
    }
if( k == mapv.size() ) {    // we did not find it, so add it
    int m = -1;
    for(int j=0; j<mapv.size(); j++) {
	size_t p2 = mapv[j].mname.find_last_of('.')-1;    // one before the last dot
	size_t p1 = mapv[j].mname.find_last_of('.', p2);
        int n = atoi(mapv[j].mname.substr(p1+1,p2-p1).c_str() );
        m = max(m, n);
        }
    char file_name[256], copy[256];
    strcpy(copy, argv[2]);
    char *p = strrchr(copy, '.');
    *p = 0;  // prune off last .tif
    sprintf(file_name,"%s.map.%d.tif", copy, m+1);
    map1.mname = file_name;
    mapv.push_back(map1);
    }

// Now write it out
FILE	*fm = FileOpenOrDie( argv[4], "w", flog );

fprintf(fm, "<local_maps>\n");
for(k=0; k<mapv.size(); k++) {
    // write the image name
    fprintf(fm, " <entry overlap=\"%s\" map=\"%s\">\n", mapv[k].fname.c_str(), mapv[k].mname.c_str() );
    // write the transforms
    for(int j=0; j<mapv[k].transforms.size(); j++) {
        TForm a = mapv[k].transforms[j];
	fprintf(fm, "  <map id=\"%d\" transform=\"%f %f %f %f %f %f\"/>\n",
         j+10, a.t[0], a.t[1], a.t[3], a.t[4], a.t[2]*sc, a.t[5]*sc);
	}
    fprintf(fm, " </entry>\n");
    }
fprintf(fm, "</local_maps>\n");
fclose(fm);

// If scale > 1, then we need to make a full scale map to match the image
// If scale = 1, this is just a very complex copy.
uint8 *save_ids = ids;
save_ids = (uint8*)RasterAlloc( fullw * fullh * sizeof(uint8) );
for(int ix=0; ix<w; ix++) {
    for(int iy=0; iy<h; iy++) {
	int jx = ix*sc;
	int jy = iy*sc;
	for(int dx=0; dx<sc; dx++)
	    for(int dy=0; dy<sc; dy++)
		save_ids[jx+dx + fullw*(jy+dy)] = ids[ix+w*iy];
	}
    }

// write this out as (argv[2]-.tif).map.%d.tif so we can look at it.
Raster8ToTif8( map1.mname.c_str(), save_ids, fullw, fullh );
RasterFree(save_ids);








/* --------- */
/* ABOverlay */
/* --------- */

// Here we write out a diagnostic image.  This requires several steps:
xmin = -1000, xmax = 3000, ymin = -1000, ymax = 3000; // just for testing
w2 = uint32(xmax-xmin+1), h2=uint32(ymax-ymin+1);   // These describe the comparison image

npixels2 = w2 * h2;

raster2 = (uint32*)RasterAlloc( npixels2 * sizeof(uint32) );      // monochrome - specify plane or white

rt  = (float*)RasterAlloc( npixels2 * sizeof(float) );            // a temporary layer.  It's made of floats
									 // since we will write partial pixel
									 // values here, and so should not truncate

for(int k=0; k<npixels2; k++) {
    raster2[k] = 0xFF000000;       // set the alpha channel
    rt[k] = 0.0;
    }

// for the image above.  For each pixel, look at the map. If the map is defined, write the data into that spot
// of the image.  If it is not defined, then find the closest center and use that transform.  If using that transform
// you fall within the below image, then do nothing (we wish to see through).  Otherwise copy it (we wish to see the
// surroundings).
//
for(int i=0; i<npixels; i++) {
    int iy = i / w;
    int ix = i - w * iy;
    int mv = ids[i];  // get the map value
    Point p(ix, iy);
    if (mv >= 10)            // transform is defined, see where it goes
	map1.transforms[mv-10].Transform( p );
    else { // not strictly defined...
	double dbest = 1.0E30;
        int best = -1;
        for(int q=0; q < map1.centers.size(); q++) {
	    double dx = ix-map1.centers[q].x;
            double dy = iy-map1.centers[q].y;
	    double d = dx*dx + dy*dy;   // no need for sqrt, it's monotonic
	    if( d < dbest ) {
		dbest = d;
	        best = q;
	        }
	    }
        // now see, if using the best transform, if it's in the target image
        map1.transforms[best].Transform( p );
        if( p.x >= 0.0 && p.x < below[0].w && p.y >= 0.0 && p.y < below[0].h )
	    continue; // do not want to set it.
	}
    // so now p is the transformed point.  Subtract (xmin, ymin) to get coordinate in image
    p.x -= xmin;
    p.y -= ymin;
    DistributePixel( p.x, p.y, raster[ix+w*iy], rt, w2, h2 );
}
// transfer temp to the RED channel.  Invert and stretch contrast
for(int i=0; i<npixels2; i++) {
    if( rt[i] > 0.0 ) {
        int pix = int(127 + (127-rt[i])*2);
		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;
        raster2[i] |= pix;
        }
    }
// now just copy the image below to the GREEN channel.  This will need to be changed to use maps
for(int i=0; i<npixels; i++) {
    int iy = i / w;
    int ix = i - w * iy;
    int pix = 127 + (127 - below[0].raster[i]) * 2;  // invert and contrast stretch
	if( pix < 0 )
		pix = 0;
	else if( pix > 255 )
		pix = 255;
    ix -= xmin;
    iy -= ymin;
    if( ix >= 0 && ix < w2 && iy >= 0 && iy < h2 )  // it's in the output image
	raster2[ix + w2*iy] |= (pix << 8);
    }

// Write it out
Raster32ToTifRGBA( "comp.tif", raster2, w2, h2 );

exit(-1);  // temporary








// here we align between layers
TForm tr_guess; // on exit from the loop, this will be our best guess at the transform relating
                // pixels in 'above' -> 'below'
for(int j=0; j<below.size(); j++) {
    for(int k=0; k<above.size(); k++) {
        // transfers above's bounding box to below's coordinate system
        printf("orig size %d %d\n", above[k].w, above[k].h);
        Point bb1(0.0,0.0);       Point bb2(0.0 , above[k].h-1);  // corners of image
	Point bb3(above[k].w-1,0.0); Point bb4(above[k].w-1, above[k].h-1);  // corners of image
	printf("orig pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	above[k].tr.Transform( bb1 ); above[k].tr.Transform( bb2 ); // into global space
	above[k].tr.Transform( bb3 ); above[k].tr.Transform( bb4 );
	printf(" global pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	below[j].Inverse.Transform( bb1 ); below[j].Inverse.Transform( bb2 );
	below[j].Inverse.Transform( bb3 ); below[j].Inverse.Transform( bb4 );
	printf("transformed pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	double xll = min(bb1.x, min(bb2.x, min(bb3.x, bb4.x)));  // probably not needed,
	double yll = min(bb1.y, min(bb2.y, min(bb3.y, bb4.y)));  // since transforms are similar
	double xur = max(bb1.x, max(bb2.x, max(bb3.x, bb4.x)));
	double yur = max(bb1.y, max(bb2.y, max(bb3.y, bb4.y)));
	double oxmax = min(double(below[j].w-1), xur);
	double oxmin = max(0.0, xll);
	double oymax = min(double(below[j].h-1), yur);
	double oymin = max(0.0, yll);
	if( oxmax > oxmin && oymax > oymin ) {
	    printf( "\nOverlap: xrange [%.1f %.1f]; yrange [%.1f %.1f]\n", oxmin, oxmax, oymin, oymax);
	    printf(" for file '%s'\n", below[j].fname.c_str() );
            Point centerb((oxmin+oxmax)/2.0, (oymin+oymax)/2.0);
            Point centera(centerb);
            below[j].tr.Transform( centera );
            above[k].Inverse.Transform( centera );
            printf("Initial guess: below (%f %f) -> above (%f %f)\n", centerb.x, centerb.y, centera.x, centera.y);
	    // Compute the mean and std deviation of the first using only the 'real' pixels -
	    // those that are >0, and hence not part of the background.
	    printf("Aligning images %d and %d\n", j, k);
	    int npixels = below[j].w*below[j].h;
            MeanStd m;
	    for(int i=0; i<npixels; i++) {
		if( below[j].raster[i] > 0 )
		    m.Element(below[k].raster[i]);
		}
	    printf("below: %d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);
	    double mean2, std2;
	    m.Stats(mean2, std2);  // find existing statistics
	    printf("Of the target image,  mean= %f and std dev = %f\n", mean2, std2);

            // Make a 4Kx4K normalized copy, with image in lower left 2Kx2K, so FFT will work.
            // Background pixels map to zero
	    vector<double> image2(4096*4096, 0.0);
	    for(int i=0; i<npixels; i++) {
		int y = i / below[j].w;
		int x = i - below[j].w * y;
		double pix = below[j].raster[i];
		if (pix == 0)    // background pixels
		    pix = mean2;  // will be set to the mean value (0)
		image2[x + 4096*y] = (pix-mean2)/std2;
		}

	    // Now, for the other image, do this
	    // Find the points in image that map onto the first picture
	    vector<Point> pts;
	    vector<double> vals;
	    for(int i=0; i<above[k].w*above[k].h; i++) {
		int y = i / above[k].w;
		int x = i - above[k].w * y;
		Point pt(x,y);            // Coordinates in picture j
                Point d(pt.x-centera.x, pt.y-centera.y);
                double dist = sqrt(d.x*d.x + d.y*d.y);
		if( dist <= RADIUS ) { // for now
		    //printf("pt: x,y=%f %f\n", pt.x, pt.y);
		    pts.push_back(pt);
		    vals.push_back(above[k].raster[i]);
		    //vals.push_back(below[j].raster[i+6050]);  // just testing
		    }
                //else  // so we can look at it...
                    //above[k].raster[i] = 0;
		}
	    printf("There are %d overlap pixels\n", pts.size() );
	    if( pts.size() < 10000 ) {
		printf("Not enough pixels.\n");
		fprintf(flog,"Not enough pixels.\n");
		exit( 42 );
		}
	    Normalize(vals);
	    double	dx, dy;
	    double	c = CorrPatchToImage( dx, dy, pts, vals, image2, 0, 0, 4000, true );
	    printf(" c = %f, deltas are %f %f\n", c, dx, dy);
	    c=ImproveCorrelation(pts, vals, image2, dx, dy, tr_guess, flog);
            Point cb(centera.x+dx, centera.y+dy);
            Point save(centera);
            printf("Disk around spot (%f %f) in above should map to (%f %f) in below\n",
              centera.x, centera.y, cb.x, cb.y);
            below[j].tr.Transform( cb );  // move to global space
            above[k].tr.Transform( centera );
            above[k].tr.AddXY( cb.x - centera.x, cb.y - centera.y );
            // try this instead.   We have a better transformation tr, that maps B to A
            // so A^-1(B(pt)) = tr;
            MultiplyTrans(above[k].tr, below[j].tr, tr_guess);
	    InvertTrans(above[k].Inverse, above[k].tr);  // Recompute j's inverse
            }
	}
    }
printf("New and improved: "); above[0].tr.PrintTransform();

// Now we have a rough alignment, close to exact at one spot (but not necessarily in the center).
// Compute 5 spots and their best local transforms.
for(int j=0; j<below.size(); j++) {
    for(int k=0; k<above.size(); k++) {
        // transfers above's bounding box to below's coordinate system
        printf("orig size %d %d\n", above[k].w, above[k].h);
        Point bb1(0.0,0.0);       Point bb2(0.0 , above[k].h-1);  // corners of image
	Point bb3(above[k].w-1,0.0); Point bb4(above[k].w-1, above[k].h-1);  // corners of image
	printf("orig pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	above[k].tr.Transform( bb1 ); above[k].tr.Transform( bb2 ); // into global space
	above[k].tr.Transform( bb3 ); above[k].tr.Transform( bb4 );
	printf(" global pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	below[j].Inverse.Transform( bb1 ); below[j].Inverse.Transform( bb2 );
	below[j].Inverse.Transform( bb3 ); below[j].Inverse.Transform( bb4 );
	printf("transformed pts %f %f %f %f %f %f %f %f\n", bb1.x, bb1.y, bb2.x, bb2.y, bb3.x, bb3.y, bb4.x, bb4.y);
	double xll = min(bb1.x, min(bb2.x, min(bb3.x, bb4.x)));  // probably not needed,
	double yll = min(bb1.y, min(bb2.y, min(bb3.y, bb4.y)));  // since transforms are similar
	double xur = max(bb1.x, max(bb2.x, max(bb3.x, bb4.x)));
	double yur = max(bb1.y, max(bb2.y, max(bb3.y, bb4.y)));
	double oxmax = min(double(below[j].w-1), xur);
	double oxmin = max(0.0, xll);
	double oymax = min(double(below[j].h-1), yur);
	double oymin = max(0.0, yll);
	if( oxmax > oxmin && oymax > oymin ) {
	    printf( "\nOverlap: xrange [%.1f %.1f]; yrange [%.1f %.1f]\n", oxmin, oxmax, oymin, oymax);
	    printf(" for file '%s'\n", below[j].fname.c_str() );
            Point centerb((oxmin+oxmax)/2.0, (oymin+oymax)/2.0);
            Point centera(centerb);
            below[j].tr.Transform( centera );
            above[k].Inverse.Transform( centera );
            printf("Initial guess: below (%f %f) -> above (%f %f)\n", centerb.x, centerb.y, centera.x, centera.y);
	    // Compute the mean and std deviation of the first using only the 'real' pixels -
	    // those that are >0, and hence not part of the background.
	    printf("Aligning images %d and %d\n", j, k);
	    int npixels = below[j].w*below[j].h;
            MeanStd m;
	    for(int i=0; i<npixels; i++) {
		if( below[j].raster[i] > 0 )
		    m.Element(below[k].raster[i]);
		}
	    printf("below: %d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);
	    double mean2, std2;
	    m.Stats(mean2, std2);  // find existing statistics
	    printf("Of the target image,  mean= %f and std dev = %f\n", mean2, std2);

            // Make a 4Kx4K image, with the data of at most 2Kx2K in the lower left corner.
            // This prevents wrap-around, and is needed for the FFT matching.
            // The data is normalized, and the background pixels map to zero
	    vector<double> image2(4096*4096, 0.0);
	    for(int i=0; i<npixels; i++) {
		int y = i / below[j].w;
		int x = i - below[j].w * y;
		double pix = below[j].raster[i];
		if (pix == 0)    // background pixels
		    pix = mean2;  // will be set to the mean value (0)
		image2[x + 4096*y] = (pix-mean2)/std2;
		}
            vector<double>image3;
            if( below[j].scale > 1 ) { // then we have a copy with more than 2Kx2K resolution
                int np3 = npixels*below[j].scale*below[j].scale;
		MeanStd m3;
		for(int i=0; i<np3; i++) {
		    if( below[j].original[i] > 0 )
			m3.Element(below[j].original[i]);
		    }
		double mean3, std3;
		m3.Stats(mean3, std3);  // find existing statistics
		printf("Of the full resolution target image,  mean= %f and std dev = %f\n", mean3, std3);
		image3.push_back(0.0);
                image3.insert(image3.begin(), 4096*4096-1, 0.0); // is there a better way?
                int bigw = below[j].w*below[j].scale;
		for(int i=0; i<np3; i++) {
		    int y = i / bigw;
		    int x = i - bigw * y;
		    double pix = below[j].original[i];
		    if (pix == 0)    // background pixels
			pix = mean3;  // will be set to the mean value (0)
		    image3[x + 4096*y] = (pix-mean3)/std3;
		    }
		}

	    // now, find a number of points, in above's coordinates,
	    // that extend as far as possible but still map into 'below'
            map1.centers.push_back(centera);
	    FindMoreSpots(below, j, above, k, map1.centers, centera.x, centera.y);
            for(int q=0; q<map1.centers.size(); q++)
		printf("Center at (%f %f)\n", map1.centers[q].x, map1.centers[q].y);

	    // Now, for the other image, do this
	    // Find the points in image that map onto the first picture
	    for(int m=0; m<map1.centers.size(); m++) {
                printf("----------- Looking at (%f %f) ------------\n", map1.centers[m].x, map1.centers[m].y);
	        vector<Point> pts;
	        vector<double> vals;
                CircleFromRaster(above[k].raster, above[k].w, above[k].h, map1.centers[m], RADIUS, pts, vals);
		printf("There are %d overlap pixels\n", pts.size() );
		if( pts.size() < 10000 ) {
		    printf("Not enough pixels.\n");
		    fprintf(flog,"Not enough pixels.\n");
		    exit( 42 );
		    }
		Normalize(vals);
                vector<Point> copy = pts;
                tr_guess.Transform( copy );
                // now we would expect dx, dy to be small
		double	dx, dy;
		double	c = CorrPatchToImage( dx, dy, copy, vals, image2, 0, 0, 4000, true );
		printf(" c = %f, deltas are %f %f\n", c, dx, dy);
		TForm tr = tr_guess;
                tr.t[2] += dx;
                tr.t[5] += dy;
		c=ImproveCorrelation(pts, vals, image2, double(BIG), double(BIG), tr, flog);   // BIG->use starting tr
                if( above[k].scale > 1 ) { // improve further, if we can
                    printf("--- Full res improve:"); tr.PrintTransform();
                    int sc = above[k].scale;
                    Point nc(map1.centers[m].x*sc, map1.centers[m].y*sc);
		    CircleFromRaster(above[k].original, above[k].w*sc, above[k].h*sc, nc, RADIUS*sc, pts, vals);
                    TForm bigtr = tr;
                    bigtr.t[2] *= sc; bigtr.t[5] *= sc;  // scale up to full image
		    c=ImproveCorrelation(pts, vals, image3, double(BIG), double(BIG), bigtr, flog);   // BIG->use starting tr
                    printf("--- Full res result:"); bigtr.PrintTransform();
                    bigtr.t[2] /= sc; bigtr.t[5] /= sc;  // scale back to 2K size
                    tr = bigtr;
		    }
                map1.transforms.push_back(tr);
                fprintf(flog," %f", c);
		}
	    // We have a better transformation tr, that map1 points in 'above' to points in 'below'
	    // so B^-1(A(pt)) = tr;  so A = B*tr;
	    // use the 0th transform since that's centered.
	    MultiplyTrans(above[k].tr, below[j].tr, map1.transforms[0]);
	    InvertTrans(above[k].Inverse, above[k].tr);  // Recompute j's inverse
            }
	}
    }
// Original may have been bigger than the 2Kx2k we use for processing; scale the centers and transforms back up
sc = above[0].scale;
for(k=0; k<map1.centers.size(); k++) {
    map1.centers[k].x *= sc;
    map1.centers[k].y *= sc;
    map1.transforms[k].t[2] *= sc;
    map1.transforms[k].t[5] *= sc;
    }

// create a map image.  For every pixel in A, tell the closest control point.  0-9  reserved for 'no data',
// for cases such as folds, off the target image, piece too small, etc. .
// Otherwise if the pixel has value N, the closest point is N-10
//
fullw = above[0].w*sc;
fullh = above[0].h*sc;
ids = (uint8*)malloc(fullw*fullh*sizeof(uint8));
for(int ix=0; ix<fullw; ix++) {
    for(int iy=0; iy<fullh; iy++) {
	double dbest = 1.0E30;
        int best = -1;
        for(int q=0; q < map1.centers.size(); q++) {
	    double dx = ix-map1.centers[q].x;
            double dy = iy-map1.centers[q].y;
	    double d = dx*dx + dy*dy;   // no need for sqrt, it's monotonic
	    if( d < dbest ) {
		dbest = d;
	        best = q;
	        }
	    }
        // now see, if using the best transform, if it's in the target image
        Point p(ix, iy);
        map1.transforms[best].Transform( p );
        if( p.x >= 0.0 && p.x < (below[0].w)*sc && p.y >= 0.0 && p.y < (below[0].h)*sc )
	    ids[ix+fullw*iy] = best+10;
        else
            ids[ix+fullw*iy] = 0;  // does not map to anything
	}
    }

// Read any existing mapping files
ReadMappingFile(mapv, argv[4], flog);

// If we already have data for this image, replace it.
map1.fname = argv[1];
map1.mname = "map.tif";
for(k=0; k<mapv.size(); k++) {
    if( mapv[k].fname == argv[1] ) {  // found the name; use old map file name
        map1.mname = mapv[k].mname;
	mapv[k] = map1;
        break;
	}
    }
if( k == mapv.size() ) {    // we did not find it, so add it
    int m = -1;
    for(int j=0; j<mapv.size(); j++) {
	size_t p2 = mapv[j].mname.find_last_of('.')-1;    // one before the last dot
	size_t p1 = mapv[j].mname.find_last_of('.', p2);
        int n = atoi(mapv[j].mname.substr(p1+1,p2-p1).c_str() );
        m = max(m, n);
        }
    char file_name[256], copy[256];
    strcpy(copy, argv[2]);
    char *p = strrchr(copy, '.');
    *p = 0;  // prune off last .tif
    sprintf(file_name,"%s.map.%d.tif", copy, m+1);
    map1.mname = file_name;
    mapv.push_back(map1);
    }

// Now write it out
fm = FileOpenOrDie( argv[4], "w", flog );

fprintf( fm, "<local_maps>\n" );
for(k=0; k<mapv.size(); k++) {
    // write the image name
    fprintf(fm, " <entry overlap=\"%s\" map=\"%s\">\n", mapv[k].fname.c_str(), mapv[k].mname.c_str() );
    // write the transforms
    for(int j=0; j<mapv[k].transforms.size(); j++) {
        TForm a = mapv[k].transforms[j];
	fprintf(fm, "  <map id=\"%d\" transform=\"%f %f %f %f %f %f\"/>\n",
         j+10, a.t[0], a.t[1], a.t[3], a.t[4], a.t[2], a.t[5]);
	}
    fprintf(fm, " </entry>\n");
    }
fprintf(fm, "</local_maps>\n");
fclose(fm);

// write this out as (argv[2]-.tif).map.%d.tif so we can look at it.
Raster8ToTif8( map1.mname.c_str(), ids, fullw, fullh );

// Now write the composite image for visulization.  This is done with the 2Kx2K copies, so scale transforms back
for(int k=0; k<map1.centers.size(); k++) {
    map1.centers[k].x /= sc;
    map1.centers[k].y /= sc;
    map1.transforms[k].t[2] /= sc;
    map1.transforms[k].t[5] /= sc;
    }

xmin = -1000; xmax = 3000; ymin = -1000; ymax = 3000; // just for testing
w2 = uint32(xmax-xmin+1), h2=uint32(ymax-ymin+1);   // These describe the comparison image
npixels2 = w2 * h2;

int nvpix = 0; // number of valid pixels

raster2 = (uint32*)RasterAlloc( npixels2 * sizeof(uint32) );              // monochrome - specify plane or white
uint32 *craster = (uint32*)RasterAlloc( npixels2 * sizeof (uint32) );      // color
uint32 *craster2 = (uint32*)RasterAlloc( npixels2 * sizeof (uint32) );      // color
for(int k=0; k<npixels2; k++) {
    raster2[k] = 0xFF000000;
    craster[k] = 0xFF000000;
    craster2[k] = 0xFF000000;
    }

// for creating a composite image, we want each piece to map global to above (A^-1)
for(int q=0; q<map1.centers.size(); q++) {
    TForm temp;
    MultiplyTrans(temp, below[0].tr, map1.transforms[q]); // map1 A to global
    InvertTrans(map1.transforms[q], temp);
    }

ffmap empty;
nvpix = CreateCompositeImage(int(xmin), int(ymin), w2, h2, above, raster2, 'R', craster, ids, map1.transforms);
nvpix = CreateCompositeImage(int(xmin), int(ymin), w2, h2, below, raster2, 'G', craster2, NULL, empty.transforms);
// write the synthesised circles into the BLUE plane
for(int j=0; j<map1.centers.size(); j++) {
    // map1.transforms[k] contains mapping global->above.  We need the inverse here
    TForm temp;
    InvertTrans(temp, map1.transforms[j]);
    for(int k=0; k<2*3.14159*RADIUS; k++) { // write points
	double ang=double(k)/RADIUS;
        for(int r = RADIUS-1; r<=RADIUS+1; r++) {
	    double dx = r*cos(ang);
	    double dy = r*sin(ang);
	    Point p(map1.centers[j].x + dx, map1.centers[j].y + dy);
	    temp.Transform( p );
	    int ix = int(p.x+0.5 - xmin);
	    int iy = int(p.y+0.5 - ymin);
	    if( ix >= 0 && ix < w2 && iy >= 0 && iy < h2 )
		raster2[ix + w2*iy] |= (255 << 16);
	    }
	}
    }
printf("Synthesized image has %d valid pixels\n", nvpix);

// write out the synthesized image in black and white
//uint8 *nb = (uint8*)malloc( w2 * h2 * sizeof(uint8) );
//for(int j=0; j<w2*h2; j++)
    //nb[j] = raster2[j] & 0xFF;

// Write it out
Raster32ToTifRGBA( "comp.tif", raster2, w2, h2 );

if( nvpix < 100000 ) {
    printf("Not enough pixels (%d) in the synthesised image\n", nvpix);
    fprintf(flog,"Not enough pixels (%d) in the synthesised image\n", nvpix);
    exit( 42 );
    }

// Write it out in color to check alignment
Raster32ToTifRGBA( "ccomp.tif", craster, w2, h2 );

if( nvpix < 100000 ) {
    printf("Not enough pixels (%d) in the synthesised image\n", nvpix);
    fprintf(flog,"Not enough pixels (%d) in the synthesised image\n", nvpix);
    exit( 42 );
    }
RasterFree(craster);

fprintf(flog,"\n"); fclose(flog);
return 0;
}
