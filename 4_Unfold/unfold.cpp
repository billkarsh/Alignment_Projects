

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"CTemplate.h"
#include	"CQueueEl.h"
#include	"CCorrImages.h"
#include	"CCorrCand.h"

#include	"ls_svd.h"
#include	"tinyxml.h"

#include	<string.h>

#include	<algorithm>
#include	<queue>
using namespace std;


int OVERLAP = 75;  //expect about this much overlap, in pixels
int RADIUS = 250;  // correction for alignment should be found within this distance
double THRESHOLD = 0.25;  // lowest correlation considered a match


class Picture : public PicBase {

public:
   double	overlap;	// approximate amount of overlap

public:
	Picture()	{overlap = 0.0;};

	bool operator < (const Picture &rhs) const
		{return overlap > rhs.overlap;};
};


class ConnRegion_old { // a connected region
  public:
   vector<Point> pts;  // pixels within the region
   vector<Point> all;  // all the points, since the list gets trimmed but we need it later
		       // to reconstruct the image
   int dx, dy;          // deltas to make it line up with image 2
   IBox B; // bounding box in original image
   TAffine trans;  // transform that maps into image2
   double corr;  // the correlation we got on this piece
   };






// Improve the correlation, if possible, by tweaking the transform.  pts are the points
// in the original, and pts are their values.  image is a 4Kx4K matrix of doubles, already
// normalized.  dx and dy are the initial estimates of how to map the original points into
// the array image2.
// returns the best correlation obtained.
double ImproveCorrelation(vector<Point> &Plist, vector<double> &spv, vector<double>image2,
 int dx, int dy, TAffine &t, FILE *flog)
{
printf("Contains %d pixels\n", Plist.size() );
Normalize(spv);
t.t[0]=1; t.t[1]=0; t.t[2]=dx; t.t[3]=0; t.t[4]=1; t.t[5]=dy;
double best_so_far = 0.0;

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
    TAffine tbest;  // best transform found
    int bdir = -1;  // flag to tell what the best direction is

    // make steps of this size until no more improvement
    double prev = start;
    double nc = prev;     // new correlation
    do {
        prev = nc;
	double norm = sqrt(dvx*dvx + dvy*dvy);
	double sx = step*dvx/norm;  // direction sin()
	double sy = step*dvy/norm;
	printf("step is %f %f\n", sx, sy);
	TAffine t2( t );
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
	    tbest.CopyIn( t2 );
	    bdir = 1;
	    }
        }
    while(nc > prev);

    for(int rot = -1; rot<2; rot += 2) {  // try two rotations, too

    TAffine	t2, R;
    R.SetCWRot( rot * step, cog );
    t2 = R * t;

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
	    tbest.CopyIn( t2 );
	    }
	}
    // Now tried 8 directions and two rotations; pick the best if any were better
    if( bdir >= 0 ) { //we found a better transform
	 t.CopyIn( tbest );
	 }
    else // nothing was better; reduce step size
	step = step/2;
    }
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
double	c = CorrPatchToImage( dx, dy, pts, vals, image2, 0, 0, RADIUS, false );
if( pass == 2 && c < 0.25 ) {
    printf("Dubious alignment improvement: %d pts, corr=%f\n", pts.size(), c);
    printf("file 1: %s\n", vp[i].fname.c_str() );
    printf("file 2: %s\n", vp[j].fname.c_str() );
    fprintf(flog,"Dubious alignment improvement: %d pts, corr=%f\n", pts.size(), c);
    exit( 42 );
    }
printf(" deltas are %d %d\n", dx, dy);
Point zero(0.0, 0.0);
Point delta(dx, dy);      // what we want in space of image i
vp[i].tr.Transform( delta );        vp[i].tr.Transform( zero );
// make the fix in global space, since this is where the tr[2] and [5] are applied
vp[j].tr.AddXY( delta.x - zero.x, delta.y - zero.y );
printf("deltas in global image space %f %f\n", delta.x-zero.x, delta.y-zero.y);
vp[j].Inverse.InverseOf( vp[j].tr );
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
for(int k=0; k<np; k++) {
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
TAffine tf(
	scale*cos(ai-aj),
	scale*(-sin(ai-aj)),
	0.0,
	-tf.t[1],
	tf.t[0],
	0.0 );
// now transform the first point and make sure it ends up in the right place
Point p(j1.x,j1.y);
tf.Transform( p );
tf.AddXY( i1.x - p.x, i1.y - p.y );
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
vp[j].tr.CopyIn( tf );
vp[j].Inverse.InverseOf( tf );
}

// Point Pi in image i should align with point Pj in image j.  i < j.  Direction 'unconstrained'
// is not tighly constrained by the data, and needs an extra constraint to make sure it does
// not contract to 0.
void WriteEqns(FILE *feq, vector<Picture> &vp, int i, Point Pi, int j, Point Pj,
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
        printf(" Initial transform %d: ", i); vp[i].tr.TPrint();
        printf(" Initial transform %d: ", j); vp[j].tr.TPrint();
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

        //if( nc == 0 ) {
            //FindCorrPoints(vp, i, j, cpi, cpj);
            // are there any correspondence points?  If so add to list
            //if( cpi.size() > 0 )
	        //CI->Add( vp[i].fname, vp[j].fname, cpi, cpj );
            //}
        // generate equations from points, if any
        for(int k=0; k<cpi.size(); k++)
	    WriteEqns(feq, vp, i, cpi[k], j, cpj[k], eqns, rhs);
	}

    // for all except the first, add some preference for a square soln (rotation and scaling only)
    if( i != 0 ) {
        int id = (i-1)*6;
        double strength = 50.0;  // how strong is this preference?  pixel error vs matrix asymmetry
                                  // these are both true for pure rotation + scaling
        fprintf(feq, "%f X%d  %f X%d = 0\n", strength, id+1, -strength, id+5);  // a11 = a22
        fprintf(feq, "%f X%d  %f X%d = 0\n", strength, id+2,  strength, id+4);  // a12 = -a21
        vector<double>eq1(id+5,0.0);
        eq1[id] = strength; eq1[id+4] = -strength;
        eqns.push_back(eq1);  rhs.push_back(0.0);
        vector<double>eq2(id+4,0.0);
        eq2[id+1] = strength; eq2[id+3] = strength;
        eqns.push_back(eq2);  rhs.push_back(0.0);
        }
    }
fclose( feq );
//CI->Write( "Corr.xml" );
delete CI;
//system("../svdfit eqns -print");
//feq = FileOpenOrDie( "coeff", "r", flog );

vector<double>x;
vector<double>rslt;
int ns = LeastSquaresBySVD(eqns, x, rhs, rslt);
if( ns > 0 ) {
    printf("Singular values in least-squares fit!\n");
    fprintf(flog,"Singular values in least-squares fit!\n");
    exit( 42 );
    }
for(int j=1; j<vp.size(); j++) {
    vp[j].tr.CopyIn( &x[(j-1)*6] );
    vp[j].Inverse.InverseOf( vp[j].tr );
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




int main(int argc, char* argv[])
{
FILE	*flog = FileOpenOrDie( "unfold.log", "a" );

if( argc < 3 ) {
    printf("Usage: unfold <xml-file> <string>, where <string> is a substring of the file name\n");
    exit( 42 );
    }
bool WriteBack = false;   // overwrite the original?
bool FromAbove = false;   // get comparison from above?

time_t t0 = time(NULL);
char atime[32];
strcpy(atime, ctime(&t0));
atime[24] = '\0';  // remove the newline
fprintf(flog,"Unfold: %s ", atime);
for(int i=1; i<argc; i++) {
    fprintf(flog,"%s ", argv[i]);
    if( strcmp(argv[i],"-o") == 0 && i+1 < argc )
	OVERLAP = atoi(argv[i+1]);
    if( strcmp(argv[i],"-t") == 0 && i+1 < argc )
	THRESHOLD = atof(argv[i+1]);
    if( strcmp(argv[i],"-r") == 0 && i+1 < argc )
	RADIUS = atoi(argv[i+1]);
    if( strcmp(argv[i],"-f") == 0 )
	WriteBack = true;
    if( strcmp(argv[i],"-a") == 0 )
	FromAbove = true;
    }
fflush(flog);
printf("Overlap %d, radius %d, threshold %f\n", OVERLAP, RADIUS, THRESHOLD);
// Read the trakEM xml file to see what's there.
vector<Picture> vp;
TiXmlDocument doc(argv[1]);
bool loadOK = doc.LoadFile();
printf("XML load gives %d\n", loadOK);
if (!loadOK)  {
    printf("Could not open XML file '%s'\n", argv[1]);
    fprintf(flog,"Could not open XML file '%s'\n", argv[1]);
    exit( 42 );
    }
TiXmlHandle hDoc(&doc);
TiXmlElement* child;
TiXmlHandle hRoot(0);

// block: should be <trakem2>
TiXmlNode*	node=0;
node =doc.FirstChild();
if( node == NULL )
    printf("No node??\n");
child=hDoc.FirstChild("trakem2").FirstChild("t2_layer_set").FirstChild("t2_layer").ToElement();
// should always have a valid root but handle gracefully if it does
if( !child ) {
    printf("No first child??\n");
    return 42;
    }
printf("child element value %s\n", child->Value() );
for( child; child; child=child->NextSiblingElement() ) {
    //printf("Got a <t2_layer>\n");
    const char *attr  = child->Attribute("z");
    //printf("z= %s\n", attr);
    // now look through all the <t2_patch> elements in each
    TiXmlElement*	c2;
    c2 = child->FirstChildElement("t2_patch");
    for( c2; c2; c2=c2->NextSiblingElement() ) {
	//printf("got a <t2_patch>\n");
        const char *tf  = c2->Attribute("transform");
        const char *fp  = c2->Attribute("file_path");
        //printf("File is '%s'\n, transform is '%s'\n", fp, tf);
        Picture p;
        p.tr.ScanTrackEM2( tf );
        p.z = int(atof(attr)+0.5);
        p.fname = fp;
        vp.push_back(p);
	}
    }

// OK, look for a name containing argv[2], then open that...
char fname[2048];
strcpy(fname, "No match for file name");
int i;
for(i=0; i<vp.size(); i++){
    int j = vp[i].fname.find(argv[2]);
    if( j != string::npos ) {
	strcpy(fname, vp[i].fname.c_str());
	break;
	}
    }

// Read the first file
uint32 w, h;
uint8* raster = Raster8FromTif(fname, w, h, flog);
size_t npixels = w*h;

// See if we need to process it...
//
// Step 1 - compute the mean and std deviation using only the 'real' pixels -
// those that are >0, and hence not part of the background.
    vector<double> v;
    for(int i=0; i<npixels; i++)
        if( raster[i] > 0 )
            v.push_back(raster[i]);
    printf("%d real pixels, %f percent\n", v.size(), v.size()*100.0/npixels);
    double mean, std;
    Stats(v, mean, std);  // find existing statistics
    printf("Of the image points, mean= %f and std dev = %f\n", mean, std);

// Make two normalized copies.  Do not copy the pixels on the 'edge', since they are averages
// of the background and 'real' pixels, and can falsely assume intermediate values.
// Make 2 copies since the first (v) will be destroyed during connected region evaluation.
    v.clear();
    vector<double> vorig;
    for(int i=0; i<npixels; i++) {
	int y = i / w;
	int x = i - w * y;   // do not copy pixels on edge or next to background pixels
        int pix = raster[i];
        if (x-1 <  0 || (x-1 >= 0 && raster[i-1] == 0) ) pix = 0;
        if (x+1 >= w || (x+1 <  w && raster[i+1] == 0) ) pix = 0;
        if (y-1 <  0 || (y-1 >= 0 && raster[i-w] == 0) ) pix = 0;
        if (y+1 >= h || (y+1 <  h && raster[i+w] == 0) ) pix = 0;
        v.push_back((pix-mean)/std);
        vorig.push_back((pix-mean)/std);
        }

// Now find the connected regions
    int nbig = 0;
    vector<ConnRegion_old> cr;
    //printf("Suggest %f\n", (183-mean)/std);
    double thresh = 4.0;
    int start = 0;
    // find all connected regions of 'reasonable' values.  The background should be about
    // -4 or -5 on this scale.
    for(int k=0;;k++){
        int i;
        for(i=start; i<npixels; i++)
             if( fabs(v[i]) < thresh )
		break;
        if( i >= npixels )
            break;
        // found at least one pixel.  Find all connected ones.
        //printf("found v[%d] = %f\n", i, v[i]);
        start = i;  // next time, start here...
        ConnRegion_old c;
        cr.push_back(c);
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
                cr[k].pts.push_back(p);
                v[j] = -5.0;
                if (x-1 >= 0) st.push(j-1);
                if (x+1 < w)  st.push(j+1);
                if (y-1 >= 0) st.push(j-w);
                if (y+1 < h)  st.push(j+w);
                }
            }
        printf("Connected region of %d pixels\n", cr[k].pts.size());
        BBoxFromPoints( cr[k].B, cr[k].pts );
        printf("region size is [%d %d] in x, [%d %d] in y\n",
         cr[k].B.L, cr[k].B.R, cr[k].B.B, cr[k].B.T);
        if( cr[k].pts.size() > 10000 )
	    nbig++;
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

    // If there is only one connected region big enough to process, or one dominates, then no sense going further
    if( nbig == 1 || bc_size > 0.98*w*h ) {
          printf("One connected region dominates: %d of %d, %f%%\n", bc_size, w*h,
	   double(bc_size)/w/h*100.0);
          fprintf(flog,"One connected region dominates: %d of %d, %f%%\n", bc_size, w*h,
	   double(bc_size)/w/h*100.0);
          exit( 42 );
          }

// find the width of the fold.  Make a map, color by connected parts.  Then for each pixel on the edge, look for
// the nearest colored pixel.
vector<double> ds; // all the widths, measured from each pixel on the edge
int *map = (int*)malloc( w * h * sizeof(int) );
memset( map, 0, w * h * sizeof(int) );

for(int j=0; j<cr.size(); j++) // for each region
    for(int k=0; k<cr[j].pts.size(); k++) {
       int ix = int(cr[j].pts[k].x);
       int iy = int(cr[j].pts[k].y);
       map[ix+w*iy] = (j == bc) ? 1 : 2; // biggest region is 1, all else is two
       }
for(int ix=30; ix<w-30; ix++) {       // don't look at edge since opposite edge may not be there.
    for(int iy=30; iy<h-30; iy++) {
	int j = ix + w*iy;
        //printf("ix , iy, j, map[j] = %d %d %d %d\n", ix, iy, j, map[j]);
        if( map[j] == 1 && (map[j-1]==0 || map[j+1]==0 || map[j+w] == 0 || map[j-w] == 0) ) {
	    // we found a pixel on the edge of a fold
	    double dmin = BIG;
	    for(int x=ix-30; x<=ix+30; x++) {
		for(int y=iy-30; y<=iy+30; y++) {
		     if( map[x + w*y] == 2 ) {
			int d = (x-ix)*(x-ix) + (y-iy)*(y-iy);
			if( d < dmin )
			    dmin = d;
			}
		    }
		}
            dmin = sqrt(dmin);
            if( dmin < 50 ) {  // only record if a neighbor was found
	        //printf("edge at %d, %d.  Nearest %f\n", ix, iy, dmin );
                ds.push_back(dmin);
                }
	    }
	}
    }
double avgf, stdf;
Stats(ds, avgf, stdf);
printf("Average fold width is %f, std dev %f\n", avgf, stdf);
//exit( 42 );
//
//We only want stuff 'close' to a fold.  Here we call close 200 pixels.
//Push the fold pixels into the queue, set rest to +BIG
priority_queue<QueueEl> q;
for(int k=0; k<w*h; k++) {
    int y = k / w;
    int x = k - w * y;
    if( map[k] == 0 && x>0 && y>0 && x<w-1 && y<h-1 ) {
        q.push(QueueEl(x,y,0));
        }
    else
        map[k] = BIG;
    }
// now do the expansion
while(!q.empty()) {
    int x= q.top().x, y=q.top().y, c=q.top().cost;
    //printf("%d entries: Can get to %d %d for cost %d\n", q.size(), x, y, c);
    q.pop();
    int z = x + w*y;
    //printf("z, map[z] %d %d\n", z, map[z]);
    if( x > 0 && map[z-1] > c+1 ) {
	map[z-1] = c+1;
        q.push(QueueEl(x-1,y,c+1));
        }
    if( x < w-1 && map[z+1] > c+1 ) {
	map[z+1] = c+1;
        q.push(QueueEl(x+1,y,c+1));
        }
    if( y > 0 && map[z-w] > c+1 ) {
	map[z-w] = c+1;
        q.push(QueueEl(x,y-1,c+1));
        }
    if( y < h-1 && map[z+w] > c+1 ) {
	map[z+w] = c+1;
        q.push(QueueEl(x,y+1,c+1));
        }
    }
vector<int>histo(w+h,0);
for(int k=0; k<w*h;k++)
    histo[map[k]]++;
int cumulative = 0;
for(int k=0; k<300; k++){
    cumulative += histo[k];
    //printf("%7d pixels of cost %3d, cumulative %7d\n", histo[k], k, cumulative);
    }

// now, for every connected region, get rid of the pixels with cost >200, if any
for(int k=0; k<cr.size(); k++) {
    // make a copy before we trim the list; we need the full list later to re-create image
    cr[k].all.clear();
    int j=0; // the write pointer
    for(int i=0; i<cr[k].pts.size(); i++) {
        cr[k].all.push_back(cr[k].pts[i]);
        int ix = int(cr[k].pts[i].x);
        int iy = int(cr[k].pts[i].y);
        if( map[ix+w*iy] <= 200 )
	    cr[k].pts[j++] = cr[k].pts[i];
	}
    int del = cr[k].pts.size() - j;
    printf("Connected region %d, deleting %d points\n", k, del);
    cr[k].pts.erase(cr[k].pts.end()-del, cr[k].pts.end());
    }

// If we get here, we need to process the image.
// Now go through all the images from the XML file, looking for the ones above and below
vector<Picture> above;
vector<Picture> below;
Point pt(1020.0,1020.0); // a temporary point - BLOOT - should center in image
vp[i].tr.Transform( pt );
printf("Location in global space is %f %f\n", pt.x, pt.y);
double bd = 1.0E30;  // best distance
int bp = -1;      // best picture
for(int j=0; j<vp.size(); j++) {
    vp[j].overlap = -1.0;
    Point p(1020.0, 1020.0);  //bloot
    vp[j].tr.Transform( p );   // map this to global space
    double oxmax = min(p.x+1020, pt.x+1020);
    double oxmin = max(p.x-1020, pt.x-1020);
    double oymax = min(p.y+1020, pt.y+1020);
    double oymin = max(p.y-1020, pt.y-1020);
    if( vp[j].z == vp[i].z+1 || vp[j].z == vp[i].z-1 ) { // adjacent layer
        if( oxmax > oxmin && oymax > oymin ) {
            vp[j].overlap = (oxmax-oxmin)*(oymax-oymin);
            if( vp[j].z > vp[i].z )
		above.push_back(vp[j]);
	    else
                below.push_back(vp[j]);
            printf(
             "\nOverlap, layer delta %d; xrange [%.1f %.1f]; yrange [%.1f %.1f]; area %.1f\n",
             vp[j].z - vp[i].z, oxmin, oxmax, oymin, oymax, vp[j].overlap );
            printf(" for file '%s'\n", vp[j].fname.c_str() );
	    double d = sqrt(pow(p.x-pt.x,2.0) + pow(p.y-pt.y, 2.0));
            if( d < bd ) {
	        bd = d;
	        bp = j;
	        }
	    }
	}
    }

// If the user specified a second file, use that.  Otherwise use best match
char fname2[2048];
if( argc >= 4 )
    strcpy(fname2, argv[3]);
else {
    if( bp == -1 ) {
        fprintf(flog, "No overlapping candidate found\n");
        printf("No overlapping candidate found\n");
	exit( 42 );
	}
    printf("OK, will open '%s' on layer %d\n\n", vp[bp].fname.c_str(), vp[bp].z);
    vp[i].tr.TPrint(  stdout, "First  " );
    vp[bp].tr.TPrint( stdout, "Second " );
    strcpy(fname2, vp[bp].fname.c_str());
    }

// sort the vector by overlap.  Read each of them in
sort(above.begin(), above.end());
sort(below.begin(), below.end());
for(int j=0; j<above.size(); j++) {
    printf("Above: z=%d, overlap=%f\n", above[j].z, above[j].overlap);
    above[j].raster = Raster8FromTif(above[j].fname.c_str(), above[j].w, above[j].h, flog);
    above[j].Inverse.InverseOf( above[j].tr );
    }
if( above.size() > 1 ) {
    printf("--------- Aligning image above ------\n");
    NewAlign(above, flog);
    //ImproveAlignment(above, flog,1);
    //printf("------------------------------ Pass 2 --------------------------------------\n");
    //ImproveAlignment(above, flog,2);
    }
for(int j=0; j<below.size(); j++) {
    printf("Below: z=%d, overlap=%f\n", below[j].z, below[j].overlap);
    below[j].raster = Raster8FromTif(below[j].fname.c_str(), below[j].w, below[j].h, flog);
    below[j].Inverse.InverseOf( below[j].tr );
    }
if( below.size() > 1 ) {
    printf("--------- Aligning image below ------\n");
    NewAlign(below, flog);
    //ImproveAlignment(below, flog,1);
    //printf("------------------------------ Pass 2 --------------------------------------\n");
    //ImproveAlignment(below, flog,2);
    }

uint32 w2, h2;   // These describe the comparison image
uint8*	raster2;
size_t npixels2;

// One way to get the image - commented out for now
    //raster2 = Raster8FromTif( fname2, w2, h2 );
    //npixels2 = w2*h2;
//
//The other way = construct a master image.
//
vector<Picture> *ref;
const char *pstr;
ref = FromAbove ? &above : &below;
if( (*ref).size() == 0 ) {
    printf("No overlapping image on plane %s\n", FromAbove?"above":"below");
    fprintf(flog,"No overlapping image on plane %s\n", FromAbove?"above":"below");
    exit( 42 );
    }

int nvpix = 0; // number of valid pixels
w2 = w;
h2 = h;
npixels2 = w2*h2;
raster2 = (uint8*)RasterAlloc( npixels2 * sizeof(uint8) );
uint32 *craster = (uint32*)RasterAlloc( npixels2 * sizeof(uint32) );
for(int j=0; j<npixels2; j++) {
    raster2[j] = 0;
    craster[j] = 0xFF000000;
    int y = j/w;
    int x = j-w*y;
    Point pt(x,y);
    vp[i].tr.Transform( pt );  // transform to global coordinates
    //printf("in global space %f %f, %d candidates\n", pt.x, pt.y, (*ref).size() );
    bool first = true; // for composite image, will take first match, since they
                       // are sorted by size.  For color image, combine all
    for(int k=0; k<(*ref).size(); k++) {
	Point p2(pt.x, pt.y);
	(*ref)[k].Inverse.Transform( p2 );
        //printf("Image %d: %f %f\n", p2.x, p2.y);
	// is the point within the image?
	if( p2.x >= 0 && p2.x < (*ref)[k].w-1 && p2.y >= 0 && p2.y < (*ref)[k].h-1 ) {

	   double val =
	   InterpolatePixel( p2.x, p2.y, raster, (*ref)[k].w );

	    int pix = int(val + 0.5);  // rounding
	    //printf("set [%d,%d] to %f\n", x, y, val);
	    if( first ) {
	        raster2[x+w*y] = pix;
            first = false;
	        nvpix++;
                }
            int color = k % 3;
            // invert the pixels for the color image
            pix = 255 - pix;
	    craster[x+w*y] |= (pix << (8*color));
            }
	}
    }
printf("Synthesized image has %d valid pixels\n", nvpix);
// write out the synthesized image in black and white
uint8 *nb = (uint8*)malloc( w * h * sizeof(uint8) );
for(int j=0; j<w*h; j++)
    nb[j] = raster2[j];

// Write it out
Raster8ToTif8( "comp.tif", nb, w, h );

if( nvpix < 100000 ) {
    printf("Not enough pixels (%d) in the synthesised image\n", nvpix);
    fprintf(flog,"Not enough pixels (%d) in the synthesised image\n", nvpix);
    exit( 42 );
    }

// Write it out in color to check alignment
Raster32ToTifRGBA( "ccomp.tif", craster, w, h );

if( nvpix < 100000 ) {
    printf("Not enough pixels (%d) in the synthesised image\n", nvpix);
    fprintf(flog,"Not enough pixels (%d) in the synthesised image\n", nvpix);
    exit( 42 );
    }
RasterFree(craster);

bool invert = getenv("INVERT") != NULL;
if (invert)  {
    printf("Inverting images...\n");
    for(int i=0; i<w*h; i++) {
        if( raster[i] )
			raster[i] = 255 - raster[i];
	}
    for(int i=0; i<w2*h2; i++) {
        if( raster2[i] )
			raster2[i] = 255 - raster2[i];
        }
    }

// Compute the mean and std deviation of image 2 using only the 'real' pixels -
// those that are >0, and hence not part of the background.
    vector<double> v2;
    for(int i=0; i<npixels; i++)
        if( raster2[i] > 0 )
            v2.push_back(raster2[i]);
    printf("%d real pixels, %f percent\n", v2.size(), v2.size()*100.0/npixels2);
    double mean2, std2;
    Stats(v2, mean2, std2);  // find existing statistics
    printf("Of the second image points, mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy.  Background pixels map to zero
    vector<double> image2(4096*4096, 0.0);
    for(int i=0; i<npixels2; i++) {
        int y = i / w2;
        int x = i - w2 * y;
        double pix = raster2[i];
        if (pix == 0)    // background pixels
            pix = mean2;  // will be set to the mean value (0)
        image2[x + 4096*y] = (pix-mean2)/std2;
        }
// Now we find the folds, dust, and tears
    //
    // Start by finding the correlation of the biggest region
    vector<double> vals;
    for(int k=0; k<cr[bc].pts.size(); k++) {
        int j = int(cr[bc].pts[k].x) + w*int(cr[bc].pts[k].y);
	vals.push_back(vorig[j]);
        }

    double	_dx, _dy;
    CorrPatchToImage( _dx, _dy, cr[bc].pts, vals, image2, 0, 0, 4000, false );
    cr[bc].dx = (int)_dx;
    cr[bc].dy = (int)_dy;

    cr[bc].corr =
      ImproveCorrelation(cr[bc].pts, vals, image2, cr[bc].dx, cr[bc].dy, cr[bc].trans, flog);
    fprintf(flog, "%6.3f ", cr[bc].corr);
    for(int k=0; k<cr.size(); k++) {
        if( k == bc )
            continue;
        if(  cr[k].pts.size() > 10000 ) {
            vals.clear();
	    for(int l=0; l<cr[k].pts.size(); l++) {
		int j = int(cr[k].pts[l].x) + w*int(cr[k].pts[l].y);
                if( fabs(vorig[j]) > 10 )
                     printf("Bogus %f %d\n", v[j], j);
		vals.push_back(vorig[j]);
		}
             printf("Before corr %d %d\n", cr[k].pts.size(), vals.size());

             double	_dx, _dy;
             CorrPatchToImage( _dx, _dy, cr[k].pts, vals, image2, cr[bc].dx, cr[bc].dy, 200, false );
             cr[k].dx = (int)_dx;
             cr[k].dy = (int)_dy;

             cr[k].corr =
	       ImproveCorrelation(cr[k].pts, vals, image2, cr[k].dx, cr[k].dy, cr[k].trans, flog);
             fprintf(flog, "%6.3f ", cr[k].corr);
             }
        else { // too small to correlate; just set dx and dy to 0
             cr[k].dx = 0;
             cr[k].dy = 0;
             }
        }
    int oxmin = BIG, oxmax = -BIG, oymin = BIG, oymax = -BIG;
    // copy each chunk to the output, using the relationship to bc
    TAffine inv_big;
    inv_big.InverseOf( cr[bc].trans );
    printf(
    "Inverse transform is %9.4f %9.4f %10.2f\n"
    "                     %9.4f %9.4f %10.2f\n",
    inv_big.t[0], inv_big.t[1], inv_big.t[2], inv_big.t[3], inv_big.t[4], inv_big.t[5] );
    for(int k=0; k<cr.size(); k++) {
        if (cr[k].pts.size() <= 10000)  // if too small, don't bother
	    continue;
        int ox = cr[k].dx - cr[bc].dx;
        int oy = cr[k].dy - cr[bc].dy;
        printf("Will copy chunk %d with offset [%d,%d]\n", k, ox, oy);
        oxmin = min(oxmin, cr[k].B.L + ox);
        oxmax = max(oxmax, cr[k].B.R + ox);
        oymin = min(oymin, cr[k].B.B + oy);
        oymax = max(oymax, cr[k].B.T + oy);
        }
    printf("New image size is [%d %d] in X, [%d %d] in Y\n", oxmin, oxmax, oymin, oymax);

    // allocate a buffer
    int new_w = oxmax - oxmin + 1;
    int new_h = oymax - oymin + 1;
    vector<double> out_image(new_w*new_h, 0.0);
    vector<int> who(new_w*new_h, -1);  // which fragment wrote this pixel?
    int nMultiWrite=0;                     // number of pixels with 2 or more writers
    //
    // Copy the data
    //
    for(int k=0; k<cr.size(); k++) { // for each fragment
        if (cr[k].pts.size() <= 10000)  // if too small, don't bother
	    continue;
	int ox = cr[k].dx - cr[bc].dx;   // the offsets
	int oy = cr[k].dy - cr[bc].dy;
        printf("Copy chunk %d with offset [%d,%d]\n", k, ox, oy);
        TAffine tf = inv_big * cr[k].trans;
        printf(
        "True transform is %9.4f %9.4f %10.2f\n"
        "                  %9.4f %9.4f %10.2f\n",
        tf.t[0], tf.t[1], tf.t[2], tf.t[3], tf.t[4], tf.t[5] );
        for(int l=0; l<cr[k].all.size(); l++) {
	    int x0 = int(cr[k].all[l].x);  // pixel coordinates in original buffer
	    int y0 = int(cr[k].all[l].y);
	    int x1 = x0 + ox - oxmin;      // pixel coordinates in new buffer
	    int y1 = y0 + oy - oymin;
            double nx = tf.t[0]*x0 + tf.t[1]*y0 + tf.t[2] - oxmin;
            double ny = tf.t[3]*x0 + tf.t[4]*y0 + tf.t[5] - oymin;
            //if( k != bc )
                //printf("x1, y1 %d %d,   nx, ny %8.2f %8.2f\n", x1, y1, nx, ny);
            if( x1 < 0 || x1 >= new_w || y1 < 0 || y1 >= new_h )
		printf("Very odd %d %d (%d %d) \n", x1, y1, new_w, new_h);
            // may need to invert the Y
            // invert the Y.  There is something odd going on, since reading a file with
            // the tiff libary, then writing the exact same buffer out, results in a Y reversed
            // file as shown by gimp.
            if( nx >= 0.0 && nx < new_w-1 && ny >= 0.0 && ny < new_h-1 ) {
				x1 = int(nx);  // will round down, since numbers are positive
				y1 = int(ny);
				int base = x1 + new_w * y1;
				if( who[base] != k && who[base] != -1 )
					nMultiWrite++;
				who[base] = k; // last person to write this pixel
				DistributePixel( nx, ny, raster[x0 + w*y0], out_image, new_w, new_h );
            }
	    //nb[x1+(new_h-1-y1)*new_w] = raster[x0 + w*y0];
	    }
         }

    // Make one more pass over the image.  We compute two output images - one with the
    // fold filled in, and one subtracted from the reference image.
    //
    vector<double>diffs;
    // For image subtraction purposes, find the mean and std dev of the non-zero pixels
    // of the output image;  Use 'diffs' as a temporary vector
    for(int k=0; k<new_w*new_h; k++) {
	if( out_image[k] != 0 )
		diffs.push_back(out_image[k]);
        }
    double out_mean, out_std;
    Stats(diffs, out_mean, out_std);
    diffs.clear();

    for(int k=0; k<new_w*new_h; k++) {
	int y = k/new_w;       // coordinates in the buffer
	int x = k-new_w*y;
	x += oxmin;            // convert the coordinates to the space of the biggest chunk
	y += oymin;            // (also output, since the biggest one always has an identity
			       // transform.
	Point pt(x,y);
	cr[bc].trans.Transform( pt );  // now in the space of image2
	double val = 0.0; // value of pixel in reference image
	// does this fall within the reference image?
	if( pt.x >= 0 && pt.x < w2-1 && pt.y >= 0 && pt.y < h2-1 ) { // ok, it does

		val = InterpolatePixel( pt.x, pt.y, image2, 4096 );

		if( !out_image[k] )
			out_image[k] = 127 + 25 * val;
     }
        // since the transformed image is already 0 mean and 1 std dev, just adjust it to
        // to the same standard deviation, then subtract it.
	diffs.push_back(out_image[k] - val * out_std);
        }

    // round them back to integers
    nb = (uint8*)malloc(sizeof(uint8)*new_w*new_h);
    uint8 *nd = (uint8*)malloc(sizeof(uint8)*new_w*new_h);
    for(int j=0; j<new_w*new_h; j++) {
       nb[j] = int(out_image[j]+0.5);
       nd[j] = int(diffs[j]+0.5);
       }

// Write it out the unfolded image as 'u.tif'
	Raster8ToTif8( "u.tif", nb, new_w, new_h );



    // Write out the residual image image as 'res.tif'
	Raster8ToTif8( "res.tif", nd, new_w, new_h );

    printf("nMultiWrite = %d\n", nMultiWrite);
    if( nMultiWrite > 3*new_w ) {
	printf("Lots (%d) of overlap pixels.  This seems bad\n", nMultiWrite);
	fprintf(flog,"Lots (%d) of overlap pixels.  This seems bad\n", nMultiWrite);
	exit( 42 );
	}
// If we get here, everything looks OK.  So write the results out...
// Write the unfolded image (again) as 'u32.tif', in the floating point format, or
// if we are overwriting, put it back where we got it....(but with .new extension)
    if( WriteBack )
	strcat(fname,".new");    // add .new to original name
    else
	strcpy(fname,"u32.tif"); // just write a local file

    Raster8ToTifFlt( fname, nb, new_w, new_h );

// Write the information to the file as floats

    fprintf(flog,"\n"); fclose(flog);
    return 0;
}
