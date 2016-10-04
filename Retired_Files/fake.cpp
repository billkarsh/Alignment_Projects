

#include	"File.h"
#include	"Maths.h"
#include	"TAffine.h"

#include	<math.h>


// Random numbers, Numerical recipes with superficial changes
#define SQR(u) ((u)*(u))






struct Ranq1 {
    unsigned long long v;
    Ranq1(unsigned long long j): v(4101842887655102017LL) {
    v ^= j;
    v = int64();
    }
    inline unsigned long long int64() {
    v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
    return v * 2685821657736338717LL;
    }
    inline double doub(){return 5.42101086242752217E-20 * int64();}
    };

struct NormalDev : Ranq1 {
    double mu, sig;
    NormalDev(double mmu, double ssig, unsigned long long i) :
    Ranq1(i), mu(mmu), sig(ssig){}
    double dev() {
    double u,v,x,y,q;
        do {
        u = doub();
        v = 1.7156*(doub() - 0.5);
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = SQR(x) + y*(0.19600*y-0.25472*x);
    } while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
    return mu + sig*v/u;
    }
    };


int main(int argc, char **argv)
{
printf("Welcome to the fake alignment problem generator\n");
printf("Enter the mosaic size - x by y: ");
int nx, ny;
scanf("%d %d", &nx, &ny);
printf("Enter the number of sections (layers): ");
int nz;
scanf("%d", &nz);
printf("Enter image size and overlap): ");
int imagesize, overlap;
scanf("%d %d", &imagesize, &overlap);

const double TILE_OFF = 10.0;          // how far off from perfect tiling
const double BS_ANGLE = 0.2;           // angle between sections


vector<TAffine> tfs;
NormalDev n(0.0, 1.0, 17);  // create a generator of normally distributed random numbers
for(int z=0; z<nz; z++) {
    double theta =  z == 0 ? 0.0 : BS_ANGLE*n.dev(); // first layer at angle 0 by definition
    double c = cos(theta);
    double s = sin(theta);
    // for each layer, create a mosaic
    for(int y=0; y<ny; y++) {
    for(int x=0; x<nx; x++) {
        double x0 = x*(imagesize-overlap);
        double y0 = y*(imagesize-overlap);
            if( !(z == 0 && y == 0 && x == 0) ) {
            x0 += TILE_OFF*n.dev();  // except for very first tile, add random offset
            y0 += TILE_OFF*n.dev();
        }
            double x1 = c*x0 - s*y0;
            double y1 = s*x0 + c*y0;
            TAffine t(c, -s, x1, s, c, y1);
            tfs.push_back(t);
        }
        }
    }
vector<TAffine> inv(tfs.size());  // create an array of inverse transforms
// also print out the right answer, so we can check our work
FILE	*fp = FileOpenOrDie( "correct", "w" );

for(int k=0; k<tfs.size(); k++) {
    tfs[k].TPrint();
    tfs[k].TPrint( fp );
    inv[k].InverseOf( tfs[k] );
    inv[k].TPrint();
    }
fclose(fp);

vector<int> uses(tfs.size(),0);

fp = FileOpenOrDie( "pts", "w" );

const double MARGIN = 3.0;  // BLOOT - should compute
Ranq1 r(112);
printf("nx, ny= %d %d\n", nx, ny);
int npts = 10*(nx*ny-1); // find this many points to stitch a plane together.  If nx and ny both one, this is not possible
for(int z=0; nx*ny > 1 && z<nz; z++) {
    // Now pick random points.  Invert with each transform, and see if they fall inside.  Pick those with 2 or more
    int base = nx*ny*z;  // first transform on this layer
    for(int i=0; i<npts; ) {
       double gx = (nx+2*MARGIN)*imagesize*r.doub() - MARGIN*imagesize;
       double gy = (ny+2*MARGIN)*imagesize*r.doub() - MARGIN*imagesize;
       //printf("Same plane %d, i=%d, gx,gy=%f %f\n", z, i, gx, gy);
       int nin = 0;  // number of images point is inside of
       for(int j = 0; j < nx*ny; j++) {
        Point p(gx,gy);
            inv[base+j].Transform( p );
            nin += (0 <= p.x && p.x < imagesize-1 && 0 <= p.y && p.y <= imagesize-1);
            }
       if( nin == 2 ) {
           fprintf(fp,"POINT ");
           for(int j = 0; j < nx*ny; j++) {
            Point p(gx,gy);
            inv[base+j].Transform( p );
            if( 0 <= p.x && p.x < imagesize-1 && 0 <= p.y && p.y <= imagesize-1 ) {
                 fprintf(fp, "%d %f %f ", base+j, p.x, p.y);
                     uses[base+j]++;
             }
                }
            fprintf(fp,"\n");
        i++;
        }
        }
    }
// Now do the points between the layers
const int N_inter_plane = 20;
const double INTER_ERROR = 3.0;

for(int z=0; z<nz-1; z++) {
   for(int npts = 0; npts < N_inter_plane; ) {
       double gx = (nx+2*MARGIN)*imagesize*r.doub() - MARGIN*imagesize;  // pick a point in global space
       double gy = (ny+2*MARGIN)*imagesize*r.doub() - MARGIN*imagesize;  // allows somewhat negative values
       //printf("z=%d, gx,gy=%f %f\n", z, gx, gy);

       vector<Point> savept;  // save points here
       vector<int>   savet;   // corresponding image
       int t0 = nx*ny*z;  // first tile on layer z
       //printf("start looking\n");
       for(int t=t0; t<t0+nx*ny && savet.size() == 0; t++) {
        Point p(gx,gy);
        inv[t].Transform( p );
        int inz = (0 <= p.x && p.x < imagesize-1 && 0 <= p.y && p.y <= imagesize-1);
            //printf("inz=%d\n", inz);
        if( inz ) {
        savept.push_back(p);
        savet.push_back(t);
        }
        }
       for(int t=t0+nx*ny; t < t0+2*nx*ny && savet.size() == 1; t++) {
        Point p(gx,gy);
        inv[t].Transform( p );
        int inz1 = (0 <= p.x && p.x < imagesize-1 && 0 <= p.y && p.y <= imagesize-1);
            //printf("inz1=%d\n", inz1);
        if( inz1 ) {
        savept.push_back(p);
        savet.push_back(t);
        }
        }
       // if it maps into both layer z and z+1, keep it
       if( savet.size() == 2 ) {
            if( savet[0]/nx/ny == savet[1]/nx/ny ) {
            printf("Bogons! %d %d, nx,ny %d %d\n", savet[0], savet[1], nx, ny);
                return 1;
                }
        // Add some error to the correspondence, then print it
        Point p(savept[1].x + INTER_ERROR*n.dev(), savept[1].y + INTER_ERROR*n.dev());
        fprintf(fp, "POINT %d %f %f %d %f %f\n", savet[0], savept[0].x, savept[0].y, savet[1], p.x, p.y);
        npts++;
            uses[savet[0]]++; uses[savet[1]]++;
        }
    }
    }
fclose(fp);
for(int i=0; i<tfs.size(); i++) {
    if( uses[i] <= 3 ) {
    printf("Transform %d is used %d times\n", i, uses[i]);
        tfs[i].TPrint();
        }
    }
double C = N_inter_plane;  // three points can be taken out exactly by the fit.
double e = INTER_ERROR*INTER_ERROR*2.0;
double Np = (nx*ny == 1) ? imagesize : nx*imagesize-(nx-1)*overlap;  // number of pixels in flat image
printf("Np=%f\n", Np);
vector<double> SF(nz);
vector<double> alpha(nz);
double prod = 1.0;
for(int n=2; n<nz; n++) {
    double cSF = 1.0;   // cumulative scaling factor
    for(int j=2; j<n; j++)
    cSF = 1.0 + alpha[j]*cSF;
    SF[n] = cSF;
    alpha[n] = 1.0/(1.0 + 6*e/(Np*Np)*SF[n]);
    printf("SF[%d] = %f, alpha[%d]=%f\n", n, SF[n], n, alpha[n]);
    prod = prod * alpha[n];
    }
printf("Expect to shrink by %e\n", prod);
double a = 1.0/(1.0+6*e/(Np*Np));
printf("analytic model (small z) gives %f\n", 1.0-(1.0-a)*nz*nz/2);
// find the limiting alpha in the limit of a very large number of layers
double la = 1.0+3*e/Np/Np - sqrt(6*e/Np/Np + 9*e*e/Np/Np/Np/Np);
printf("Large layer limiting alpha=%f, strength %f\n", la, 1.0/(1.0-la) );
return 0;
}
