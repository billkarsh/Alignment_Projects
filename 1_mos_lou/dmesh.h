
#ifndef _DMESH_
#define _DMESH_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <stack>
#include <complex>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <fftw3.h>
using namespace std;

// this declaration used the "vector" type, so must go after "using"

typedef complex<double> CD;


class tform {
    public:
    tform(double a, double b, double c, double d, double e, double f)
	{t[0]=a; t[1]=b; t[2]=c; t[3]=d; t[4]=e; t[5]=f;}
    tform() { t[0]=1.0; t[1] = 0.0; t[2] = 0.0; t[3] = 0.0; t[4]=1.0; t[5]=0.0;}
    void ToMatlab()  {double b=t[1], c=t[2], d=t[3], e=t[4]; t[1]=d; t[2]=b; t[3]=e; t[4]=c;} // from order here to one Matlab prefers
    void FromMatlab(){double b=t[2], c=t[4], d=t[1], e=t[3]; t[1]=b; t[2]=c; t[3]=d; t[4]=e;} // and vice versa
    double det() {return t[0]*t[4] - t[1]*t[3];}
    double t[6];
    };


class ffmap {  // set of maps from oen file to another
    public:
    string fname;  // name of the image mapped to
    string mname;  // name of mapping file
    vector<Point> centers;     // centers of sub-maps
    vector<tform> transforms;  // transform for each portion
    };

// classes for creating a mesh
class vertex {
  public:
    vertex(){x = y = dir = orig = 0;}
    vertex(int x0, int y0, char d0){x = x0; y = y0; dir = d0;}
    vertex(const vertex &v){x = v.x; y = v.y; dir = v.dir; orig = v.orig;}
    friend bool operator==(const vertex &a, const vertex &b);
    friend bool operator<(const vertex &a, const vertex &b);
    int x,y;  // location of vertex
    char dir; // direction to next vertex
    int orig;  // where was this in the original array?
    };

class lineseg {
  public:
    lineseg(vertex &v1, vertex &v2){v[0] = v1; v[1] = v2;}
    lineseg(int x1, int y1, int x2, int y2){v[0].x = x1; v[0].y=y1; v[1].x = x2, v[1].y=y2;}
    lineseg(){v[0].x = 0; v[0].y = 0; v[1].x = 0; v[1].y = 0;}
    bool operator==(const lineseg &a) const {return a.v[0] == this->v[0] && a.v[1] == this->v[1];}
    vertex v[2];
  };

class triangle {
   public:
    int v[3];  // indices into the vector of control points
    double a[3][3];   // matrix that turns a point (x,y,1) into a barycentric representation (linear combination of vertices).
   };

class ConnRegion { // a connected region
  public:
   ConnRegion(){corr=0.0;}
   vector<Point> pts;  // pixels within the region
   vector<Point> all;  // all the points, since the list gets trimmed but we need it later
		       // to reconstruct the image
   double dx, dy;      // deltas to make it line up with image 2
   int xmin, xmax, ymin, ymax; // bounding box in original image
   tform trans;  // transform that maps into image2
   double corr;  // the correlation we got on this piece
   int id;       // pixel value in the mask that generated this connected region
   };

class Picture {
   public:
     Picture(){original = NULL; raster = NULL;}
     ~Picture(){if(raster!=NULL && raster!=original) free(raster);}

   void DownsampleIfNeeded();
   void MakeDoGExist( double r1, double r2, vector<CD> &cache);
   tform tr;  // the transform
   int z;     // the Z layer
   string fname; // the file name
   double dist_from_center;  // distance from center
   bool operator<(const Picture &rhs) const {return this->dist_from_center < rhs.dist_from_center;}
   uint32 w,h;     // These describe the image
   uint8 *raster;  // This is the one used in most computations (usually 2K by 2K)
   vector<uint8> DoG;  // Difference of Gaussians, if it exists
   tform Inverse;      // the inverse transform
   CD *fft_of_frame;   // the FFT of the frame.
   void MakeFFTExist( int i);
   uint8 *original;    // original raster
   int scale;          // scale factor between original and processed (typically 1, 2, 4)
   };

class Template {
   private:
   int nx, ny; // size of fft
   int M;      // number of complex values in FFT
   int x0,y0;  // origin of data.  Needs to be added to all results
   CD* fft;    // points to the FFT of the data
   public:
   Template(vector<Picture> &vp, int i, int xmin, int xmax, int ymin, int ymax);
   ~Template(){fftw_free(fft);}
   Point Match(vector<Picture> &vp, int j, int xmin, int xmax, int ymin, int ymax);
   };

class QueueEl {
  public:
   int x,y;       // we can get to (x,y)
   int cost;      // for this cost
   QueueEl(int xx, int yy, int cc){x = xx; y = yy; cost = cc;}
   bool operator<(const QueueEl &rhs) const {return this->cost > rhs.cost;} // backwards since we want lowest first
   };

class CorrCand { // a correlation candidate
  public:
    int x,y;
    double val;
    CorrCand(){x=0; y=0; val=-BIG;}
    CorrCand(int xx, int yy, double v){x = xx; y = yy; val=v;}
    bool operator<(const CorrCand &rhs) const {return this->val < rhs.val;} // sort
    };

// -------------- for storing correlations among images ---------------------------
class CorrImages;

class CorrPair {
  friend class CorrImages;
  friend CorrImages *ReadImageCorrelations(const char *fname, FILE *flog);
  public:
   CorrPair(int ii1, int ii2){i1 = ii1; i2 = ii2;}
   bool operator<(const CorrPair &rhs) const {return this->i1 < rhs.i1 || (this->i1 == rhs.i1 && this->i2 < rhs.i2);} // sort
  private:
   int i1, i2;  // image 1 and image 2
   vector<Point> p1s;  // points in image 1
   vector<Point> p2s;  // corresponding points in image 2
   };

class CorrImages {
   friend CorrImages *ReadImageCorrelations(const char *fname, FILE *flog);
   public:
    int FindImageCorr(string name1, string name2, vector<Point>&p1s, vector<Point>&p2s);
    int AddImageCorr( string name1, string name2, vector<Point>&p1s, vector<Point>&p2s);
    int WriteImageCorrelations(const char *fname);
   private:
    map<string,int> names;
    set<CorrPair> corrs;
    };

class Decimate {  // used for decimating Point lists
  public:
    int x,y;
    double v;
    bool operator<(const Decimate &rhs) const {return this->y < rhs.y || (this->y == rhs.y && this->x < rhs.x);} // sort
    };

// used for seeing how well triangles match.
class Match {
  public:
    vector<Point> pts;  // the points in image a
    vector<double> a;   // the value in image a
    vector<double> b;   // the value in the (mapped) image b
    };

void Transform(Point &p, tform &t);
void PrintTransform(FILE *f, tform &tr);
void InvertTrans(tform &inv, tform &b);
void MultiplyTrans(tform &r, tform &a, tform &b);
double FindNormCorrelation(vector<Point> &pts, vector<double> &vals, vector<Point> &p2, vector<double>&v2,
                double &dx, double &dy, int tx, int ty, int radius, FILE *flog,
		bool (*LegalRegion)(int sx, int sy, void *v), void *arg,
		bool (*LegalCounts)(int sx, int sy, void *v), void *arg2,
	        vector<CD> &ftc);
void SetDMeshParams(const char *arg);

#endif
