

/* --------------------------------------------------------------- */
/* JustTesting --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void JustTesting(vector<Point> ap, vector<double> av, vector<Point> bp, vector<double> bv, FILE *of)
{
long size = 100000;
int mult = 1;
vector<CD> ftc;
for(int i=0; i<8; i++) {
    double dx, dy;
    ftc.clear(); // no caching wanted here
    double c = CorrPatches(
                    of, false, dx, dy,
                    ap, av, bp, bv, 0, 0, 4000,
                    BigEnough, (void *)size,
                    NULL, NULL, ftc );

    fprintf(of, " JT: %4d by %4d, size %6d, corr %f, dx, dy = %9.2f, %9.2f\n",
     int(sqrt(ap.size())), int(sqrt(bp.size())),
     size, c, dx*mult, dy*mult);
    DecimateVector(ap, av, 2);
    DecimateVector(bp, bv, 2);
    size = size/4;
    mult = mult * 2;
    }
}

/* --------------------------------------------------------------- */
/* TryNewOptimizer ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void TryNewOptimizer(vector<Point> &plist, vector<double> &spv, vector<double> &image2, TAffine &t, FILE *flog)
{
// compute the bounding box
fprintf(of,"\n---------- Try new optimizer on %d points----------------\n", plist.size());
DBox	B;
BBoxFromPoints( B, plist );
fprintf(of,"region size is [%f %f] in x, [%f %f] in y\n", B.L, B.R, B.B, B.T);

// create 3 control points
vector<Point> cpts;
cpts.push_back(Point(B.L, B.B));
cpts.push_back(Point(B.R, B.B));
cpts.push_back(Point((B.L+B.R)/2, B.T));
fprintf(of,"Control points are (%f %f) (%f %f) (%f %f)\n", cpts[0].x, cpts[0].y, cpts[1].x, cpts[1].y, cpts[2].x, cpts[2].y);

// find each point as a linear combination of control points
double a[3][3];
a[0][0] = cpts[0].x; a[0][1] = cpts[1].x; a[0][2] = cpts[2].x;
a[1][0] = cpts[0].y; a[1][1] = cpts[1].y; a[1][2] = cpts[2].y;
a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
double inv[3][3];
Print3x3Matrix( of, a );
Invert3x3Matrix( inv, a );
Print3x3Matrix( of, inv );
vector<vector<double> > lambda;
for(int j=0; j<plist.size(); j++) {
    //fprintf(of," Point is (%f %f)\n", plist[j].x, plist[j].y);
    vector<double> lam;
    lam.push_back(inv[0][0]*plist[j].x + inv[0][1]*plist[j].y + inv[0][2]*1.0);
    lam.push_back(inv[1][0]*plist[j].x + inv[1][1]*plist[j].y + inv[1][2]*1.0);
    lam.push_back(inv[2][0]*plist[j].x + inv[2][1]*plist[j].y + inv[2][2]*1.0);
    //fprintf(of," lambdas are %f %f %f\n", lam[0], lam[1], lam[2]);
    lambda.push_back(lam);
    }

// Transform the control points to the target frame
vector<Point> orig = cpts;
t.Transform( cpts );

// call the optimizer
ImproveControlPts(cpts, lambda, spv, image2, 4096, flog, "new opt", 1.0);

// Now, find a transformation that maps ORIG into the new cpts
// first, create a transform that maps a unit right triangle to the original pts
TAffine o(orig[1].x-orig[0].x, orig[2].x-orig[0].x, orig[0].x,
        orig[1].y-orig[0].y, orig[2].y-orig[0].y, orig[0].y);
// now one that maps the final optimized control points to the unit right triangle
TAffine c(cpts[1].x-cpts[0].x, cpts[2].x-cpts[0].x, cpts[0].x,
        cpts[1].y-cpts[0].y, cpts[2].y-cpts[0].y, cpts[0].y);
// now, to get from the original to the final, apply o^-1, then c;
TAffine oi;
oi.InverseOf( o );
//TAffine temp;
t = c * oi;
t.TPrint( of );
}

/* --------------------------------------------------------------- */
/* ImproveCorrelation -------------------------------------------- */
/* --------------------------------------------------------------- */

// Improve the correlation, if possible, by tweaking the transform.  pts are the points
// in the original, and pts are their values.  image is a 4Kx4K matrix of doubles, already
// normalized.  dx and dy are the initial estimates of how to map the original points into
// the array image2.
// returns the best correlation obtained.
static double ImproveCorrelation(vector<Point> &Plist, vector<double> &spv, vector<double>image2,
 double dx, double dy, TAffine &t, FILE *flog)
{
fprintf(of,"Contains %d pixels\n", Plist.size() );
Normalize(spv);
if (dx != BIG)  // if dx == BIG, start improving from the transform we have
    t = TAffine(1.0, 0.0, dx, 0.0, 1.0, dy);  // otherwise, create a transform with just dx, dy

double best_so_far = 0.0;

TryNewOptimizer(Plist, spv, image2, t, flog);

// Now t is the best transform we can find.
fprintf( of,
"Best transform is %9.4f %9.4f %10.2f\n"
"                  %9.4f %9.4f %10.2f\n",
t.t[0], t.t[1], t.t[2], t.t[3], t.t[4], t.t[5] );

return best_so_far;
}


