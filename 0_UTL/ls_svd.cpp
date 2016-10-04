
#include	"jh6.h"

#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#include	<stdlib.h>

#include	<vector>
using namespace std;


double ** matrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;
    m = (double **)malloc((nrh-nrl+1)*sizeof(double *));
    if (m == NULL) {
        printf("Odd...\n");
    exit(-1);
    }
    m -= nrl;
    for(i=nrl; i<=nrh; i++) {
    m[i] = (double *)malloc((nch-ncl+1)*sizeof(double));
    m[i] -= ncl;
    }
    return m;
}

double * NRvector(int nl, int nh)
{
    double *v;
    v = (double *)malloc((nh-nl+1)*sizeof(double));
    return v-nl;
}
void nrerror(const char *c)
{
    printf("%s\n", c);
    exit(1);
}

void free_vector(double *v, int nl, int nh)
{
    free((char *)(v+nl));
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
for(int i=nrl; i<=nrh; i++)
    free((char*)(m[i] + ncl));  // free each row
free((char *)(m+nrl));         // then free the array of pointers
}

double fabs(double f);
double fmin(double a, double b) { return a < b ? a : b;}
double fmax(double a, double b) { return a > b ? a : b;}
int     max(int    a, int    b) { return a > b ? a : b;}
//int     max(int    a, int    b);

// Solve A*x = b by least squares, using SVD from Numerical Recipes, modified for
// double precision and C++.  Also return the approximate solution bprime.
int LeastSquaresBySVD(vector<vector<double> > &aa, vector<double> &xx, vector<double> &bb, vector<double> &bp)
{
// step 1, find how may equations and how many unknowns.
int m = aa.size();  // number of equations
int realm;          // real number of equations, in case there are not enough and we need to add some
int n=0;            // number of unknowns
for(int i=0; i<aa.size(); i++)
    n = max(n, aa[i].size());
printf("there are %d equations in %d unknowns\n", m, n);
realm = m;
if (m < n)  /* make up some fake rows */
    m = n;

// OK, make up the stuff that NR needs
double **   a = matrix(1,m,1,n);
double **save = matrix(1,m,1,n);
double **   v = matrix(1,n,1,n);
double *    w = NRvector(1,n);
double *    x = NRvector(1,n);	/* get final solution here */
double *    b = NRvector(1,m);	/* constant vector of equation */

for(int i=1; i<=m; i++)
    for(int j=1; j<=n; j++)
    save[i][j] = a[i][j] = 0.0;
for(int i=0; i<aa.size(); i++) {
    for(int j=0; j<aa[i].size(); j++){
        save[i+1][j+1] = a[i+1][j+1] = aa[i][j];
        //if (aa[i][j] != 0.0)
             //printf("%f X%d  ", aa[i][j], j+1);
        }
    b[i+1] = bb[i];
    //printf("= %f\n", bb[i]);
    }

svdcmp(a,m,n,w,v);

// remove singular values
int nsingular = 0;
{
   double wmax = 0.0, wmin;
       for(int i=1; i<=n; i++) if (w[i] > wmax) wmax = w[i];
       //for(int i=1; i<=n; i++) {
       //printf("w[%d] = %g\n",i,w[i]);
       //for(int j=1; j<=n; j++)
           //printf("%10g ", v[i][j]);
       //printf("\n");
       //}
       wmin = wmax * 1.0E-10;
       for(int i=1; i<=n; i++)
      if (w[i] < wmin) {
         printf("Singular value! Index %d\n",i);
             w[i] = 0.0;
             nsingular++;
         }
    }

// Now do the back substitution
svbksb(a,w,v,m,n,b,x);
//printf("\n");
//printf("Solution: ");
xx.clear();
for(int i=1; i<=n; i++) {
    //printf("%11g ", x[i]);
    xx.push_back(x[i]);
    }
//printf("\n");
bp.clear();
for(int i=1; i<=realm; i++) {
    double sum = 0.0;
    for(int j=1; j<=n; j++)
       sum += save[i][j]*x[j];
    bp.push_back(sum);
    }
// free everything we allocated
free_matrix(a,   1,m,1,n);
free_matrix(save,1,m,1,n);
free_matrix(v,   1,n,1,n);
free_vector(w,   1,n);
free_vector(x,   1,n);
free_vector(b,   1,m);
return nsingular;
}
