
#include	"File.h"

#include	<math.h>
#include	<string.h>

#include	<vector>
using namespace std;


double ** matrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;
    m = (double **)malloc((nrh-nrl+1)*sizeof(double *));
    if( m == NULL ) {
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

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
void svdcmp(double **a,int m, int n, double* w, double **v);

// Solve A*x = b by least squares, using SVD from Numerical Recipes, modified for
// double precision and C++.  Also return the approximate solution bprime.
void LeastSquaresBySVD(vector<vector<double> > &aa, vector<double> &xx, vector<double> &bb, vector<double> &bp)
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
    for(int j=0; j<aa[i].size(); j++)
        save[i+1][j+1] = a[i+1][j+1] = aa[i][j];
    b[i+1] = bb[i];
    }

svdcmp(a,m,n,w,v);

// remove singular values
{
   double wmax = 0.0, wmin;
       for(int i=1; i<=n; i++) if (w[i] > wmax) wmax = w[i];
       for(int i=1; i<=n; i++) {
       printf("w[%d] = %g\n",i,w[i]);
       for(int j=1; j<=n; j++)
           printf("%10g ", v[i][j]);
       printf("\n");
       }
       wmin = wmax * 1.0E-10;
       for(int i=1; i<=n; i++)
      if( w[i] < wmin ) {
         printf("Singular value! Index %d\n",i);
             w[i] = 0.0;
         }
    }

// Now do the back substitution
svbksb(a,w,v,m,n,b,x);
printf("\n");
printf("Solution: ");
xx.clear();
for(int i=1; i<=n; i++) {
    printf("%11g ", x[i]);
    xx.push_back(x[i]);
    }
printf("\n");
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
}





int main(int argc, char **argv)
{
    double **a, *w, **v, **t, **save;
    double *x, *b;
    int i,j,  m,n;
    int realm;	/* real # of equations */
    FILE *f;
    char **eqn_names;

    m = n = 0;
    if( argc < 2 ) {
    printf("Usage: svdfit <file of equations> -print\n");
    exit(-1);
    }

    // scan through the file.  Count number of '=' to get equations,
    // biggest number after 'X' to get variables.
    // lines beginning with '#' are comments

    CLineScan	*ls = new CLineScan;

    f = FileOpenOrDie( argv[1], "r" );

    for(;;) {
        if( ls->Get( f ) <= 0 )
            break;
        if (*ls->line == '#') continue;  // skip comments
        for( char *c = ls->line; *c != '\0'; c++ ) {
            if( *c == '=' )
            m++;
        if( *c == 'X' )
        n = max(n, atoi(c+1));
        }
    }

    fclose( f );

    printf("there are %d equations in %d unknowns\n", m, n);
    realm = m;
    if (m < n)  /* make up some fake rows */
       m = n;

    eqn_names = (char **)malloc(sizeof(char *)*m);
    a = matrix(1,m,1,n);
    save = matrix(1,m,1,n);
    v = matrix(1,n,1,n);
    t = matrix(1,m,1,n);  /* temporary matrix */
    w = NRvector(1,n);
    x = NRvector(1,n);	/* get final solution here */
    b = NRvector(1,m);	/* constant vector of equation */
    f = fopen( argv[1],"r" );

    for(i=1;i<=realm; i++) {  /* read equation i */
        if( ls->Get( f ) <= 0 )
            break;
        printf("line='%s' a[92]=%d\n", ls->line, a[92]);
        if( *ls->line == '#' ) {
        i--;
        continue;
            }
    double term = 1.0;
    eqn_names[i] = NULL;
    for(int j=1; j<=n; j++)  /* zero the row */
        a[i][j] = 0.0;
    char *terminate = " \n\t\0";
    for( char *t = strtok(ls->line, terminate);
        t != NULL;
        t = strtok(NULL, terminate) ) {  /* read each term */

        if( t[0] == '=' ) {
            b[i] = atof(strtok(NULL, terminate));  // then b[i] is the next token
            break;                               // and we are done with this line
        }
            else if( t[0] == 'X' ) {
        int indx = atoi(t+1);
        a[i][indx] = term;
        }
        else if( t[0] == '"' )
        eqn_names[i] = strdup(t);
        else
        term = atof(t);
        }
    }

    fclose( f );
    delete ls;

    /* augment # of rows, if necessary */
    if( realm < m )
        for(i=realm+1; i <= m; i++) {
        for(j=1; j<=n; j++)  /* zero the row */
        a[i][j] = 0.0;
        b[i] = 0.0;
        }

    for(i=1; i<=m; i++) {
       for(j=1; j<=n; j++)
      printf("a[%d][%d] = %f  ", i, j, a[i][j]);
       printf("\n");
       }
    for(i=1; i<=m; i++)
       for(j=1; j<=n; j++)
      save[i][j] = a[i][j]; /* svdcmp will kill array a */
    svdcmp(a,m,n,w,v);
    {
       double wmax = 0.0, wmin;
       for(i=1; i<=n; i++) if (w[i] > wmax) wmax = w[i];
       for(i=1; i<=n; i++) {
       printf("w[%d] = %g\n",i,w[i]);
       for(j=1; j<=n; j++)
           printf("%10g ", v[i][j]);
       printf("\n");
       }
       wmin = wmax * 1.0E-10;
       for(i=1; i<=n; i++)
      if( w[i] < wmin ) {
         printf("Singular value! Index %d\n",i);
             w[i] = 0.0;
         }
    }
    /* for(i=1; i<=m; i++)
    printf("%10g %10g\n", a[i][1], a[i][2]);
    printf("\n");
    for(i=1; i<=n; i++)
    printf("   %10g        %10g %10g\n", w[i], v[i][1], v[i][2] ); */

    for(i=1; i<=m; i++)
    for(j=1; j<=n; j++)
        t[i][j] = a[i][j] * w[j];

    /* for(i=1; i<=m; i++)
    printf("%10g %10g\n", t[i][1], t[i][2]);
    for(i=1; i<=m; i++) {
    for(j=1; j<=n; j++) {
        double s = 0.0;
        int k;
        for(k=1; k<=n; k++)
           s += t[i][k] * v[j][k];
        printf("%12g", s);
        }
    printf("\n");
    } */
    svbksb(a,w,v,m,n,b,x);
    printf("\n");
    printf("Solution: ");
    FILE *fpc = fopen("coeff", "w"); // also write coefficients to a file
    if( fpc == NULL ) {
        printf("Cannot open 'coeff' for write\n");
    exit(-1);
    }
    for(i=1; i<=n; i++) {
        printf("%11g ", x[i]);
        fprintf(fpc,"%11g\n", x[i]);
    }
    printf("\n");
    fclose(fpc);
    {
    double err = 0.0;
    double err2 = 0.0;
    double avg = 0.0;
    double min_err = 1.0E30;
    double max_err = -1.0E30;
    double min_data =  1.0E30;
    double max_data = -1.0E30;
    double var; // variance
    for(i=1; i<=realm; i++) {
    double sum = 0.0;
    double e;
    double e_percent;
    for(j=1; j<=n; j++)
       sum += save[i][j]*x[j];
    min_data = fmin(min_data, b[i]);
    max_data = fmax(max_data, b[i]);
    e = b[i]-sum;
    min_err = fmin(min_err, e);
    max_err = fmax(max_err, e);
    if( argc >= 3 ) {
       if( eqn_names[i] != NULL )
           printf("%s: ", eqn_names[i]);
       else
           printf("Eqn %4d: ", i);
       e_percent = e/b[i]*100.0;
       printf("expect %12g, got %12g, residual %8g, pct %.3f\n", b[i], sum, e, e_percent);
       }
    avg += e;
    err += fabs(e);
    err2 += e*e;
    }
    printf("Data range [%g,%g], error range [%g,%g]\n",
     min_data, max_data, min_err, max_err);
    avg = avg/m;
    var = (err2 - m*avg*avg)/(m-1);
    printf("Average of errors is %g, variance %g\n", avg, var);
    printf(" Average error %g, RMS err %g\n",
     err/m, sqrt(err2/m) );
    }
return 0;
}


