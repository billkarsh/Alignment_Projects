

#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#include	<stdlib.h>


double **
matrix(nrl,nrh,ncl,nch)
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

double *
vector(nl, nh)
{
    double *v;
    v = (double *)malloc((nh-nl+1)*sizeof(double));
    return v-nl;
}
nrerror(c)
 char *c;
{
    printf("%s\n", c);
    exit(1);
}

free_vector(v, nl, nh)
 double *v;
{
    free((char *)(v+nl));
}

double fabs(double f);
double fmin(double a, double b) { return a < b ? a : b;}
double fmax(double a, double b) { return a > b ? a : b;}

int main(argc, argv)
 int argc;
 char **argv;
{
    double **a, *w, **v, **t, **save;
    double *x, *b;
    int i,j,  m,n;
    int realm;	/* real # of equations */
    FILE *f;
    char **eqn_names;

    m = n = 0;
    if( argc < 3 ) {
	printf("Usage: svdfit <file of equations> -print\n");
	exit(-1);
	}
    f = fopen(argv[1],"r");
    if( f == NULL ) {
	printf("Can't open file '%s'\n", argv[1]);
	exit(-1);
	}
    for(;;) {
        char c;
        if (fscanf(f,"%c", &c) != 1) break;
        if( c == '=' )
	    m++;
	if( c == 'X' ) {
	    int indx;
	    fscanf(f, "%d", &indx);
	    if( indx > n )
		n = indx;
	    }
	}
    fclose(f);
    printf("there are %d equations in %d unknowns\n", m, n);
    realm = m;
    if (m < n)  /* make up some fake rows */
       m = n;

    eqn_names = (char **)malloc(sizeof(char *)*m);
    a = matrix(1,m,1,n);
    save = matrix(1,m,1,n);
    v = matrix(1,n,1,n);
    t = matrix(1,m,1,n);  /* temporary matrix */
    w = vector(1,n);
    x = vector(1,n);	/* get final solution here */
    b = vector(1,m);	/* constant vector of equation */
    f = fopen(argv[1],"r");
    for(i=1;i<=realm; i++) {  /* read equation i */
	int j;
	double term;
	eqn_names[i] = NULL;
	for(j=1; j<=n; j++)  /* zero the row */
	    a[i][j] = 0.0;
	for(term = 1.0;;) {  /* read each term */
	    char token[256];
            if (fscanf(f,"%s", token) != 1) break;
            if( token[0] == '=' ) {
		fscanf(f, "%lf", &b[i]);
		break;
		}
            else if( token[0] == 'X' ) {
		int indx;
		sscanf(token+1, "%d", &indx);
		a[i][indx] = term;
		}
	    else if( token[0] == '"' )
		eqn_names[i] = strdup(token);
	    else
		sscanf(token, "%lf", &term);
	    }
	}
    fclose(f);

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


