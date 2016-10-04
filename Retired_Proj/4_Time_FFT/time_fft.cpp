

#include	"GenDefs.h"

#include	<fftw3.h>

#include	<stdlib.h>






int main(int argc, char **argv)
{
int N = (argc >= 2) ? atoi(argv[1]) : 64;     // size of image
int n1 = (argc >= 3) ? atoi(argv[2]) : 1000;  // times to try


// input is NxN real numbers
int M = N*(N/2+1);  // output is this many complex numbers
double *in1 = (double *)fftw_malloc(sizeof(double) * N*N);
double *in2 = (double *)fftw_malloc(sizeof(double) * N*N);
double *out = (double *)fftw_malloc(sizeof(double) * N*N);
CD *out1 = (CD*) fftw_malloc(sizeof(CD) * M);
CD *out2 = (CD*) fftw_malloc(sizeof(CD) * M);

// make up some data
for(int i=0; i<N*N; i++) {
    in1[i] = rand()/1000000000.0;
    in2[(i+2*N+3)%(N*N)] = in1[i];    // offsets should be 3 and 2 pixels
    //in2[i] = in1[i];
    }

fftw_plan p1, p2, p3;
p1 = fftw_plan_dft_r2c_2d( N, N, in1, (double (*)[2])out1, FFTW_ESTIMATE );
p2 = fftw_plan_dft_r2c_2d( N, N, in2, (double (*)[2])out2, FFTW_ESTIMATE );
p3 = fftw_plan_dft_c2r_2d( N, N, (double (*)[2])out1, out, FFTW_ESTIMATE );

double big;
int bigj;
printf("Starting %d passes\n", n1);
for(int i=0; i<n1; i++) {
    fftw_execute(p1); /* repeat as needed */
    fftw_execute(p2);
    for(int j=0; j<M; j++)
    out1[j] = out1[j]*conj(out2[j]);
    fftw_execute(p3);
    big=-1.0E30;
    bigj = -1;
    for(int j=0; j<N*N; j++) {
    if( out[j] > big ) {
       big = out[j];
           bigj = j;
           }
        }
    }
int x = bigj/N;
int y = bigj - N*x;
if (x > N/2) x -= N;
if (y > N/2) y -= N;
printf("x,y %d %d\n", x, y);
fftw_destroy_plan(p1);
fftw_destroy_plan(p2);
fftw_destroy_plan(p3);
return 0;
}
