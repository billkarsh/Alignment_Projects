

#include	"GenDefs.h"

#include	<fftw3.h>

#include	<stdlib.h>






void DFT(int N, CD* in, CD *out, CD*table)
{
    for(int k=0; k<N; k++) {

        CD sum(0.0,0.0);

        for(int i=0; i<N; i++) {
            int indx = (i*k) & (N-1);  // and is like mod for 2^N
            sum += in[i]*table[indx];
        }

        out[k] = sum;
    }
}


int main(int argc, char **argv)
{
    const int		N = 1024;
    complex<double>	table[N];

    CD *in = (CD*) fftw_malloc(sizeof(CD) * N);

    for(int i=0; i<N; i++) {
        in[i] = complex<double>(rand()/1000000000.0, rand()/1000000000.0);
        double ang = 2*3.14159265358*double(i)/N;
        table[i] = complex<double>(cos(ang), -sin(ang));
    }

    CD *outFFT = (CD*) fftw_malloc(sizeof(CD) * N);
    CD *outDFT = (CD*) fftw_malloc(sizeof(CD) * N);

    fftw_plan p;
    p = fftw_plan_dft_1d( N, (double (*)[2])in, (double (*)[2])outFFT, FFTW_FORWARD, FFTW_ESTIMATE );

    int	n1 = (argc >= 2) ? atoi(argv[1]) : 1;
    int	n2 = (argc >= 3) ? atoi(argv[2]) : 1;

    for(int i=0; i<n1; i++)
        fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);

    for(int i=0; i<n2; i++)
        DFT(N, in, outDFT, table);

    for(int i=0; i<10; i++)
        printf("%d (%f %f) (%f %f)\n", i, outFFT[i].real(), outFFT[i].imag(),  outDFT[i].real(), outDFT[i].imag() );

    fftw_free(in); fftw_free(outDFT); fftw_free(outFFT);
    return 0;
}
