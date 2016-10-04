

// Notes
// -----
// MKL is preferred for FFT operations over fftw if available:
//
// (1) MKL executes transforms in ~ 60% to 75% of the time.
// (2) MKL is thread-safe allowing better MT efficiency.
//

#include	"mkl.h"


/* --------------------------------------------------------------- */
/* MKLCheck ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void MKLCheck( MKL_LONG status, FILE *flog )
{
    if( status && !DftiErrorClass( status, DFTI_NO_ERROR ) ) {

        fprintf( flog,
        "MKL error [%s].\n", DftiErrorMessage( status ) );
        exit( 44 );
    }
}

/* --------------------------------------------------------------- */
/* _FFT_2D ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void _FFT_2D(
    vector<CD>				&out,
    const vector<double>	&in,
    int						Nfast,
    int						Nslow,
    FILE					*flog )
{
    int	Nhlf = (Nfast/2 + 1), M = Nslow * Nhlf;

    out.resize( M );

    DFTI_DESCRIPTOR_HANDLE	h;
    MKL_LONG				dim[2]  = {Nslow, Nfast},
                            stro[3] = {0, Nhlf, 1},
                            status;

    status = DftiCreateDescriptor( &h,
                DFTI_DOUBLE,
                DFTI_REAL,
                2, dim );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_NUMBER_OF_USER_THREADS,
                1 );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_CONJUGATE_EVEN_STORAGE,
                DFTI_COMPLEX_COMPLEX );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_OUTPUT_STRIDES,
                stro );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_PLACEMENT,
                DFTI_NOT_INPLACE );
    MKLCheck( status, flog );

    status = DftiCommitDescriptor( h );
    MKLCheck( status, flog );

    status = DftiComputeForward( h,
                (double*)&in[0],
                &out[0] );
    MKLCheck( status, flog );

    status = DftiFreeDescriptor( &h );
    MKLCheck( status, flog );
}

/* --------------------------------------------------------------- */
/* FFT_2D -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Forward FFT of 2D data (real to complex).
//
// Assumes input data in row-major order. That is,
// ordered like a C-array: in[Nslow][Nfast].
//
int FFT_2D(
    vector<CD>				&out,
    const vector<double>	&in,
    int						Nfast,
    int						Nslow,
    bool					cached,
    FILE					*flog )
{
    int	M = Nslow * (Nfast/2 + 1);

    if( cached ) {

        pthread_mutex_lock( &mutex_fft );

        if( out.size() != M )
            _FFT_2D( out, in, Nfast, Nslow, flog );

        pthread_mutex_unlock( &mutex_fft );
    }
    else
        _FFT_2D( out, in, Nfast, Nslow, flog );

    return M;
}

/* --------------------------------------------------------------- */
/* IFT_2D -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Inverse FFT of 2D data (complex to real).
//
// Creates output data in row-major order. That is,
// ordered like a C-array: out[Nslow][Nfast].
//
void IFT_2D(
    vector<double>			&out,
    const vector<CD>		&in,
    int						Nfast,
    int						Nslow,
    FILE					*flog )
{
    int	Nhlf = (Nfast/2 + 1);

    out.resize( Nslow * Nfast );

    DFTI_DESCRIPTOR_HANDLE	h;
    MKL_LONG				dim[2]  = {Nslow, Nfast},
                            stri[3] = {0, Nhlf,  1},
                            stro[3] = {0, Nfast, 1},
                            status;

    status = DftiCreateDescriptor( &h,
                DFTI_DOUBLE,
                DFTI_REAL,
                2, dim );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_NUMBER_OF_USER_THREADS,
                1 );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_CONJUGATE_EVEN_STORAGE,
                DFTI_COMPLEX_COMPLEX );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_INPUT_STRIDES,
                stri );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_OUTPUT_STRIDES,
                stro );
    MKLCheck( status, flog );

    status = DftiSetValue( h,
                DFTI_PLACEMENT,
                DFTI_NOT_INPLACE );
    MKLCheck( status, flog );

    status = DftiCommitDescriptor( h );
    MKLCheck( status, flog );

    status = DftiComputeBackward( h,
                (double*)&in[0],
                &out[0] );
    MKLCheck( status, flog );

    status = DftiFreeDescriptor( &h );
    MKLCheck( status, flog );
}


