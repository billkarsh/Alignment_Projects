

#include	"fftw3.h"


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* _FFT_2D ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void _FFT_2D(
	vector<CD>				&out,
	const vector<double>	&in,
	int						Nfast,
	int						Nslow )
{
	int	M = Nslow * (Nfast/2 + 1);

	out.resize( M );

	fftw_plan	p;

	p = fftw_plan_dft_r2c_2d( Nslow, Nfast, (double*)&in[0],
			(double (*)[2])&out[0], FFTW_ESTIMATE );

	fftw_execute( p );
	fftw_destroy_plan( p );
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
	bool					cached )
{
	int	M = Nslow * (Nfast/2 + 1);

	pthread_mutex_lock( &mutex_fft );

	if( !cached || out.size() != M )
		_FFT_2D( out, in, Nfast, Nslow );

	pthread_mutex_unlock( &mutex_fft );

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
	int						Nslow )
{
	int	N = Nslow * Nfast;

	out.resize( N );

	pthread_mutex_lock( &mutex_fft );

	fftw_plan	p;

	p = fftw_plan_dft_c2r_2d( Nslow, Nfast, (double (*)[2])&in[0],
			&out[0], FFTW_ESTIMATE );

	fftw_execute( p );
	fftw_destroy_plan( p );

	pthread_mutex_unlock( &mutex_fft );
}


