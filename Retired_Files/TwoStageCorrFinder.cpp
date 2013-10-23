/* --------------------------------------------------------------- */
/* FFTPeaks ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static int FFTPeaks(
	vector<double>			&topnr,
	vector<Point>			&topnp,
	int						w1,
	int						h1,
	int						w2,
	int						h2,
	vector<double>			&rslt,
	int						N,
	int						Ox,
	int						Oy,
	int						radius,
	FILE*					flog )
{
	const double	pkrad	= 5.0;
	int				tnmax	= 100;

	int		i, j, Nsqr = N * N;

	topnr.resize( tnmax, -1E30 );
	topnp.resize( tnmax, Point(-1E30, -1E30) );

// Max lag

	for( i = 0; i < Nsqr; ++i ) {

		// calculate signed test coordinates

		int	y = i / N;
		int	x = i - N * y;

		if( y > N - h1 ) {				// negative zone

			y -= N;

			if( Oy - y > radius )
				goto skip;
		}
		else if( y >= h2 )				// dead zone
			goto skip;
		else if( y - Oy > radius )		// positive zone
			goto skip;

		if( x > N - w1 ) {				// negative zone

			x -= N;

			if( Ox - x > radius )
				goto skip;
		}
		else if( x >= w2 )				// dead zone
			goto skip;
		else if( x - Ox > radius ) {	// positive zone
skip:
			rslt[i] = 0.0;
			continue;
		}

		// associate value with listed peak?

		for( j = 0; j < tnmax; ++j ) {

			Point	p( x, y );

			if( abs( p.x - topnp[j].x ) <= pkrad &&
				abs( p.y - topnp[j].y ) <= pkrad ) {

				// belongs to this peak

				if( rslt[i] >= topnr[j] ) {

					// improves this peak

					topnr[j] = rslt[i];
					topnp[j] = p;
				}

				break;
			}
			else if( rslt[i] >= topnr[j] ) {

				// new best peak - add to list

				topnr.resize( tnmax - 1 );
				topnp.resize( tnmax - 1 );

				topnr.insert( topnr.begin() + j, rslt[i] );
				topnp.insert( topnp.begin() + j, p );

				break;
			}
		}
	}

// Consolidate peaks

	for( j = tnmax - 1; j >= 0; --j ) {

		// remove unused slots

		if( topnr[j] == -1E30 || topnp[j].x == -1E30 ) {

			topnr.erase( topnr.begin() + j );
			topnp.erase( topnp.begin() + j );
			--tnmax;
			continue;
		}

		// join close peaks with their betters at lower index

		for( i = j - 1; i >= 0; --i ) {

			if( abs( topnp[j].x - topnp[i].x ) <= pkrad &&
				abs( topnp[j].y - topnp[i].y ) <= pkrad ) {

				topnr.erase( topnr.begin() + j );
				topnp.erase( topnp.begin() + j );
				--tnmax;
				break;
			}
		}
	}

// Report best peaks

	fprintf( flog, "\nTop N best cross-prod: %d\n", tnmax );
	for( j = 0; j < tnmax; ++j ) {
		fprintf( flog, "[%.3e] ", topnr[j] );
		if( (j + 1) % 10 == 0 ) fprintf( flog, "\n" );
	}
	fprintf( flog, "\n" );

	fprintf( flog, "\nTop N best cross-prod: %d\n", tnmax );
	for( j = 0; j < tnmax; ++j ) {
		fprintf( flog, "[%d,%d] ", (int)topnp[j].x, (int)topnp[j].y );
		if( (j + 1) % 10 == 0 ) fprintf( flog, "\n" );
	}
	fprintf( flog, "\n\n" );

	return tnmax;
}

/* --------------------------------------------------------------- */
/* CLinCorr ------------------------------------------------------ */
/* --------------------------------------------------------------- */

class CLinCorr {

private:
	vector<double>	i1sum, i1sum2,
					i2sum, i2sum2;
	vector<int>		i1nz,  i2nz;
	int				w1,  h1,
					w2,  h2,
					N,   Nsq;
	IBox			OL1, OL2;
	double*			rslt;
	int				ir,
					dx,  dy,
					olw, olh,
					i1c, i2c;

public:
	void Initialize(
		const vector<double>	&I1,
		int						w1,
		int						h1,
		const vector<double>	&I2,
		int						w2,
		int						h2,
		int						N );

	int CheckSize(
		vector<double>	&rslt,
		int				irslt,
		EvalType		LegalRgn,
		void*			arglr,
		int				dx,
		int				dy );

	int CheckDensity(
		EvalType		LegalCnt,
		void*			arglc );

	double LinearCorr();
	int    SizeIndex();
};


void CLinCorr::Initialize(
		const vector<double>	&I1,
		int						w1,
		int						h1,
		const vector<double>	&I2,
		int						w2,
		int						h2,
		int						N )
{
	this->w1	= w1;
	this->h1	= h1;
	this->w2	= w2;
	this->h2	= h2;
	this->N		= N;
	Nsq			= N * N;

	IntegrateImage( i1sum, i1sum2, i1nz, w1, h1, I1, N );
	IntegrateImage( i2sum, i2sum2, i2nz, w2, h2, I2, N );
}


int CLinCorr::CheckSize(
	vector<double>	&rslt,
	int				irslt,
	EvalType		LegalRgn,
	void*			arglr,
	int				dx,
	int				dy )
{
	int		ok = true;

	this->rslt	= &rslt[0];
	ir			= irslt;
	this->dx	= dx;
	this->dy	= dy;

	BoxesFromShifts( OL1, OL2, w1, h1, w2, h2, dx, dy );

// Large enough overlap?

	olw = OL1.R - OL1.L + 1;
	olh = OL1.T - OL1.B + 1;

	if( LegalRgn && !LegalRgn( olw, olh, arglr ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


int CLinCorr::CheckDensity(
	EvalType		LegalCnt,
	void*			arglc )
{
	int		ok = true;

	i1c = IntegralTable( i1nz, w1, OL1 );
	i2c = IntegralTable( i2nz, w2, OL2 );

	if( LegalCnt && !LegalCnt( i1c, i2c, arglc ) ) {

		rslt[ir]	= 0.0;
		ok			= false;
	}

	return ok;
}


double CLinCorr::LinearCorr()
{
	double	n = olw * olh;

	double	im1sum	= IntegralTable( i1sum,  w1, OL1 );
	double	im1sum2	= IntegralTable( i1sum2, w1, OL1 );
	double	im2sum	= IntegralTable( i2sum,  w2, OL2 );
	double	im2sum2	= IntegralTable( i2sum2, w2, OL2 );

	double	num	= n * rslt[ir] / Nsq - im2sum * im1sum;
	double	d1	= n * im1sum2 - im1sum * im1sum;
	double	d2	= n * im2sum2 - im2sum * im2sum;
	double	d	= d1 * d2;
	double	r	= (d < n * n * 1.0E-9 ? 0.0 : num / sqrt( d ));

	return rslt[ir] = r;
}


int CLinCorr::SizeIndex()
{
	return (int)(10.0 * log( max( 1, min(i1c,i2c) ) ));
}

/* --------------------------------------------------------------- */
/* CorrPatches --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return the normalized cross-correlation and (dx, dy), which
// must be added to points in patch1 to match patch2.
//
// Search confined to disc: (origin, radius) = {Ox, Oy, radius}.
//
// 'fft2' is a cache of the patch2 FFT. On entry, if fft2 has
// the correct size it is used. Otherwise recomputed here.
//
double CorrPatches(
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	double					&dx,
	double					&dy,
	int						Ox,
	int						Oy,
	int						radius,
	FILE*					flog,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 )
{
	of = flog;

// Bounding boxes of point lists

	IBox	B1, B2;
	int		w1, h1, w2, h2;

	BBoxFromPoints( B1, ip1 );
	BBoxFromPoints( B2, ip2 );

	w1 = B1.R - B1.L + 1;
	h1 = B1.T - B1.B + 1;

	w2 = B2.R - B2.L + 1;
	h2 = B2.T - B2.B + 1;

	fprintf( of,
	"NormCorr: region size is [%d %d] in x, [%d %d] in y.\n",
	B1.L, B1.R, B1.B, B1.T );

	fprintf( of,
	"NormCorr: target size is [%d %d] in x, [%d %d] in y.\n",
	B2.L, B2.R, B2.B, B2.T );

// Find N, the size we need for FFTs.
// We will always use square FFTs.

	int		N = FFTSizeSQR( w1, h1, w2, h2 );
	fprintf( of, "NormCorr: N = %d\n", N );

// Number of complex values in FFT of 2D real

	int		M = N*(N/2+1);

// Create images from point lists.

	vector<double>	i1, i2;

	ImageFromValuesAndPoints( i1, N, N, iv1, ip1, B1.L, B1.B );
	ImageFromValuesAndPoints( i2, N, N, iv2, ip2, B2.L, B2.B );

// FFTs and lags

	vector<double>	rslt;
	vector<CD>		fft1;

	if( fft2.size() != M )
		FFT_2D( fft2, i2, N, N );

	FFT_2D( fft1, i1, N, N );

	for( int i = 0; i < M; ++i )
		fft1[i] = fft2[i] * conj( fft1[i] );

	IFT_2D( rslt, fft1, N, N );

// Cross-correlation peak candidates

	vector<double>	topnr;
	vector<Point>	topnp;
	int				tnmax;

	tnmax = FFTPeaks( topnr, topnp, w1, h1, w2, h2,
				rslt, N, Ox, Oy, radius, flog );

// Create array indexed by 'size': int( 10 * log( overlap_size ) ).
// Each element contains the best rslt[i] at that size index.
// The idea is to look for correlation peaks that are nearly
// as good as the best, but that derive from larger overlaps.

	int lmax = (int)(10.0*log(double(N*N))) + 1;
	vector<double>	max_by_size( lmax, 0.0 );

// Refine peaks
//
// For each of the candidate peaks, examine a small neighborhood
// (2D+1)^2 of pixels about the peak and seek the best value of
// the linear correlation coefficient r = sxy/sqrt(sxx * syy).
//
// At the same time we store the r values back into the rslt
// image for sake of the PrintCorLandscape diagnostic.

	CLinCorr	LC;

	LC.Initialize( i1, w1, h1, i2, w2, h2, N );

	double	biggest = -1.0E30;
	int		bigx	= -1,
			bigy	= -1;
	int		D		= 3;

	for( int it = 0; it < tnmax; ++it ) {

		// coordinate names:
		// s -> signed
		// i -> unsigned indices into rslt[]

		// y-neighbors

		for( int ky = -D; ky <= D; ++ky ) {

			int	iy, sy;

			iy = sy = (int)topnp[it].y + ky;

			if( sy < 0 )
				iy += N;

			// x-neighbors

			for( int kx = -D; kx <= D; ++kx ) {

				int	i, ix, sx;

				ix = sx = (int)topnp[it].x + kx;

				if( sx < 0 )
					ix += N;

				i = ix + iy * N;

				if( sy <= -h1 || sy >= h2 || sx <= -w1 || sx >= w2 ) {

					if( i >= 0 && i < N*N )
						rslt[i] = 0.0;

					continue;
				}

				// Refine pixel {sx, sy}

				if( !LC.CheckSize( rslt, i, LegalRgn, arglr, sx, sy ) )
					continue;

				if( !LC.CheckDensity( LegalCnt, arglc ) )
					continue;

				double	r = LC.LinearCorr();

				if( r > biggest ) {
					biggest	= r;
					bigx	= sx;
					bigy	= sy;
				}

				// Update max_by_size

				int im = LC.SizeIndex();

				if( r > max_by_size[im] )
					max_by_size[im] = r;
			}
		}
	}

// Reports

	if( biggest < -2.0 ) {

		fprintf( of, "NormCorr: No legal subregions at all...\n");

		// for sake of grep searches

		fprintf( of,
		"NormCorr: Maximum correlation of %f at [%d,%d].\n",
		0.0, 0, 0 );

		return 0.0;
	}

	PrintCorLandscape( biggest, bigx, bigy, Ox, Oy, radius,
		D, 1, &rslt[0], N, N, 1.0, of );

	//fprintf( of, "----v-center-v---\n" );
	//PrintCorLandscape( biggest, 32, 32, Ox, Oy, radius,
	//	D, 1, &rslt[0], N, N, 1.0, of );

// For debugging, print the results of correlation by region size.
// If we find a larger region with a correlation almost as good,
// that's potentially a better match.

	double	bc		= -10.0;		// biggest correlation
	int		we		= 0;			// which entry had it?
	int		nmax	= max_by_size.size();

	for( int i = 0; i < nmax; ++i ) {

		if( max_by_size[i] > bc ) {

			bc = max_by_size[i];
			we = i;
		}
	}

// Now print the entries with bigger size and comparable correlation
// Mpy by 64 to get into pixel counts in the 2K working image.

	double PT = 0.8 * bc;	// Print threshold

	for( int i = we + 1; i < nmax; ++i ) {

		if( max_by_size[i] >= PT ) {

			int i1 = (int)ceil(  64.0 * exp( i/10.0 ) );
			int i2 = (int)floor( 64.0 * exp( (i+1)/10.0 ) );

			fprintf( of,
			"NormCorr: Possible bigger area match:"
			" %8d - %8d : %8.2f\n", i1, i2, max_by_size[i] );
		}
	}

	fprintf( of,
	"NormCorr: Maximum correlation of %f at [%d,%d].\n",
	biggest, bigx, bigy );

	//double limit = 0.9*biggest;
	//for(int x=-w1+1; x<w2-1; x++) {
	//for(int y=-h1+1; y<h2-1; y++) {
	//double a = LkFFT(rslt, N, N, x, y);
	//if (a > limit && a > LkFFT(rslt, N, N, x-1, y) && a > LkFFT(rslt, N, N, x+1, y) &&
	//a > LkFFT(rslt, N, N, x, y-1) && a > LkFFT(rslt, N, N, x, y+1) ){
	//fprintf(of, "Local max at %d %d, value %f\n", x, y, a);
	//}
	//}
	//}

// Interpolate peak

	dx = bigx;
	dy = bigy;
	ParabPeakFFT( dx, dy, 1, &rslt[0], N, N );
	dx += B2.L - B1.L;
	dy += B2.B - B1.B;

	return biggest;
}


