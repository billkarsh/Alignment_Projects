

#include	"CTemplate.h"
#include	"Maths.h"
#include	"Correlation.h"






/* --------------------------------------------------------------- */
/* FillFrWithNormValues ------------------------------------------ */
/* --------------------------------------------------------------- */

void Template::FillFrWithNormValues(
    vector<double>	&fr,
    const PicBase	&P,
    int				i,
    int				xmin,
    int				xmax,
    int				ymin,
    int				ymax )
{
    int	Nx = xmax - xmin + 1,
        Ny = ymax - ymin + 1;

// Get stats

    vector<double>	px( Nx * Ny );
    double			avg, std;

    for( int y = ymin, k = 0; y <= ymax; ++y ) {

        for( int x = xmin; x <= xmax; ++x )
            px[k++] = P.raster[x+P.w*y];
    }

    Stats( px, avg, std );
    printf( "Template: %d, mean %f, std %f\n", i, avg, std );

// Fill with normalized values

    fr.resize( nx * ny, 0.0 );

    for( int y = ymin; y <= ymax; ++y ) {

        for( int x = xmin; x <= xmax; ++x )
            fr[(x-xmin)+nx*(y-ymin)] = (P.raster[x+P.w*y] - avg)/std;
    }
}

/* --------------------------------------------------------------- */
/* Template::Template -------------------------------------------- */
/* --------------------------------------------------------------- */

Template::Template(
    const PicBase	&P,
    int				i,
    int				xmin,
    int				xmax,
    int				ymin,
    int				ymax )
{
    int	NX = xmax - xmin + 1,
        NY = ymax - ymin + 1;

    nx = CeilPow2( 2 * NX );
    ny = CeilPow2( 2 * NY );

    M = ny*(nx/2+1);

    printf( "Template: Transform is size nx=%d, ny=%d\n", nx, ny );

// for the same cost, we can have a larger region.

    int	nmax, bigger;

    nmax	= nx - NX;
    bigger	= (nmax - NX)/2;
    xmin	= max( 0, xmin - bigger );
    xmax	= min( P.w-1, xmax + bigger );

    nmax	= ny - NY;
    bigger	= (nmax - NY)/2;
    ymin	= max( 0, ymin - bigger );
    ymax	= min( P.h-1, ymax + bigger );

    printf(
    "Template: New improved data size x=[%d,%d], y=[%d,%d].\n",
    xmin, xmax, ymin, ymax );

    x0 = xmin;
    y0 = ymin;

// Get frame values

    vector<double>	fr;

    FillFrWithNormValues( fr, P, i, xmin, xmax, ymin, ymax );

// Now fft it

    FFT_2D( fft, fr, nx, ny, false );
}

/* --------------------------------------------------------------- */
/* Template::Match ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Match a chunk of the jth picture to the stored template.
//
Point Template::Match(
    const PicBase	&P,
    int				j,
    int				xmin,
    int				ymin,
    int				xmax,
    int				ymax )
{
    int				NX = xmax - xmin + 1,
                    NY = ymax - ymin + 1;

    printf(
    "Template: Matching x=[%d %d], y=[%d %d] of image %d\n",
    xmin, xmax, ymin, ymax, j );

// Get frame values

    vector<double>	fr;

    FillFrWithNormValues( fr, P, j, xmin, xmax, ymin, ymax );

// FFTs and lags

    vector<CD>	tfft;

    FFT_2D( tfft, fr, nx, ny, false );

    for( int k = 0; k < M; ++k )
        tfft[k] = fft[k] * conj( tfft[k] );

    IFT_2D( fr, tfft, nx, ny );

// Now find the maximum value

    double	norm	= NX * NY * double(nx) * double(ny);
    double	biggest	= -1.0E30;
    int		bigx	= -1;
    int		bigy	= -1;

    for( int iy = 0; iy < ny; ++iy ) {

        for( int ix = 0; ix < nx; ++ix ) {

            int	k = ix + nx*iy;

            if( fr[k] > biggest ) {

                int	x = ix,
                    y = iy;

                if( x > nx/2 ) x -= nx;
                if( y > ny/2 ) y -= ny;

                biggest = fr[k];
                bigx	= x;
                bigy	= y;
            }
        }
    }

    biggest /= norm;

    PrintCorLandscape( biggest, bigx, bigy, 0, 0, 0,
        6, 2, &fr[0], nx, ny, norm, stdout );

    printf( "Template: Match: Maximum %f at (%d, %d).\n",
        biggest, bigx, bigy );

    Point	pt( bigx, bigy );
    ParabPeakFFT( pt.x, pt.y, 1, &fr[0], nx, ny );
    pt.x += x0 - xmin;
    pt.y += y0 - ymin;

    printf( "Template: Match: Final at (%f, %f).\n",
        pt.x, pt.y );

    return pt;
}


