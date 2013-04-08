

#include	"lsq_MDL.h"
#include	"lsq_Types.h"

#include	"GenDefs.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* Magnitude ----------------------------------------------------- */
/* --------------------------------------------------------------- */

double MDL::Magnitude( const vector<double> &X, int itr )
{
	Point	p0( 0, 0 ), p1( 1, 1 );

	L2GPoint( p0, X, itr );
	L2GPoint( p1, X, itr );

	return sqrt( p1.DistSqr( p0 ) / 2 );
}

/* --------------------------------------------------------------- */
/* PrintMagnitude ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MDL::PrintMagnitude( const vector<double> &X )
{
	int	ntr = X.size() / NX;

	if( ntr )
		printf( "Final magnitude is %g\n", Magnitude( X, ntr-1 ) );

	fflush( stdout );
}

/* --------------------------------------------------------------- */
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void MDL::Bounds(
	double					&xbnd,
	double					&ybnd,
	vector<double>			&X,
	const vector<double>	&lrbt,
	double					degcw,
	FILE					*FOUT )
{
	printf( "---- Global bounds ----\n" );

// Transform each included regions's rectangle to global
// space, including any global rotation (degcw) and find
// bounds over whole set.

	double	xmin, xmax, ymin, ymax;
	int		nr = vRgn.size();

	if( lrbt.size() ) {

		xmin = lrbt[0];
		xmax = lrbt[1];
		ymin = lrbt[2];
		ymax = lrbt[3];
	}
	else {

		xmin =  BIGD;
		xmax = -BIGD;
		ymin =  BIGD;
		ymax = -BIGD;

		if( degcw )
			RotateAll( X, degcw );

		for( int i = 0; i < nr; ++i ) {

			int	itr = vRgn[i].itr;

			if( itr < 0 )
				continue;

			vector<Point>	cnr( 4 );

			cnr[0] = Point(  0.0, 0.0 );
			cnr[1] = Point( gW-1, 0.0 );
			cnr[2] = Point( gW-1, gH-1 );
			cnr[3] = Point(  0.0, gH-1 );

			L2GPoint( cnr, X, itr );

			for( int k = 0; k < 4; ++k ) {

				xmin = fmin( xmin, cnr[k].x );
				xmax = fmax( xmax, cnr[k].x );
				ymin = fmin( ymin, cnr[k].y );
				ymax = fmax( ymax, cnr[k].y );
			}
		}
	}

	printf( "Propagate bounds with option -lrbt=%f,%f,%f,%f\n\n",
	xmin, xmax, ymin, ymax );

// Translate all transforms to put global origin at ~(0,0).

	NewOriginAll( X, xmin, ymin );

	xmax = ceil( xmax - xmin + 1 );
	ymax = ceil( ymax - ymin + 1 );
	xmin = 0;
	ymin = 0;

// Open GNUPLOT files for debugging

	FILE	*fEven		= FileOpenOrDie( "pf.even", "w" ),
			*fOdd		= FileOpenOrDie( "pf.odd", "w" ),
			*fLabEven	= FileOpenOrDie( "pf.labels.even", "w" ),
			*fLabOdd	= FileOpenOrDie( "pf.labels.odd", "w" );

// Write rects and labels

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		vector<Point>	cnr( 4 );
		double			xmid = 0.0, ymid = 0.0;

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		L2GPoint( cnr, X, itr );

		for( int k = 0; k < 4; ++k ) {
			xmid += cnr[k].x;
			ymid += cnr[k].y;
		}

		xmid /= 4.0;
		ymid /= 4.0;

		// select even or odd reporting

		FILE	*f, *flab;
		int		color;

		if( vRgn[i].z & 1 ) {
			f		= fOdd;
			flab	= fLabOdd;
			color	= 1;
		}
		else {
			f		= fEven;
			flab	= fLabEven;
			color	= 2;
		}

		// transformed rect

		for( int k = 0; k < 5; ++k )
			fprintf( f, "%f %f\n", cnr[k%4].x, cnr[k%4].y );

		fprintf( f, "\n" );

		// label

		fprintf( flab, "set label \"%d:%d.%d \" at %f,%f tc lt %d\n",
		vRgn[i].z, vRgn[i].id, vRgn[i].rgn, xmid, ymid, color );
	}

// Close files

	fclose( fLabOdd );
	fclose( fLabEven );
	fclose( fOdd );
	fclose( fEven );

// Report

	fprintf( FOUT, "BBOX %f %f %f %f\n", xmin, ymin, xmax, ymax );

	printf( "Bounds of global image are x=[%f %f] y=[%f %f].\n\n",
	xmin, xmax, ymin, ymax );

	xbnd = xmax;
	ybnd = ymax;
}

/* --------------------------------------------------------------- */
/* EqvAffine ----------------------------------------------------- */
/* --------------------------------------------------------------- */

TAffine MDL::EqvAffine(
	const vector<double>	&X,
	int						itr )
{
	Point	v0( 0, 0 ), v1( 1, 0 ), v2( 0, 1 );

	L2GPoint( v0, X, itr );
	L2GPoint( v1, X, itr );
	L2GPoint( v2, X, itr );

	return TAffine(	v1.x - v0.x, v2.x - v0.x, v0.x,
				    v1.y - v0.y, v2.y - v0.y, v0.y );
}


