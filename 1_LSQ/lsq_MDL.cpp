

#include	"lsq_MDL.h"
#include	"lsq_Types.h"

#include	"GenDefs.h"
#include	"File.h"

#include	<math.h>
#include	<string.h>


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
/* DisplayStrings ------------------------------------------------ */
/* --------------------------------------------------------------- */

// TrakEM2-style title buffer is optional.
//
// path is not a buffer, it's a char pointer!!
//
void MDL::DisplayStrings(
	char			*title,
	const char*		&path,
	const RGN		&I )
{
	const Til2Img	*m;

	RGN::GetMeta( &m, NULL, I, I );
	path = m->path.c_str();

	if( title ) {

		if( m->col != -999 ) {

			sprintf( title, "%d.%d:1_%d.%d.%d",
				I.z, m->id, m->col, m->row, m->cam );
		}
		else
			sprintf( title, "%d.%d:1", I.z, m->id );
	}
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

	vector<Point>	cnr;
	double			xmin, xmax, ymin, ymax;
	int				nr = vRgn.size();

	Set4Corners( cnr, gW, gH );

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

			vector<Point>	c( 4 );
			memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
			L2GPoint( c, X, itr );

			for( int k = 0; k < 4; ++k ) {
				xmin = fmin( xmin, c[k].x );
				xmax = fmax( xmax, c[k].x );
				ymin = fmin( ymin, c[k].y );
				ymax = fmax( ymax, c[k].y );
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

		vector<Point>	c( 4 );
		double			xmid = 0.0, ymid = 0.0;

		memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
		L2GPoint( c, X, itr );

		for( int k = 0; k < 4; ++k ) {
			xmid += c[k].x;
			ymid += c[k].y;
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

	if( FOUT )
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


