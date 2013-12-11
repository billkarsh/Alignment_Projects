

#include	"CPoint.h"

#include	<math.h>
#include	<stdio.h>






/* --------------------------------------------------------------- */
/* Dist ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

double Point::Dist( const Point& rhs ) const
{
	double	dx = x - rhs.x,
			dy = y - rhs.y;

	return sqrt( dx*dx + dy*dy );
}

/* --------------------------------------------------------------- */
/* DistSqr ------------------------------------------------------- */
/* --------------------------------------------------------------- */

double Point::DistSqr( const Point& rhs ) const
{
	double	dx = x - rhs.x,
			dy = y - rhs.y;

	return dx*dx + dy*dy;
}

/* --------------------------------------------------------------- */
/* FindCOG ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Center of gravity of a point vector.
// Use only the non-zero entries, so we only get the values
// that exist in the target image.
//
Point FindCOG( const vector<Point> &v, const vector<double> &vals )
{
	double	sumx = 0.0, sumy = 0.0;
	int		n = v.size(), m = 0;

	for( int i = 0; i < n; ++i ) {

		if( fabs( vals[i] ) > 1.0E-8 ) {
			sumx += v[i].x;
			sumy += v[i].y;
			++m;
		}
	}

	return Point( sumx/m, sumy/m );
}


// Center of gravity of a point vector.
//
Point FindCOG( const vector<Point> &v )
{
	double	sumx = 0.0, sumy = 0.0;
	int		n = v.size();

	for( int i = 0; i < n; ++i ) {
		sumx += v[i].x;
		sumy += v[i].y;
	}

	return Point( sumx/n, sumy/n );
}

/* --------------------------------------------------------------- */
/* MakeZeroBasedPoints ------------------------------------------- */
/* --------------------------------------------------------------- */

void MakeZeroBasedPoints( vector<Point> &P, int w, int h )
{
	int		np = w * h;

	P.resize( np );

	for( int i = 0; i < np; ++i ) {

		int	y = i / w,
			x = i - w * y;

		P[i] = Point( x, y );
	}
}

/* --------------------------------------------------------------- */
/* Mangle -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experimental code to mangle a point with a non-affine transform.
//
void Mangle( Point &p, int w, int h )
{
	double	alpha	= p.x/(w-1);
	double	beta	= p.y/(w-1);
	double	scale	= 1.01;
	Point	p0		= p;
	Point	p1, p2;	// intermediate points

// first vector goes from (0,0) to (w-1,0)
	p1.x = (1-alpha)*0.0 + alpha*(w-1);	// could be simplified..
	p1.y = (1-alpha)*0.0 + alpha*0.0;

// top vector goes from (0,h-1) to (w-1,h-1)*scale
	p2.x = (1-alpha)*0.0   + alpha*((w-1)*scale);
	p2.y = (1-alpha)*(h-1) + alpha*((h-1)*scale);

// now interpolate among these using beta
	p.x = (1-beta)*p1.x + beta*p2.x;
	p.y = (1-beta)*p1.y + beta*p2.y;

	printf( "Mangle: Was (%f %f), now (%f %f).\n",
		p0.x, p0.y, p.x, p.y );
}

/* --------------------------------------------------------------- */
/* Set4Corners --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Convenience to initialize a vector with the four corners
// of an image (or other rectangle).
//
void Set4Corners( vector<Point> &cnr, int w, int h )
{
	cnr.resize( 4 );
	cnr[0] = Point( 0.0, 0.0 );
	cnr[1] = Point( w-1, 0.0 );
	cnr[2] = Point( w-1, h-1 );
	cnr[3] = Point( 0.0, h-1 );
}


