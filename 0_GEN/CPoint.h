

#pragma once


#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* class Point --------------------------------------------------- */
/* --------------------------------------------------------------- */

class Point {

public:
	double	x, y;

public:
	Point() : x(0.0), y(0.0) {};
	Point( double x, double y ) : x(x), y(y) {};

	double Dist( const Point& rhs ) const;
	double DistSqr( const Point& rhs ) const;

	inline double RSqr() const		{return x*x + y*y;};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

Point FindCOG( const vector<Point> &v, const vector<double> &vals );
Point FindCOG( const vector<Point> &v );

void MakeZeroBasedPoints( vector<Point> &P, int w, int h );

void Mangle( Point &p, int w, int h );

void Set4Corners( vector<Point> &cnr, int w, int h );


