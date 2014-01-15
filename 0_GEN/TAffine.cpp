

#include	"GenDefs.h"
#include	"TAffine.h"

#include	<math.h>






/* --------------------------------------------------------------- */
/* Discussion ---------------------------------------------------- */
/* --------------------------------------------------------------- */

/*
	Coordinates
	-----------
	     ^+Z
	     |
	     |
	     /--->+X   +theta rotates +X toward +Y (clockwise)
	  +Y/

	Individual images lie in the XY-plane, but unlike usual math
	conventions the Y-axis is negated, so increases 'downward'.

	Z-coordinates usually start at zero and increase upward. The
	standard perspective is looking down on a stack.

	Because the Y-axis is negated, positive angles cause a CW
	rather than a CCW rotation.

	Some portions of this code body refer to boxes or bounds or
	bounding boxes wherein the sides are {L,R,B,T}. X increases
	from left to right in the usual way. Y increases from bottom
	to top, but again, the bottom of a box appears higher on a
	computer screen because Y runs downward.

	2D transformations map a region (such as a whole image) from
	local region coords to an abstract global coordinate space.

	In our code, the 6 affine transform elements are interpreted
	as follows:

		x'   |t0  t1|   x     t2
		   = |      | *    +
		y'   |t3  t4|   y     t5

	In Matlab, our 6 elements are filled by column into one array
	as follows:

		{t0, t3, t1, t4, t2, t5}.

	In ImageJ (including TrakEM2) the elements in source code
	are the following:

		x'   |m00  m01|   x     m02
		   = |        | *    +
		y'   |m10  m11|   y     m12

	In a TrackEM2 XML file transform attributes are written
	"matrix(m00,m10,m01,m11,m02,m12)".

	We use a standard order for compounding operations. Let
	X=xskew, Y=yskew, S=scaling, P=pretweak, R=rot+trans.
	Then, D = S(YX), and T = R(DP).
*/

/* --------------------------------------------------------------- */
/* Near unity transforms ----------------------------------------- */
/* --------------------------------------------------------------- */

// pos a enlarges
void TAffine::NUSetScl( double a )
{
	t[0] = a; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = a; t[5] = 0;
}


// pos a enlarges
void TAffine::NUSetXScl( double a )
{
	t[0] = a; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = 1; t[5] = 0;
}


// pos a enlarges
void TAffine::NUSetYScl( double a )
{
	t[0] = 1; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = a; t[5] = 0;
}


// pos a tips right
void TAffine::NUSetXSkw( double a )
{
	t[0] = 1; t[1] = a; t[2] = 0;
	t[3] = 0; t[4] = 1; t[5] = 0;
}


// pos a tips up
void TAffine::NUSetYSkw( double a )
{
	t[0] = 1; t[1] = 0; t[2] = 0;
	t[3] = a; t[4] = 1; t[5] = 0;
}


// pos r rotates CW
void TAffine::NUSetRot( double r )
{
	double	c = cos( r ), s = sin( r );

	t[0] = c; t[1] = -s; t[2] = 0;
	t[3] = s; t[4] =  c; t[5] = 0;
}


void TAffine::NUSelect( int sel, double a )
{
	switch( sel ) {
		case tafnuXScl:	NUSetXScl( a );	break;
		case tafnuYScl:	NUSetYScl( a );	break;
		case tafnuXSkw:	NUSetXSkw( a );	break;
		case tafnuYSkw:	NUSetYSkw( a );	break;
		case tafnuRot:	NUSetRot( a );	break;
		default:		NUSetScl( a );	break;
	}
}

/* --------------------------------------------------------------- */
/* ComposeDfm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Conventionalized compounding of deformations.
// Tdfm = Scale.(Yskew.Xskew)
//
void TAffine::ComposeDfm(
	double	scl,
	double	xscl,
	double	yscl,
	double	xskw,
	double	yskw )
{
	TAffine	S, Y, X;

	X.NUSetXSkw( xskw );
	Y.NUSetYSkw( yskw );

	S.t[0] = xscl*scl; S.t[1] = 0;        S.t[2] = 0;
	S.t[3] = 0;        S.t[4] = yscl*scl; S.t[5] = 0;

	*this = S * (Y * X);
}

/* --------------------------------------------------------------- */
/* SetCWRot ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Set a transform that performs a CW rotation of deg degrees
// about the given pivot point (remember +deg is CW).
//
// That is, A' = piv + R(A-piv) = R(A) + piv - R(piv).
//
void TAffine::SetCWRot( double deg, const Point &pivot )
{
	Point	Prot = pivot;

	NUSetRot( deg*PI/180 );
	Apply_R_Part( Prot );
	AddXY( pivot.x - Prot.x, pivot.y - Prot.y );
}

/* --------------------------------------------------------------- */
/* FromAToB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Compute TAffine [atob = binv*a] from global TForms a, b.
//
void TAffine::FromAToB( const TAffine &a, const TAffine &b )
{
	TAffine	binv;

	binv.InverseOf( b );
	*this = binv * a;
}

/* --------------------------------------------------------------- */
/* InverseOf ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void TAffine::InverseOf( const TAffine &a )
{
	double	det = 1 / a.det();

// simple inverse of matrix part
	t[0] =  a.t[4]*det;
	t[1] = -a.t[1]*det;
	t[3] = -a.t[3]*det;
	t[4] =  a.t[0]*det;

// apply inverse to translation and negate
	t[2] = -(t[0]*a.t[2] + t[1]*a.t[5]);
	t[5] = -(t[3]*a.t[2] + t[4]*a.t[5]);
}

/* --------------------------------------------------------------- */
/* operator_* ---------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	L(i)	this->t[i]
#define	R(i)	rhs.t[i]

TAffine TAffine::operator * ( const TAffine& rhs ) const
{
	TAffine	r;

	r.t[0] = L(0)*R(0) + L(1)*R(3);
	r.t[1] = L(0)*R(1) + L(1)*R(4);
	r.t[2] = L(0)*R(2) + L(1)*R(5) + L(2);
	r.t[3] = L(3)*R(0) + L(4)*R(3);
	r.t[4] = L(3)*R(1) + L(4)*R(4);
	r.t[5] = L(3)*R(2) + L(4)*R(5) + L(5);

	return r;
}

/* --------------------------------------------------------------- */
/* ScanTrackEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

void TAffine::ScanTrackEM2( const char *s )
{
	sscanf( s, "matrix(%lf,%lf,%lf,%lf,%lf,%lf",
		&t[0], &t[3], &t[1], &t[4], &t[2], &t[5] );
}

/* --------------------------------------------------------------- */
/* TPrint -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void TAffine::TPrint( FILE *f, const char *s ) const
{
	if( !f )
		f = stdout;

	fprintf( f, "%s%7.4f %7.4f %8.2f   %7.4f %7.4f %8.2f\n",
		(s ? s : ""), t[0], t[1], t[2], t[3], t[4], t[5] );
}

/* --------------------------------------------------------------- */
/* TPrintAsParam ------------------------------------------------- */
/* --------------------------------------------------------------- */

void TAffine::TPrintAsParam( FILE *f, bool newline ) const
{
	fprintf( f, " -TRA=%.4f,%.4f,%.2f,%.4f,%.4f,%.2f%c",
		t[0], t[1], t[2], t[3], t[4], t[5],
		(newline ? '\n' : ' ') );
}

/* --------------------------------------------------------------- */
/* EffArea ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return area of polygon obtained from: T( unit square ).
//
// Method follows Geometry::AreaOfPolygon().
//
double TAffine::EffArea() const
{
	vector<Point>	v;

	v.push_back( Point(0,0) );
	v.push_back( Point(1,0) );
	v.push_back( Point(1,1) );
	v.push_back( Point(0,1) );

	Transform( v );

	double	A = 0;

	for( int i = 0; i < 4; ++i ) {

		const Point&	a = v[i];
		const Point&	b = v[(i+1)%4];

		A += (a.x - b.x) * (a.y + b.y);
	}

	return fabs( A * 0.5 );
}

/* --------------------------------------------------------------- */
/* GetRadians ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// A general affine may include (ignoring reflection here),
// rotation, skew and scale. Using polar decomposition we
// attempt to extract the (roughly) pure rotation matrix R
// from general affine A and return its associated angle
// in radians.
//
// Method (adapted from "Graphics Gems IV", section III.4,
// by Ken Shoemake).
//
// From TAffine's matrix-part P = R( a ), iteratively compute
// N = (P + Trp(Inv(P)) / 2, until max element-wise change
// in (N - P) is tiny.
//
double TAffine::GetRadians() const
{
	double	P[4] = {t[0], t[1], t[3], t[4]};

	for( int iter = 0; iter < 100; ++iter ) {

		// N = Trp(Pinv)
		double	d = 1/(P[0]*P[3] - P[1]*P[2]);
		double	N[4] = {P[3]*d, -P[2]*d, -P[1]*d, P[0]*d};

		// N = (P + N)/2
		N[0] = (P[0] + N[0]) * 0.5;
		N[1] = (P[1] + N[1]) * 0.5;
		N[2] = (P[2] + N[2]) * 0.5;
		N[3] = (P[3] + N[3]) * 0.5;

		// max( N - P )
		d = fabs( N[0] - P[0] );
		d = fmax( d, fabs( N[1] - P[1] ) );
		d = fmax( d, fabs( N[2] - P[2] ) );
		d = fmax( d, fabs( N[3] - P[3] ) );

		// converged?
		if( d < 1e-7 ) {

			//printf( "i=%d, %f  %f  %f  %f\n",
			//iter, N[0], N[1], N[2], N[3] );

			return atan2( N[2], N[0] );
		}

		// P = N;
		P[0] = N[0]; P[1] = N[1]; P[2] = N[2]; P[3] = N[3];
	}

// fall back on standard crude estimator

	return atan2( t[3], t[0] );
}

/* --------------------------------------------------------------- */
/* Squareness ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Apply transform to a right angle and return |cosA|, which
// is the same as |sin(90-A)|, that is, return (sin of) the
// deviation from a right angle.
//
// Efficient version of same method in THmgphy class.
// Here, the two unit vectors are just the two columns
// of the R-part.
//
double TAffine::Squareness() const
{
	double	c = (t[0]*t[1] + t[3]*t[4]) /
		sqrt(   (t[0]*t[0] + t[3]*t[3]) *
		        (t[1]*t[1] + t[4]*t[4]) );

	return (c >= 0.0 ? c : -c);
}

/* --------------------------------------------------------------- */
/* Transform ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void TAffine::Transform( Point &p ) const
{
	double	x = p.x*t[0] + p.y*t[1] + t[2];
	double	y = p.x*t[3] + p.y*t[4] + t[5];

	p.x = x;
	p.y = y;
}


void TAffine::Transform( vector<Point> &v ) const
{
	int		N = v.size();

	for( int i = 0; i < N; ++i )
		Transform( v[i] );
}

/* --------------------------------------------------------------- */
/* Apply_R_Part -------------------------------------------------- */
/* --------------------------------------------------------------- */

// If transform T(p) is represented as R(p) + V,
// then here we apply R(p) without V.
//
void TAffine::Apply_R_Part( Point &p ) const
{
	double	x = p.x*t[0] + p.y*t[1];
	double	y = p.x*t[3] + p.y*t[4];

	p.x = x;
	p.y = y;
}


void TAffine::Apply_R_Part( vector<Point> &v ) const
{
	int		N = v.size();

	for( int i = 0; i < N; ++i )
		Apply_R_Part( v[i] );
}


