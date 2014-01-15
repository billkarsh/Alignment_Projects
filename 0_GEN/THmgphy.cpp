

#include	"GenDefs.h"
#include	"THmgphy.h"

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

	In our code, the 8 homography tform elements are interpreted
	as follows:

		u    |t0  t1  t2|    x
		     |          |          x'     u / w
		v  = |t3  t4  t5| *  y ;       =
		     |          |          y'     v / w
		w    |t6  t7   1|    1

	In ImageJ (including TrakEM2) the elements in source code
	are arranged by columns.

	In a TrackEM2 XML file transform attributes are written
	"matrix(m00,m10,m20,m01,m11,m21,m02,m12)".

	We use a standard order for compounding operations. Let
	X=xskew, Y=yskew, S=scaling, P=pretweak, R=rot+trans.
	Then, D = S(YX), and T = R(DP).

	Homographies are a superset of the affines, and reduce
	to affines when {t6 t7 t8} are {0 0 1}.
*/

/* --------------------------------------------------------------- */
/* Near unity transforms ----------------------------------------- */
/* --------------------------------------------------------------- */

// pos a enlarges
void THmgphy::NUSetScl( double a )
{
	t[0] = a; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = a; t[5] = 0;
	Zero67();
}


// pos a enlarges
void THmgphy::NUSetXScl( double a )
{
	t[0] = a; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = 1; t[5] = 0;
	Zero67();
}


// pos a enlarges
void THmgphy::NUSetYScl( double a )
{
	t[0] = 1; t[1] = 0; t[2] = 0;
	t[3] = 0; t[4] = a; t[5] = 0;
	Zero67();
}


// pos a tips right
void THmgphy::NUSetXSkw( double a )
{
	t[0] = 1; t[1] = a; t[2] = 0;
	t[3] = 0; t[4] = 1; t[5] = 0;
	Zero67();
}


// pos a tips up
void THmgphy::NUSetYSkw( double a )
{
	t[0] = 1; t[1] = 0; t[2] = 0;
	t[3] = a; t[4] = 1; t[5] = 0;
	Zero67();
}


// pos r rotates CW
void THmgphy::NUSetRot( double r )
{
	double	c = cos( r ), s = sin( r );

	t[0] = c; t[1] = -s; t[2] = 0;
	t[3] = s; t[4] =  c; t[5] = 0;
	Zero67();
}


void THmgphy::NUSelect( int sel, double a )
{
	switch( sel ) {
		case thgnuXScl:	NUSetXScl( a );	break;
		case thgnuYScl:	NUSetYScl( a );	break;
		case thgnuXSkw:	NUSetXSkw( a );	break;
		case thgnuYSkw:	NUSetYSkw( a );	break;
		case thgnuRot:	NUSetRot( a );	break;
		default:		NUSetScl( a );	break;
	}
}

/* --------------------------------------------------------------- */
/* ComposeDfm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Conventionalized compounding of deformations.
// Tdfm = Scale.(Yskew.Xskew)
//
void THmgphy::ComposeDfm(
	double	scl,
	double	xscl,
	double	yscl,
	double	xskw,
	double	yskw )
{
	THmgphy	S, Y, X;

	X.NUSetXSkw( xskw );
	Y.NUSetYSkw( yskw );

	S.t[0] = xscl*scl; S.t[1] = 0;        S.t[2] = 0;
	S.t[3] = 0;        S.t[4] = yscl*scl; S.t[5] = 0;
	S.Zero67();

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
void THmgphy::SetCWRot( double deg, const Point &pivot )
{
	Point	Prot = pivot;

	NUSetRot( deg*PI/180 );
	Transform( Prot );
	AddXY( pivot.x - Prot.x, pivot.y - Prot.y );
	Zero67();
}

/* --------------------------------------------------------------- */
/* FromAToB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Compute THmgphy [atob = binv*a] from global TForms a, b.
//
void THmgphy::FromAToB( const THmgphy &a, const THmgphy &b )
{
	THmgphy	binv;

	binv.InverseOf( b );
	*this = binv * a;
}

/* --------------------------------------------------------------- */
/* InverseOf ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void THmgphy::InverseOf( const THmgphy &a )
{
	double	t8;

	t[0] =  (a.t[4]        - a.t[5]*a.t[7]);
	t[3] = -(a.t[3]        - a.t[5]*a.t[6]);
	t[6] =  (a.t[3]*a.t[7] - a.t[4]*a.t[6]);
	t[1] = -(a.t[1]        - a.t[2]*a.t[7]);
	t[4] =  (a.t[0]        - a.t[2]*a.t[6]);
	t[7] = -(a.t[0]*a.t[7] - a.t[1]*a.t[6]);
	t[2] =  (a.t[1]*a.t[5] - a.t[2]*a.t[4]);
	t[5] = -(a.t[0]*a.t[5] - a.t[2]*a.t[3]);
	t8   =  (a.t[0]*a.t[4] - a.t[1]*a.t[3]);

	t8 = 1 / t8;

	for( int i = 0; i < 8; ++i )
		t[i] *= t8;
}

/* --------------------------------------------------------------- */
/* operator_* ---------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	L(i)	this->t[i]
#define	R(i)	rhs.t[i]

THmgphy THmgphy::operator * ( const THmgphy& rhs ) const
{
	THmgphy	r;
	double	t8;

	r.t[0] = L(0)*R(0) + L(1)*R(3) + L(2)*R(6);
	r.t[1] = L(0)*R(1) + L(1)*R(4) + L(2)*R(7);
	r.t[2] = L(0)*R(2) + L(1)*R(5) + L(2);
	r.t[3] = L(3)*R(0) + L(4)*R(3) + L(5)*R(6);
	r.t[4] = L(3)*R(1) + L(4)*R(4) + L(5)*R(7);
	r.t[5] = L(3)*R(2) + L(4)*R(5) + L(5);
	r.t[6] = L(6)*R(0) + L(7)*R(3) +      R(6);
	r.t[7] = L(6)*R(1) + L(7)*R(4) +      R(7);
	t8     = L(6)*R(2) + L(7)*R(5) +         1;

	t8 = 1 / t8;

	for( int i = 0; i < 8; ++i )
		r.t[i] *= t8;

	return r;
}

/* --------------------------------------------------------------- */
/* ScanTrackEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

void THmgphy::ScanTrackEM2( const char *s )
{
	sscanf( s, "matrix(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		&t[0], &t[3], &t[6], &t[1], &t[4], &t[7], &t[2], &t[5] );
}

/* --------------------------------------------------------------- */
/* TPrint -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void THmgphy::TPrint( FILE *f, const char *s ) const
{
	if( !f )
		f = stdout;

	fprintf( f, "%s%f %f %f  %f %f %f  %.12g %.12g\n",
		(s ? s : ""),
		t[0], t[1], t[2],
		t[3], t[4], t[5],
		t[6], t[7] );
}

/* --------------------------------------------------------------- */
/* TPrintAsParam ------------------------------------------------- */
/* --------------------------------------------------------------- */

void THmgphy::TPrintAsParam( FILE *f, bool newline ) const
{
	fprintf( f, " -TRH=%f,%f,%f,%f,%f,%f,%.12g,%.12g%c",
		t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
		(newline ? '\n' : ' ') );
}

/* --------------------------------------------------------------- */
/* EffArea ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return area of polygon obtained from: T( unit square ).
//
// Method follows Geometry::AreaOfPolygon().
//
double THmgphy::EffArea() const
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
// From THmgphy's matrix-part P = R( a ), iteratively compute
// N = (P + Trp(Inv(P)) / 2, until max element-wise change
// in (N - P) is tiny.
//
double THmgphy::GetRadians() const
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
double THmgphy::Squareness() const
{
	Point	po, px( 1, 0 ), py( 0, 1 );

	Transform( po );
	Transform( px );
	Transform( py );

	px.x -= po.x;
	py.x -= po.x;
	px.y -= po.y;
	py.y -= po.y;

	double	c = (px.x*py.x + px.y*py.y) /
		sqrt(   (px.x*px.x + px.y*px.y) *
		        (py.x*py.x + py.y*py.y) );

	return (c >= 0.0 ? c : -c);
}

/* --------------------------------------------------------------- */
/* Transform ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void THmgphy::Transform( Point &p ) const
{
	double	u = p.x*t[0] + p.y*t[1] + t[2];
	double	v = p.x*t[3] + p.y*t[4] + t[5];
	double	w = p.x*t[6] + p.y*t[7] + 1;

	p.x = u / w;
	p.y = v / w;
}


void THmgphy::Transform( vector<Point> &v ) const
{
	int		N = v.size();

	for( int i = 0; i < N; ++i )
		Transform( v[i] );
}


