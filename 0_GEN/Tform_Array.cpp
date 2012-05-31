

#include	"Tform_Array.h"

#include	<math.h>
#include	<stdio.h>






/* --------------------------------------------------------------- */
/* Transform ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void Transform( Point &p, const atform &t )
{
	double	x = p.x*t[0] + p.y*t[1] + t[2];
	double	y = p.x*t[3] + p.y*t[4] + t[5];

	p.x = x;
	p.y = y;
}


void Transform( vector<Point> &v, const atform &t )
{
	int		N = v.size();

	for( int i=0; i < N; ++i )
		Transform( v[i], t );
}

/* --------------------------------------------------------------- */
/* InvertTrans --------------------------------------------------- */
/* --------------------------------------------------------------- */

void InvertTrans( atform &inv, const atform &t )
{
	double	det = t[0]*t[4] - t[1]*t[3];

	inv[0] = t[4]/det;
	inv[1] = -t[1]/det;
	inv[2] = (t[5]*t[1]-t[4]*t[2])/det;
	inv[3] = -t[3]/det;
	inv[4] = t[0]/det;
	inv[5] = (t[2]*t[3]-t[0]*t[5])/det;
}

/* --------------------------------------------------------------- */
/* MultiplyTrans ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Compute r = a*b
//
void MultiplyTrans( atform &r, const atform &a, const atform &b )
{
	r[0] = a[0]*b[0]+a[1]*b[3];
	r[1] = a[0]*b[1]+a[1]*b[4];
	r[2] = a[0]*b[2]+a[1]*b[5]+a[2];
	r[3] = a[3]*b[0]+a[4]*b[3];
	r[4] = a[3]*b[1]+a[4]*b[4];
	r[5] = a[3]*b[2]+a[4]*b[5]+a[5];
}

/* --------------------------------------------------------------- */
/* RotateAround -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fiddle with the rotation, but still map point s to point 'tar'
//
void RotateAround( atform &t, Point s, Point tar, double rad )
{
	double	co = cos(rad), si = sin(rad);
	double	a = t[0], b=t[1], c=t[3], d=t[4];

	t[0] = a*co    + b*si;
	t[1] = a*(-si) + b*co;
	t[3] = c*co    + d*si;
	t[4] = c*(-si) + d*co;

// now make point s come out where it originally did
	Point test(s.x, s.y);
	Transform( test, t );
	t[2] += tar.x - test.x;
	t[5] += tar.y - test.y;
}

/* --------------------------------------------------------------- */
/* WriteTransform ------------------------------------------------ */
/* --------------------------------------------------------------- */

void WriteTransform( const char *s, const atform &t )
{
	printf( "%s %9.4f %9.4f %10.2f\n%s %9.4f %9.4f %10.2f\n",
		s, t[0], t[1], t[2], s, t[3], t[4], t[5] );
}

/* --------------------------------------------------------------- */
/* PrintTransform ------------------------------------------------ */
/* --------------------------------------------------------------- */

void PrintTransform( const atform &t )
{
	printf( "%7.4f %7.4f %8.2f   %7.4f %7.4f %8.2f\n",
		t[0], t[1], t[2], t[3], t[4], t[5] );
}


