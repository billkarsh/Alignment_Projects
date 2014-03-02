

#include	"CRigid.h"

#include	<math.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* Discussion ---------------------------------------------------- */
/* --------------------------------------------------------------- */

/*
	Solving for a rigid transform from A to B...

	Let c=cos(theta), s=sin(theta), Sum over all point-pairs:

	E = Sum[Xb - cXa + sYa - kx]^2 + [Yb - sXa - cYa - ky]^2

	The params u = {theta,kx,ky) are determined by dE/du = 0.
	If we use notation [arg] => Sum[argi] over point-pairs,

	kx = ([Xb] - c[Xa] + s[Ya]) / N
	ky = ([Yb] - s[Xa] - c[Ya]) / N

				 [YaXb] - [XaYb] + ([Xa][Yb] - [Ya][Xb]) / N
	tan(theta) = --------------------------------------------
				 ([Xa][Xb] + [Ya][Yb]) / N - [XaXb] - [YaYb]
*/

/* --------------------------------------------------------------- */
/* CRigid::CRigid ------------------------------------------------ */
/* --------------------------------------------------------------- */

CRigid::CRigid()
{
	memset( this, 0, sizeof(CRigid) );
}

/* --------------------------------------------------------------- */
/* CRigid::Add --------------------------------------------------- */
/* --------------------------------------------------------------- */

void CRigid::Add( const Point &A, const Point &B )
{
	Xa += A.x;
	Ya += A.y;
	Xb += B.x;
	Yb += B.y;

	XaXb += A.x * B.x;
	YaYb += A.y * B.y;
	XaYb += A.x * B.y;
	YaXb += A.y * B.x;

	++N;
}

/* --------------------------------------------------------------- */
/* CRigid::Solve ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CRigid::Solve( TAffine& T )
{
	if( N >= 2 ) {

		double
		theta = atan(
			(YaXb - XaYb + (Xa*Yb - Ya*Xb)/N) /
			((Xa*Xb + Ya*Yb)/N - XaXb - YaYb)
		),
		c  = cos( theta ),
		s  = sin( theta ),
		kx = (Xb - c*Xa + s*Ya) / N,
		ky = (Yb - s*Xa - c*Ya) / N;

		T.t[0] = c; T.t[1] = -s; T.t[2] = kx;
		T.t[3] = s; T.t[4] =  c; T.t[5] = ky;
	}
	else
		T.NUSetOne();
}


