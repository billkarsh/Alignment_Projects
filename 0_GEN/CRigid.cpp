

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

    Translation is merely a simplification of the rigid such
    that {c,s} = (1,0), hence, no cross terms are needed.
*/

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

/* --------------------------------------------------------------- */
/* CRigid::Regularize -------------------------------------------- */
/* --------------------------------------------------------------- */

// Solve for a rigid transform (R) and form a weighted average
// of R and given values describing an affine or homography (T):
// T -> (1-Wr)*T + Wr*R.
//
// Number of transform values (nv) must be 6 or 8.
//
void CRigid::Regularize( double *v, int nv, double Wr )
{
    if( N < 2 || Wr == 0 )
        return;

    TAffine	R;
    Solve( R );

    double	Wv = (1.0 - Wr);

    for( int i = 0; i < 6; ++i )
        v[i] = Wv * v[i] + Wr * R.t[i];

    if( nv == 8 ) {
        v[6] *= Wv;
        v[7] *= Wv;
    }
}

/* --------------------------------------------------------------- */
/* CTrans::Add --------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTrans::Add( const Point &A, const Point &B )
{
    Xa += A.x;
    Ya += A.y;
    Xb += B.x;
    Yb += B.y;

    ++N;
}

/* --------------------------------------------------------------- */
/* CTrans::Solve ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTrans::Solve( TAffine& T )
{
    if( N >= 1 ) {
        T.t[0] = 1; T.t[1] = 0; T.t[2] = (Xb - Xa) / N;
        T.t[3] = 0; T.t[4] = 1; T.t[5] = (Yb - Ya) / N;
    }
    else
        T.NUSetOne();
}


