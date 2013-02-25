

#include	"lsq_Rigid.h"
#include	"lsq_Types.h"


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MRigid::SetPointPairs(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc,
	double			same_strength )
{
	int	nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		double	fz =
		(vRgn[C.r1].z == vRgn[C.r2].z ? same_strength : 1);

		double	x1 = C.p1.x * fz / sc,
				y1 = C.p1.y * fz / sc,
				x2 = C.p2.x * fz / sc,
				y2 = C.p2.y * fz / sc;
		int		j  = vRgn[C.r1].itr * NX,
				k  = vRgn[C.r2].itr * NX;

		// T1(p1) - T2(p2) = 0
		//
		// {X0,X1,X2,X3} = {c,s,x,y}

		double	v[6]  = { x1, -y1,  fz, -x2,  y2, -fz};
		int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
		int		i2[6] = {j+1,   j, j+3, k+1,   k, k+3};

		AddConstraint( LHS, RHS, 6, i1, v, 0.0 );

		v[1] =  y1;
		v[4] = -y2;

		AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
	}
}


