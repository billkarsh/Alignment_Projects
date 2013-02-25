

#include	"lsq_Hmgphy.h"
#include	"lsq_Types.h"


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SetPointPairs(
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

		double	va[6] = { x1,  y1,  fz, -x2, -y2, -fz};
		int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
		int		i2[6] = {j+3, j+4, j+5, k+3, k+4, k+5};

		AddConstraint( LHS, RHS, 6, i1, va, 0.0 );
		AddConstraint( LHS, RHS, 6, i2, va, 0.0 );

		double	vb[4] = { x1,  y1, -x2, -y2};
		int		i3[4] = {j+6, j+7, k+6, k+7};

		AddConstraint( LHS, RHS, 4, i3, vb, 0.0 );
	}
}


