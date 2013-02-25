

#pragma once


#include	"lsq_MDL.h"


/* --------------------------------------------------------------- */
/* Affine MDL ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class MAffine : public MDL {

public:
	MAffine() : MDL( 6, 6 ) {};

	void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength );
};


