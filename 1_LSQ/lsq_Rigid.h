

#pragma once


#include	"lsq_MDL.h"


/* --------------------------------------------------------------- */
/* Rigid MDL ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class MRigid : public MDL {

public:
	MRigid() : MDL( 6, 4 ) {};

	void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength );
};


