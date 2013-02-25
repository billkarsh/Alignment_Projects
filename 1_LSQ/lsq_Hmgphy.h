

#pragma once


#include	"lsq_MDL.h"


/* --------------------------------------------------------------- */
/* Homography MDL ------------------------------------------------ */
/* --------------------------------------------------------------- */

class MHmgphy : public MDL {

public:
	MHmgphy() : MDL( 8, 8 ) {};

	void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength );
};


