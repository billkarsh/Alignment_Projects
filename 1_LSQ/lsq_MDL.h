

#pragma once


#include	"LinEqu.h"


/* --------------------------------------------------------------- */
/* MDL base classs ----------------------------------------------- */
/* --------------------------------------------------------------- */

class MDL {

protected:
	const int	NT, NX;

public:
	MDL( int NT, int NX ) : NT(NT), NX(NX) {};

	int MinLinks()	{return NX/2;};

	virtual void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength ) = 0;
};


