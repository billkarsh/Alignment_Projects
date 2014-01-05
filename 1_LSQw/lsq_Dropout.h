

#pragma once


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Dropout {
public:
	long	iter,	// dropped as result of iterating
			pnts;	// dropped as result of initial corr-pnts
private:
	void Add( const Dropout& rhs )
		{iter += rhs.iter; pnts += rhs.pnts;};
	void GatherCounts();
public:
	Dropout() : iter(0), pnts(0) {};
	void Scan();
};


