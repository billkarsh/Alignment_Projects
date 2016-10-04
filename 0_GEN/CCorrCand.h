

#pragma once


#include	"GenDefs.h"


/* --------------------------------------------------------------- */
/* class CorrCand ------------------------------------------------ */
/* --------------------------------------------------------------- */

class CorrCand {	// a correlation candidate

public:
    double	val;
    int		x, y;

public:
    CorrCand()
    : x(0), y(0), val(-BIG) {};

    CorrCand( int x, int y, double val )
    : x(x), y(y), val(val) {};

    bool operator < ( const CorrCand &rhs ) const
        {return val < rhs.val;};
};


