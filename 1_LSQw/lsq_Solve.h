

#pragma once


#include	"lsq_XArray.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void SetRegularizer( int type, double weight );

void Solve( XArray &Xsrc, XArray &Xdst, int iters );


