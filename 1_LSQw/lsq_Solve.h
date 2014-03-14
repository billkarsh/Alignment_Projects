

#pragma once


#include	"lsq_XArray.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void SetSolveParams( int type, double inWr, double inEtol );

void Solve( XArray &Xsrc, XArray &Xdst, int iters );


