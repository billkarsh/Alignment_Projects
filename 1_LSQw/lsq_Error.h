

#pragma once


#include	"lsq_XArray.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void Error( const XArray &X, double inEtol );

void GetFinalError( double &erms, double &emax );


