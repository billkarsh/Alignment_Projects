

#pragma once


// basic types
#define _USESTDVECTOR_ 1
#include	"nr3.h"
#include	"fitmrq.h"


// function prototypes
void gaussj( MatDoub_IO &a, MatDoub_IO &b );
void fgauss( Doub x, VecDoub_I &a, Doub &y, VecDoub_O &dyda );
Doub ygauss( Doub x, VecDoub_I &a );


