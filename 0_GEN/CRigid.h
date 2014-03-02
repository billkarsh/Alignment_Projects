

#pragma once


#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* class CRigid -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CRigid {
public:
	double	Xa, Ya, Xb, Yb, XaXb, YaYb, XaYb, YaXb;
	int		N;
public:
	CRigid();
	void Add( const Point &A, const Point &B );
	void Solve( TAffine& T );
};


