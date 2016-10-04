

#pragma once


#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* class CRigid -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CRigid {
protected:
    double	Xa, Ya, Xb, Yb, XaXb, YaYb, XaYb, YaXb;
    int		N;
public:
    CRigid()
    : Xa(0), Ya(0), Xb(0), Yb(0),
      XaXb(0), YaYb(0), XaYb(0), YaXb(0), N(0)
    {};
    virtual void Add( const Point &A, const Point &B );
    virtual void Solve( TAffine& T );
    void Regularize( double *v, int nv, double Wr );
};

/* --------------------------------------------------------------- */
/* class CTrans -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CTrans : public CRigid {
public:
    void Add( const Point &A, const Point &B );
    void Solve( TAffine& T );
};


