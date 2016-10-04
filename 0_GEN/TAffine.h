

#pragma once


#include	"CPoint.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum TAffineNearUnityConsts {
// near unity transform selectors
    tafnuScl	= 0,
    tafnuXScl	= 1,
    tafnuYScl	= 2,
    tafnuXSkw	= 3,
    tafnuYSkw	= 4,
    tafnuRot	= 5
};

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given packed array of doubles (X), with 6 elems per affine,
// interpret the ith set of six as the ith affine.
//
#define	X_AS_AFF( X, i )	((TAffine*)&(X)[0])[i]

/* --------------------------------------------------------------- */
/* class TAffine ------------------------------------------------- */
/* --------------------------------------------------------------- */

class TAffine {

public:
    double	t[6];

public:
    TAffine()
        {NUSetOne();};

    TAffine( const TAffine &rhs )
        {CopyIn( rhs );};

    TAffine( const double *src )
        {CopyIn( src );};

    TAffine( double a, double b, double c,
             double d, double e, double f )
        {t[0]=a; t[1]=b; t[2]=c;
         t[3]=d; t[4]=e; t[5]=f;};

    void CopyIn( const TAffine &src )
        {for( int i = 0; i < 6; ++i ) t[i] = src.t[i];};

    void CopyIn( const double *src )
        {for( int i = 0; i < 6; ++i ) t[i] = src[i];};

    void CopyOut( double *dst ) const
        {for( int i = 0; i < 6; ++i ) dst[i] = t[i];};

    inline void SetXY( double dx, double dy )
        {t[2] = dx; t[5] = dy;};

    inline void AddXY( double dx, double dy )
        {t[2] += dx; t[5] += dy;};

    inline void MulXY( double f )
        {t[2] *= f; t[5] *= f;};

    void ToMatlab()
        {double b=t[1], c=t[2], d=t[3], e=t[4];
         t[1]=d; t[2]=b; t[3]=e; t[4]=c;};

    void FromMatlab()
        {double b=t[2], c=t[4], d=t[1], e=t[3];
         t[1]=b; t[2]=c; t[3]=d; t[4]=e;};

    double det() const
        {return t[0]*t[4] - t[1]*t[3];};

    void NUSetOne()
        {NUSetScl( 1.0 );};

    void NUSetScl( double a );
    void NUSetXScl( double a );
    void NUSetYScl( double a );
    void NUSetXSkw( double a );
    void NUSetYSkw( double a );
    void NUSetRot( double r );

    void NUSelect( int sel, double a );

    void ComposeDfm(
        double	scl,
        double	xscl,
        double	yscl,
        double	xskw,
        double	yskw );

    void SetCWRot( double deg, const Point &pivot );

    void FromAToB( const TAffine &a, const TAffine &b );

    void InverseOf( const TAffine &a );

    TAffine operator * ( const TAffine& rhs ) const;

    void ScanTrackEM2( const char *s );

    void TPrint( FILE *f = NULL, const char *s = NULL ) const;
    void TPrintAsParam( FILE *f, bool newline = false ) const;

    double EffArea() const;
    double GetRadians() const;
    double Squareness() const;

    void Transform( Point &p ) const;
    void Transform( vector<Point> &v ) const;

    void Apply_R_Part( Point &p ) const;
    void Apply_R_Part( vector<Point> &v ) const;
};


