

#pragma once


#include	"CPoint.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum THmgphyNearUnityConsts {
// near unity transform selectors
	thgnuScl	= 0,
	thgnuXScl	= 1,
	thgnuYScl	= 2,
	thgnuXSkw	= 3,
	thgnuYSkw	= 4,
	thgnuRot	= 5
};

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given packed array of doubles (X), with 8 elems per hmgphy,
// interpret the ith set of eight as the ith hmgphy.
//
#define	X_AS_HMY( X, i )	((THmgphy*)&(X)[0])[i]

/* --------------------------------------------------------------- */
/* class THmgphy ------------------------------------------------- */
/* --------------------------------------------------------------- */

class THmgphy {

public:
	double	t[8];

public:
	THmgphy()
		{NUSetOne();};

	THmgphy( const THmgphy &rhs )
		{CopyIn( rhs );};

	THmgphy( const double *src )
		{CopyIn( src );};

	THmgphy( double a, double b, double c,
			 double d, double e, double f,
			 double g, double h )
		{t[0]=a; t[1]=b; t[2]=c;
		 t[3]=d; t[4]=e; t[5]=f;
		 t[6]=g; t[7]=h;};

	void CopyIn( const THmgphy &src )
		{for( int i = 0; i < 8; ++i ) t[i] = src.t[i];};

	void CopyIn( const double *src )
		{for( int i = 0; i < 8; ++i ) t[i] = src[i];};

	void CopyOut( double *dst ) const
		{for( int i = 0; i < 8; ++i ) dst[i] = t[i];};

	inline void Zero67()
		{t[6] = 0.0; t[7] = 0.0;};

	inline void SetXY( double dx, double dy )
		{t[2] = dx; t[5] = dy;};

	inline void AddXY( double dx, double dy )
		{t[2] += dx; t[5] += dy;};

	inline void MulXY( double f )
		{t[2] *= f; t[5] *= f;};

	double det() const
		{return t[0]*t[4]
			  + t[1]*t[5]*t[6]
			  + t[3]*t[7]*t[2]
			  - t[2]*t[4]*t[6]
			  - t[1]*t[3]
			  - t[5]*t[7]*t[0];};

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

	void FromAToB( const THmgphy &a, const THmgphy &b );

	void InverseOf( const THmgphy &a );

	THmgphy operator * ( const THmgphy& rhs ) const;

	void ScanTrackEM2( const char *s );

	void TPrint( FILE *f = NULL, const char *s = NULL ) const;
	void TPrintAsParam( FILE *f, bool newline = false ) const;

	double EffArea() const;
	double GetRadians() const;
	double Squareness() const;

	void Transform( Point &p ) const;
	void Transform( vector<Point> &v ) const;
};


