


#pragma once


#include	"CPoint.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum TFNearUnityConsts {
// near unity transform selectors
	tfnuScl		= 0,
	tfnuXScl	= 1,
	tfnuYScl	= 2,
	tfnuXSkw	= 3,
	tfnuYSkw	= 4,
	tfnuRot		= 5
};

/* --------------------------------------------------------------- */
/* class TForm --------------------------------------------------- */
/* --------------------------------------------------------------- */

class TForm {

public:
	double	t[6];

public:
	TForm()
		{t[0]=1.0; t[1]=0.0; t[2]=0.0;
		 t[3]=0.0; t[4]=1.0; t[5]=0.0;};

	TForm( const TForm &rhs )
		{CopyIn( rhs );};

	TForm( const double *src )
		{CopyIn( src );};

	TForm( double a, double b, double c, double d, double e, double f )
		{t[0]=a; t[1]=b; t[2]=c;
		 t[3]=d; t[4]=e; t[5]=f;};

	void CopyIn( const TForm &src )
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

	void NUSetScl( double a );
	void NUSetXScl( double a );
	void NUSetYScl( double a );
	void NUSetXSkw( double a );
	void NUSetYSkw( double a );
	void NUSetRot( double r );

	void NUSelect( int sel, double a );

	void Transform( Point &p ) const;
	void Transform( vector<Point> &v ) const;

	void Apply_R_Part( Point &p ) const;
	void Apply_R_Part( vector<Point> &v ) const;

	void RotateAround( Point s, Point tar, double rad );

	void ScanTrackEM2( const char *s );

	void WriteTransform( FILE *f, const char *s ) const;
	void PrintTransform( FILE *f = NULL ) const;
	void PrintTransformAsParam( FILE *f, bool newline = false ) const;
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void InvertTrans( TForm &inv, const TForm &t );
void MultiplyTrans( TForm &r, const TForm &a, const TForm &b );
void AToBTrans( TForm &atob, const TForm &a, const TForm &b );
void CreateCWRot( TForm &T, double deg, const Point &pivot );
double RadiansFromAffine( const TForm &a );


