

#pragma once


#include	"CPoint.h"


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef double atform[6];

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void Transform( Point &p, const atform &t );
void Transform( vector<Point> &v, const atform &t );
void InvertTrans( atform &inv, const atform &t );
void MultiplyTrans( atform &r, const atform &a, const atform &b );
void RotateAround( atform &t, Point s, Point tar, double rad );

void WriteTransform( const char *s, const atform &t );
void PrintTransform( const atform &t );


