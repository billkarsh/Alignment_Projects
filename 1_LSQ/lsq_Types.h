

#pragma once


#include	"lsq_DIR.h"

#include	"PipeFiles.h"
#include	"CPoint.h"


/* --------------------------------------------------------------- */
/* RGN ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Connected region

class RGN {

private:
	int	iname;	// index in rgnvname vector

public:
	int	z,		// layer
		id,		// index in layer
		rgn,	// connRgn index
		itr;	// global-space transform index

public:
	RGN( const char *key );
	RGN( const char *path, const DIR &dir, int _id );
	RGN( const char *path, const char *key );

	const char* GetName() const;

	static void GetMeta(
		const Til2Img*	*t1,
		const Til2Img*	*t2,
		const RGN		&I1,
		const RGN		&I2 );
};

/* --------------------------------------------------------------- */
/* zsort --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Orders regions by z, for writing files

class zsort{
public:
	int	z, id, rgn, i;
public:
	zsort()	{};
	zsort( const RGN& R, int i )
	: z(R.z), id(R.id), rgn(R.rgn), i(i) {};

	bool operator < (const zsort &rhs) const
		{
			if( z < rhs.z )
				return true;
			if( z > rhs.z )
				return false;
			if( id < rhs.id )
				return true;
			if( id > rhs.id )
				return false;

			return rgn < rhs.rgn;
		};
};

/* --------------------------------------------------------------- */
/* CRPair -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Pair of matched connRgn

class CRPair {
public:
	int	a, b;   // always stored so a <= b
public:
	CRPair( int aa, int bb )
		{
			if( aa <= bb )
				{a = aa; b = bb;}
			else
				{a = bb; b = aa;}
		};

	bool operator < (const CRPair &rhs) const
		{return a < rhs.a || (a == rhs.a && b < rhs.b);};

	bool operator == (const CRPair &rhs) const
		{return a == rhs.a && b == rhs.b;};
};

/* --------------------------------------------------------------- */
/* Constraint ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Any matched point pair (POINT entry)

class Constraint {

public:
	Point	p1, p2;	// the two points
	int		r1, r2;	// indexes to the two regions
	bool	used,	// is this constraint used?
			inlier;	// is this a RANSAC inlier or outlier?

public:
	Constraint(
		int				r1,
		const Point		&p1,
		int				r2,
		const Point		&p2 )
		: p1(p1), p2(p2), r1(r1), r2(r2) {};
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern string				idb;		// for name lookups
extern vector<RGN>			vRgn;		// the regions
extern map<CRPair,int>		r12Idx;		// idx from region-pair
extern vector<Constraint>	vAllC;		// all Point-pairs


