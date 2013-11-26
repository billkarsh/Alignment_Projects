

#pragma once


#include	"lsq_Layers.h"

#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stdio.h>

#include	<map>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Rgns {
// The rgns for given layer
// indexed by 0-based 'idx0'
public:
	vector<vector<int> >	pts;	// ea rgn's pts
	vector<char>			used;	// which rgns used
	map<int,int>			m;		// map id -> idx0
	int						nr,		// num rgns
							z,		// common z
							cached_id,
							cached_i;
public:
	Rgns( int z );
	int Map( int id, int r );
};

class CorrPnt {
public:
	Point	p1, p2;
	int		z1, z2,	// index into vR
			i1, i2;	// index into vR[iz]
	union {
	int		used;
	struct {
	uint16	r1, r2;
	};
	};
public:
	inline bool FromFile( FILE *f )
	{
		return 10 == fscanf( f, "CPOINT2"
			" %d.%d:%hd %lf %lf"
			" %d.%d:%hd %lf %lf\n",
			&z1, &i1, &r1, &p1.x, &p1.y,
			&z2, &i2, &r2, &p2.x, &p2.y );
	};
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern string			idb;
extern vector<Layer>	vL;
extern map<int,int>		mZ;
extern vector<Rgns>		vR;
extern vector<CorrPnt>	vC;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetIDB( const char *tempdir );
void InitTablesToMaximum();
void MapZPair( int &ia, int &ib, int za, int zb );
void RemapIndices();


