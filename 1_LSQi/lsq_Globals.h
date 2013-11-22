

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
public:
	vector<char>	used;
	map<int,int>	m;
	int				nr,
					z,
					cached_id,
					cached_i;
public:
	Rgns( int z );
	int Map( int id, int r );
};

class CorrPnt {
public:
	Point	p1, p2;
	int		z1, z2,
			i1, i2;
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

void MapZPair( int &ia, int &ib, int za, int zb );


