

#pragma once


#include	"../1_LSQi/lsq_Layers.h"

#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stdio.h>

#include	<map>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

 enum RgnFlags {
	fbDead		= 0x01,
	fbRead		= 0x02,
	fbPnts		= 0x04,
	fbKill		= 0x08,
	fbCutd		= 0x10,

	fmRead		= fbRead + fbDead,
	fmPnts		= fbPnts + fbDead,
	fmKill		= fbKill + fbDead,
 };

#define	FLAG_ISUSED( f )	(((f) & fbDead) == 0)
#define	FLAG_ISREAD( f )	(((f) & fbRead) != 0)
#define	FLAG_ISPNTS( f )	(((f) & fbPnts) != 0)
#define	FLAG_ISKILL( f )	(((f) & fbKill) != 0)
#define	FLAG_ISCUTD( f )	(((f) & fbCutd) != 0)

#define	FLAG_SETUSED( f )	(f = 0)
#define	FLAG_SETPNTS( f )	(f = fmPnts)
#define	FLAG_SETDISK( f )	(f = (FLAG_ISUSED(f) ? 0 : fmRead))

#define	FLAG_ADDPNTS( f )	(f |= fmPnts)
#define	FLAG_ADDKILL( f )	(f |= fmKill)
#define	FLAG_ADDCUTD( f )	(f |= fbCutd)

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Rgns {
// The rgns for given layer
// indexed by 0-based 'idx0'
public:
	vector<vector<int> >	pts;	// ea rgn's pts
	vector<uint8>			flag;	// rgn flags
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
extern int				zLlo, zLhi,	// LHS neib needs these
						zRlo, zRhi,	// RHS neib needs these
						zilo, zihi,
						zolo, zohi,	// index into vR
						gW,   gH,	// image dims
						maxthreads;
extern vector<Layer>	vL;
extern map<int,int>		mZ;
extern vector<Rgns>		vR;
extern vector<CorrPnt>	vC;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetIDB( const char *tempdir );
void InitTables( int argzilo, int argzihi );
bool MapZPair( int &ia, int &ib, int za, int zb );
void RemapIndices();
void RealZIDR( int &z, int &id, int &r, int iz, int idx );


