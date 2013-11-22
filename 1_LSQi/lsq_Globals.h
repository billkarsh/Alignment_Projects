

#pragma once


#include	"lsq_Layers.h"

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

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern string			idb;
extern vector<Layer>	vL;
extern map<int,int>		mZ;
extern vector<Rgns>		vR;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetIDB( const char *tempdir );


