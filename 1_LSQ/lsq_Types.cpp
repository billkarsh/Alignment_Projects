

#include	"lsq_Types.h"

#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

string				idb;		// for name lookups
vector<RGN>			vRgn;		// the regions
map<CRPair,int>		r12Idx;		// idx from region-pair
vector<Constraint>	vAllC;		// all Point-pairs

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static vector<string>	rgnvname;	// tile names
static map<MZID,int>	rgnmname;	// tile names






/* --------------------------------------------------------------- */
/* RGN::RGN ------------------------------------------------------ */
/* --------------------------------------------------------------- */

RGN::RGN( const char *key )
{
	sscanf( key, "%d.%d:%d", &z, &id, &rgn );
	iname	= -1;
	itr		= -1;
}


RGN::RGN( const char *path, const DIR &dir, int _id )
{
	char	*s = strrchr( path, ':' );

	s[-1]	= 0;
	z		= dir.ZFromName( path );
	id		= _id;
	rgn		= atoi( s + 1 );
	itr		= -1;

	MZID					zid( z, id );
	map<MZID,int>::iterator	it = rgnmname.find( zid );

	if( it == rgnmname.end() ) {

		rgnmname[zid] = iname = rgnvname.size();
		rgnvname.push_back( path );
	}
	else
		iname = it->second;
}


RGN::RGN( const char *path, const char *key )
{
	sscanf( key, "%d.%d:%d", &z, &id, &rgn );
	itr		= -1;

	MZID					zid( z, id );
	map<MZID,int>::iterator	it = rgnmname.find( zid );

	if( it == rgnmname.end() ) {

		rgnmname[zid] = iname = rgnvname.size();
		rgnvname.push_back( path );
	}
	else
		iname = it->second;
}

/* --------------------------------------------------------------- */
/* RGN::GetName -------------------------------------------------- */
/* --------------------------------------------------------------- */

const char* RGN::GetName() const
{
	if( iname < 0 ) {

		Til2Img	I;

		if( !IDBTil2Img( I, idb, z, id ) )
			I.path = "__noimg.jpg";

		return I.path.c_str();
	}
	else
		return rgnvname[iname].c_str();
}


