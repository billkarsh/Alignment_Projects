

#include	"lsq_Types.h"

#include	<stdlib.h>
#include	<string.h>


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
static const char		*noimg = "__noimg.jpg";
static Til2Img			_T1,
                        _T2;






/* --------------------------------------------------------------- */
/* RGN::RGN ------------------------------------------------------ */
/* --------------------------------------------------------------- */

RGN::RGN( const char *key )
{
    sscanf( key, "%d.%d-%d", &z, &id, &rgn );
    iname	= -1;
    itr		= -1;
}


RGN::RGN( const char *path, const DIR &dir, int _id )
{
    char	*s = (char*)strrchr( path, ':' );

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
    sscanf( key, "%d.%d-%d", &z, &id, &rgn );
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

        const Til2Img	*p;

        if( IDBT2ICacheNGet1( p, idb, z, id ) )
            return p->path.c_str();
        else
            return noimg;
    }
    else
        return rgnvname[iname].c_str();
}

/* --------------------------------------------------------------- */
/* RGN::GetMeta -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Pointer t1 required, t2 optional (or NULL).
//
void RGN::GetMeta(
    const Til2Img*	*t1,
    const Til2Img*	*t2,
    const RGN		&I1,
    const RGN		&I2 )
{
    if( I1.iname < 0 ) {

        if( t2 ) {

            IDBT2ICacheNGet2( *t1, *t2, idb,
                I1.z, I1.id, I2.z, I2.id );
        }
        else
            IDBT2ICacheNGet1( *t1, idb, I1.z, I1.id );
    }
    else {	// -strings option support, predates meta data

        char	buf[2048];

        strcpy( buf, rgnvname[I1.iname].c_str() );

        _T1.path	= strtok( buf, " ':\n" );
        _T1.id		= I1.id;
        _T1.col		= -999;
        _T1.row		= -999;
        _T1.cam		= 0;
        *t1			= &_T1;

        if( t2 ) {

            strcpy( buf, rgnvname[I2.iname].c_str() );

            _T2.path	= strtok( buf, " ':\n" );
            _T2.id		= I2.id;
            _T2.col		= -999;
            _T2.row		= -999;
            _T2.cam		= 0;
            *t2			= &_T2;
        }
    }
}


