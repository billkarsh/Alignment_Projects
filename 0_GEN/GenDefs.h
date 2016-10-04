

#pragma once


#include	<complex>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const int		BIG		= 0x7FFFFFFF;	// biggest 32-bit int
const double	BIGD	= 1.0e30;
const double	PI		= 3.14159265358;

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef complex<double> CD;

typedef struct {
    int		L,
            R,
            B,
            T;
} IBox;

typedef struct {
    double	L,
            R,
            B,
            T;
} DBox;

// tiffio.h types
typedef signed char		int8;
typedef signed short	int16;
typedef signed int		int32;
typedef unsigned char	uint8;
typedef unsigned short	uint16;
typedef unsigned int	uint32;


class MZID {
// Use for mappings: map<MZID,xxx>
public:
    int	z, id;
public:
    MZID()	{};
    MZID( int z, int id ) : z(z), id(id) {};

    bool operator < ( const MZID &rhs ) const
        {
            if( z < rhs.z )
                return true;
            if( z > rhs.z )
                return false;

            return id < rhs.id;
        };
    bool operator == ( const MZID &rhs ) const
        {return z == rhs.z && id == rhs.id;};
};


class MZIDR {
// Use for mappings: map<MZIDR,xxx>
public:
    int	z, id, rgn;
public:
    MZIDR()	{};
    MZIDR( int z, int id, int rgn ) : z(z), id(id), rgn(rgn) {};

    bool operator < ( const MZIDR &rhs ) const
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
    bool operator == ( const MZIDR &rhs ) const
        {return z == rhs.z && id == rhs.id && rgn == rhs.rgn;};
};


