

#pragma once


#include	<complex>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const int		BIG	= 0x7FFFFFFF;	// biggest 32-bit int
const double	PI	= 3.14159265358;

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
typedef unsigned char	uint8;
typedef unsigned short	uint16;
typedef unsigned int	uint32;


