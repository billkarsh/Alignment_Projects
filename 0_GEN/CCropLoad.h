

#pragma once


#include	"GenDefs.h"

#include	<stdio.h>

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CCropLoad {

private:
	typedef struct {
		int	x0, y0,
			dx, dy;
	} CBox;

private:
	CBox	B[4];
	FILE	*flog;
	int		isfile;

public:
	bool ReadIDB( const string &idb, FILE* flog = stdout );

	uint8* Raster8(
		const char*	name,
		int			cam,
		uint32		&w,
		uint32		&h,
		FILE*		flog = stdout,
		bool		transpose = false );
};


