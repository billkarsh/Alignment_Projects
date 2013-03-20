

#pragma once


#include	"GenDefs.h"
#include	"TAffine.h"

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class ScpTile {

public:
	string	name;
	TAffine	t2g;	// tile to global

public:
	ScpTile()	{};

	ScpTile( const string &_name, const TAffine &_t2g )
		{name = _name; t2g = _t2g;};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

uint8* Scape(
	uint32			&ws,
	uint32			&hs,
	double			&x0,
	double			&y0,
	vector<ScpTile>	&vTile,
	int				wi,
	int				hi,
	double			scale,
	int				szmult,
	int				bkval,
	int				lgord,
	int				sdnorm,
	FILE*			flog );


