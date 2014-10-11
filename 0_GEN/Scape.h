

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

	ScpTile( const string &name, const TAffine &t2g )
	: name(name), t2g(t2g) {};
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
	bool			resmask,
	int				nthr,
	FILE*			flog );


