

#pragma once


#include	<set>
#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Layer {
public:
	int			z,
				sx, sy,
				dx, dy;
	set<int>	zdown;
public:
	Layer() : sx(-1), sy(-1), dx(-1), dy(-1) {};

	inline int Lowest() const
	{
		return (zdown.size() ? *zdown.begin() : z);
	};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void LayerCat(
	vector<Layer>	&vL,
	const char		*tempdir,
	const char		*cachedir,
	int				zolo,
	int				zohi,
	bool			catclr );


