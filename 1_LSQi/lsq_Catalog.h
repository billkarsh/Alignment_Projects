

#pragma once


#include	<set>
#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSubdirCat {
public:
	int			z,
				sx, sy,
				dx, dy;
	set<int>	zdown;
public:
	CSubdirCat()	{sx=-1; sy=-1; dx=-1; dy=-1;};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CatPoints(
	vector<CSubdirCat>	&vC,
	const char			*tempdir,
	int					zolo,
	int					zohi,
	bool				clrcat );


