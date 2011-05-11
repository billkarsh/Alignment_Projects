

#pragma once


#include	"Geometry.h"
#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

double ImproveMesh(
	vector<TForm>			&transforms,
	vector<Point>			&centers,
	vector<triangle>		&tri,
	const vector<vertex>	&ctl,
	const vector<Point>		&apts,
	const vector<double>	&av,
	const vector<double>	&bimg,
	int						w,
	int						h,
	const TForm				&tr_guess,
	double					threshold,
	FILE					*flog,
	const char				*describe );


