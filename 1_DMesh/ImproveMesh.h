

#pragma once


#include	"Geometry.h"
#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

double ImproveMesh(
    vector<TAffine>			&transforms,
    vector<Point>			&centers,
    vector<triangle>		&tri,
    const vector<vertex>	&ctl,
    const vector<Point>		&apts,
    const vector<double>	&av,
    const vector<double>	&bimg,
    int						w,
    int						h,
    const TAffine			&tr_guess,
    double					threshold,
    FILE					*flog,
    const char				*describe );


