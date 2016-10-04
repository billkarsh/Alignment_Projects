

#pragma once


#include	"Geometry.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MeshGetBounds(
    IBox				&B,
    const vector<Point>	&pts,
    FILE*				flog );

void MeshMakeSingleTri(
    vector<triangle>	&tri,
    vector<vertex>		&ctl,
    const IBox			&B,
    FILE*				flog );

int MeshCreate(
    vector<triangle>	&tri,
    vector<vertex>		&ctl,
    const vector<Point>	&pts,
    const IBox			&B,
    FILE*				flog );

int MeshCreateX(
    vector<triangle>	&tri,
    vector<vertex>		&ctl,
    const vector<Point>	&pts,
    const IBox			&B,
    FILE*				flog );


