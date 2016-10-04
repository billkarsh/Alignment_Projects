

#pragma once


#include	"CPoint.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

double EarthMoversMetric(
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    bool					write_images,
    const char				*msg,
    FILE*					flog );

double FourierMatch(
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    int						wvlen,
    bool					write_images,
    const char				*msg,
    FILE*					flog );

double PercentYellow(
    const vector<double>	&a,
    const vector<double>	&b,
    FILE*					flog );


