

#pragma once


#include	"FoldMask.h"
#include	"CPixPair.h"
#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool ApproximateMatch_NoCR(
    vector<TAffine>	&guesses,
    const PixPair	&px,
    FILE*			flog );

bool ApproximateMatch(
    vector<TAffine>		&guesses,
    const PixPair		&px,
    const ConnRegion	&acr,
    const ConnRegion	&bcr,
    FILE*				flog );


