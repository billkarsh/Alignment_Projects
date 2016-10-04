

#pragma once


#include	"FoldMask.h"
#include	"CPixPair.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Thumbs_NoCR( const PixPair &px, FILE* flog );

bool Thumbs(
    const PixPair		&px,
    const ConnRegion	&acr,
    const ConnRegion	&bcr,
    FILE*				flog );


