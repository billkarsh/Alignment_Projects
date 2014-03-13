

#pragma once


#include	"FoldMask.h"
#include	"CPixPair.h"
#include	"Cffmap.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void RegionToRegionMap(
	ffmap				&maps,
	uint16				*ids,
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	TAffine				tr_guess,
	FILE				*flog,
	FILE				*ftri );


