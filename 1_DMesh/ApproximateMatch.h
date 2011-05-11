

#pragma once


#include	"FoldMask.h"

#include	"CPixPair.h"
#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool ApproximateMatch_NoCR(
	vector<TForm>	&guesses,
	const PixPair	&px,
	FILE*			flog );

bool ApproximateMatch(
	vector<TForm>		&guesses,
	const PixPair		&px,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog );


