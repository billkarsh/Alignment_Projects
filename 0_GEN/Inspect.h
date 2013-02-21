

#pragma once


#include	"GenDefs.h"
#include	"CPixPair.h"
#include	"TAffine.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void ABOverlay(
	const PixPair	&px,
	const uint16*	rmap,
	int				Ntrans,
	const TAffine*	tfs,
	const TAffine*	ifs );

void RunCorrView(
	const PixPair	&px,
	const uint16*	rmap,
	const TAffine*	tfs,
	bool			heatmap );


