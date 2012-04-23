

#pragma once


#include	"GenDefs.h"
#include	"CPixPair.h"
#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void ABOverlay(
	const PixPair	&px,
	const uint16*	rmap,
	int				Ntrans,
	const TForm*	tfs,
	const TForm*	ifs );

void RunCorrView(
	const PixPair	&px,
	const uint16*	rmap,
	const TForm*	tfs,
	bool			heatmap );


