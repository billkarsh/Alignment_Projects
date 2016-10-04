

#pragma once


#include	"GenDefs.h"
#include	"CPixPair.h"
#include	"TAffine.h"
#include	"THmgphy.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void YellowView(
    const PixPair	&px,
    const TAffine	&T,
    FILE			*flog );

void YellowView(
    const PixPair	&px,
    const THmgphy	&T,
    FILE			*flog );

void ABOverlay(
    const PixPair	&px,
    const uint16*	rmap,
    int				Ntrans,
    const TAffine*	tfs,
    const TAffine*	ifs,
    const char		*comp_png = NULL,
    FILE			*flog = stdout );

void RunCorrView(
    const PixPair	&px,
    const uint16*	rmap,
    const TAffine*	tfs,
    bool			heatmap,
    const char		*registered_png = NULL,
    FILE			*flog = stdout );


