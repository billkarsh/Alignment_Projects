

#pragma once


#include	"GenDefs.h"
#include	"CPixPair.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void PipelineDeformableMap(
    int				&Ntrans,
    double*			&tr_array,
    uint16*			map_mask,
    const PixPair	&px,
    const uint8*	fold_mask_a,
    const uint8*	fold_mask_b,
    FILE*			flog );


