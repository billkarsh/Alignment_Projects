

#pragma once


#include	"tinyxml.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

int  NextOID( const TiXmlHandle hDoc );

int  SetOID( TiXmlNode* par, int nextoid );

void CopyDTD( const char *dtdsrc, const char *bodysrc );

void TrakEM2WriteDTD( FILE *f );

void TrakEM2WriteDTDEx( FILE *f );


