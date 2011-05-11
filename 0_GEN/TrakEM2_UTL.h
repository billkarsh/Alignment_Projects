

#pragma once


#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

int  NextOID( const TiXmlHandle hDoc );

int  SetOID( TiXmlNode* par, int nextoid );

void CopyDTD( const char *dtdsrc, const char *bodysrc );


