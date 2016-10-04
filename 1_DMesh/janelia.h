

#pragma once


#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#ifdef USE_CURL
bool GetTileSpecFromURL( PicSpec &P, const char *pat, char *argv );
#else
bool GetTileSpecFromURL( PicSpec &P, const char *pat, char *argv )
    {return false;}
#endif


