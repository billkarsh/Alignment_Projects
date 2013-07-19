

#pragma once


#include	"tinyxml.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* class XML_TKEM ------------------------------------------------ */
/* --------------------------------------------------------------- */

class XML_TKEM {
public:
	FILE*			flog;
	const char		*file;
	TiXmlDocument	doc;
public:
	XML_TKEM( const char *file, FILE* flog=stdout )
		{Open( file, flog );};
	void Open( const char *file, FILE* flog=stdout );
	void Save( const char *name, bool copyDTD );
	TiXmlNode*		GetLayerset();
	TiXmlElement*	GetFirstLayer();
	TiXmlElement*	GetLastLayer();
	int				NextOID();
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

int IDFromPatch( TiXmlElement* ptch );

void XMLSetTFVals( TiXmlElement* ptch, const double *t );

int  NextOID( const TiXmlHandle hDoc );

int  SetOID( TiXmlNode* par, int nextoid );

void CopyDTD( const char *dtdsrc, const char *bodysrc );

void TrakEM2WriteDTD( FILE *f );

void TrakEM2WriteDTDEx( FILE *f );


