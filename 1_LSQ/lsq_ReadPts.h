

#pragma once


#include	"lsq_CNX.h"
#include	"lsq_RGD.h"

#include	"GenDefs.h"

#include	<map>
using namespace std;


/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void ReadPts_StrTags(
	FILE		*FOUT,
	CNX			*cnx,
	RGD			*rgd,
	int			(*IDFromName)( const char *name ),
	const char	*dirfile,
	const char	*ptsfile );

void ReadPts_NumTags(
	FILE		*FOUT,
	CNX			*cnx,
	RGD			*rgd,
	const char	*ptsfile );

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern map<MZID,int>	nConRgn;	// # connected-rgn this image
extern int				gW, gH;		// image coords


