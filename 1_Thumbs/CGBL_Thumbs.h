

#pragma once


#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_Thumbs {

// =====
// Types
// =====

public:
	typedef struct {
		double	SCALE,
				XSCALE,
				YSCALE,
				SKEW,
				CTR;
		char	*ima,			// override idb paths
				*imb,
				*fma,
				*fmb;
		bool	Transpose,		// transpose all images
				NoFolds,		// ignore fold masks
				SingleFold;		// assign id=1 to all non-fold rgns
	} DriverArgs;

	typedef struct {
		double	SCALE,
				XSCALE,
				YSCALE,
				SKEW,
				NBMXHT,
				HFANGDN,
				HFANGPR,
				RTRSH;
		long	OLAP2D;
		int		OLAP1D,
				INPALN,
				DINPUT;
	} CntxtDep;

	typedef struct {
		int			layer,
					tile;
		Til2Img		t2i;
		const char	*file;
	} PicSpecs;

// ============
// Data members
// ============

public:
	DriverArgs	arg;
	MatchParams	mch;
	CntxtDep	ctx;
	string		idb;
	PicSpecs	A, B;

// =================
// Object management
// =================

public:
	CGBL_Thumbs();

// =========
// Interface
// =========

public:
	bool SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern CGBL_Thumbs	GBL;


