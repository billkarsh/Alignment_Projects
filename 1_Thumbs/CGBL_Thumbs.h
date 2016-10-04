

#pragma once


#include	"../1_DMesh/CThmUtil.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_Thumbs {

// =====
// Types
// =====

private:
	typedef struct {
		vector<TAffine>	Tdfm,
						Tab;
		double			SCALE,
						XSCALE,
						YSCALE,
						SKEW;
		const char		*ima,	// override idb paths
						*imb;
		int				FLD,
						MODE;
	} PrvDrvArgs;

public:
	typedef struct {
		double		CTR;
		const char	*fma,			// override idb paths
					*fmb;
		bool		Transpose,		// transpose all images
					SingleFold;		// assign id=1 to all non-fold rgns
	} DriverArgs;

	typedef struct {
		TAffine	Tdfm;
		double	XYCONF,
				NBMXHT,
				HFANGDN,
				HFANGPR,
				RTRSH;
		long	OLAP2D;
		int		FLD,
				MODE,
				THMDEC,
				OLAP1D,
				LIMXY;
	} CntxtDep;

// ============
// Data members
// ============

private:
	PrvDrvArgs	_arg;

public:
	DriverArgs	arg;
	TAffine		Tab;	// start thumbs here
	MatchParams	mch;
	CntxtDep	ctx;
	string		idb;
	PicSpec		A, B;

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


