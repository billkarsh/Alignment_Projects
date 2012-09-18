

#pragma once


#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_dmesh {

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
		bool	ForceSkew,		// apply SKEW args even if same layer
				Transpose,		// transpose all images
				WithinSection,	// overlap within a section
				Verbose,		// run inspect diagnostics
				NoFolds,		// ignore fold masks
				SingleFold,		// assign id=1 to all non-fold rgns
				Heatmap;		// run CorrView
	} DriverArgs;

	typedef struct {
		double	SCALE,
				XSCALE,
				YSCALE,
				SKEW,
				NBMXHT,
				HFANGDN,
				HFANGPR,
				RTRSH,
				DIT,
				DFA,
				DFT;
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
	DriverArgs		arg;
	vector<TForm>	Tab;	// start thumbs here
	vector<TForm>	Tusr;	// bypass thumbs, start mesh here
	vector<Point>	XYusr;	// command line expected XY
	MatchParams		mch;
	CntxtDep		ctx;
	string			idb;
	PicSpecs		A, B;

// =================
// Object management
// =================

public:
	CGBL_dmesh();

// =========
// Interface
// =========

public:
	bool SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern CGBL_dmesh	GBL;


