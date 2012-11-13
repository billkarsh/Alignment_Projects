

#pragma once


#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_dmesh {

// =====
// Types
// =====

private:
	typedef struct {
		vector<TForm>	Tdfm;
		double			SCALE,
						XSCALE,
						YSCALE,
						SKEW;
		char			*ima,	// override idb paths
						*imb;
	} PrvDrvArgs;

public:
	typedef struct {
		double	CTR;
		char	*fma,			// override idb paths
				*fmb;
		bool	Transpose,		// transpose all images
				WithinSection,	// overlap within a section
				Verbose,		// run inspect diagnostics
				NoFolds,		// ignore fold masks
				SingleFold,		// assign id=1 to all non-fold rgns
				Heatmap;		// run CorrView
	} DriverArgs;

	typedef struct {
		TForm	Tdfm;
		double	NBMXHT,
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

private:
	PrvDrvArgs		_arg;

public:
	DriverArgs		arg;
	TForm			Tab;	// start thumbs here
	vector<TForm>	Tmsh;	// bypass thumbs, start mesh here
	vector<Point>	XYexp;	// command line expected XY
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


