

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
		bool	Transpose,		// transpose all images
				WithinSection,	// overlap within a section
				Verbose,		// run inspect diagnostics
				NoFolds,		// ignore fold masks
				SingleFold,		// assign id=1 to all non-fold rgns
				Heatmap;		// run CorrView
	} DriverArgs;

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
	vector<TForm>	Tusr;	// command line guessed tform
	vector<Point>	XYusr;	// command line expected XY
	ThmParams		thm;
	MeshParams		msh;
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


