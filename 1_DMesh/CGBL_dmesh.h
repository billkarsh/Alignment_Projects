

#pragma once


#include	"CThmUtil.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_dmesh {

// =====
// Types
// =====

private:
    typedef struct {
        vector<TAffine>	Tdfm,
                        Tab,
                        Ta,
                        Tb;
        double			SCALE,
                        XSCALE,
                        YSCALE,
                        SKEW;
        const char		*matchparams,
                        *ima,			// override idb paths
                        *imb;
        int				FLD,
                        MODE;
    } PrvDrvArgs;

public:
    typedef struct {
        double		CTR;
        const char	*fma,				// override idb paths
                    *fmb,
                    *comp_png,			// override comp.png path
                    *registered_png;	// override registered.png path
        bool		Transpose,			// transpose all images
                    WithinSection,		// overlap within a section
                    SingleFold,			// assign id=1 to all non-fold rgns
                    JSON,				// output JSON format
                    Verbose,			// run inspect diagnostics
                    Heatmap;			// run CorrView
    } DriverArgs;

    typedef struct {
        TAffine	Tdfm;
        double	XYCONF,
                NBMXHT,
                HFANGDN,
                HFANGPR,
                RTRSH,
                RIT,
                RFA,
                RFT;
        long	OLAP2D;
        int		FLD,
                MODE,
                THMDEC,
                OLAP1D,
                LIMXY,
                OPT;
    } CntxtDep;

// ============
// Data members
// ============

private:
    PrvDrvArgs		_arg;

public:
    DriverArgs		arg;
    TAffine			Tab;	// start thumbs here
    vector<TAffine>	Tmsh;	// bypass thumbs, start mesh here
    vector<Point>	XYexp;	// command line expected XY
    MatchParams		mch;
    CntxtDep		ctx;
    string			idb;
    PicSpec			A, B;

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


