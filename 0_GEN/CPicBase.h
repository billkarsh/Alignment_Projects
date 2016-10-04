

#pragma once


#include	"GenDefs.h"
#include	"TAffine.h"

#include	<stdio.h>

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* class PicBase ------------------------------------------------- */
/* --------------------------------------------------------------- */

class PicBase {

public:
    string	fname;				// file name
    uint8*	raster;				// used most often (usually 2K x 2K)
    uint8*	original;			// original regular size raster
    uint8*	external;			// an original we didn't load
    uint32	w, h;				// image dims
    int		z;					// Z layer
    int		scale;				// original/raster (1, 2, 4)
    TAffine	tr;					// transform this raster to target
    TAffine	Inverse;			// tr inverse
    vector<CD>		fft_of_frame;	// FFT of the frame
    vector<uint8>	DoG;			// Difference of Gaussians

public:
    PicBase();
    virtual ~PicBase();

    void LoadOriginal( const char* name, FILE* flog, bool transpose );
    void SetExternal( const uint8* in_raster, uint32 w, uint32 h );
    void CopyOriginal();
    void DownsampleIfNeeded( FILE* flog );
    void MakeFFTExist( int i );
    void MakeDoGExist( vector<CD> &filter, int r1, int r2 );
};


