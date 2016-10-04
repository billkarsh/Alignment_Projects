

#pragma once


#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* class PixPair ------------------------------------------------- */
/* --------------------------------------------------------------- */

class PixPair {

private:
    // There are the actual pixel values for A and B.
    // There are always original full size (f) versions.
    // Scaled copies (s) are created if needed.
    // Filtered versions are created conditionally.
    vector<double>	_avf, _avs, _avfflt, _avsflt,
                    _bvf, _bvs, _bvfflt, _bvsflt;

public:
    // The public pointers are set according to usage...
    // either for the alignment phase or verify phase...
    // but they may both point to the vfy values if job
    // parameters do not specify filtering, and they may
    // point to (f) versions if scaling not required.
    const
    vector<double>	*avf_aln, *avs_aln, *avf_vfy, *avs_vfy,
                    *bvf_aln, *bvs_aln, *bvf_vfy, *bvs_vfy;
    vector<uint8>	resmska, resmskb;
    int				wf, hf,
                    ws, hs,
                    scl;

private:
    void Downsample(
        vector<double>			&dst,
        const vector<double>	&src );

public:
    bool Load(
        const PicSpec	&A,
        const PicSpec	&B,
        const string	&idb,
        bool			lens,
        bool			resmsk,
        int				order,
        int				bDoG,
        int				r1,
        int				r2,
        FILE*			flog = stdout,
        bool			transpose = false );
};


