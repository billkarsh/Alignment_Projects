

#pragma once


#include	"CPicBase.h"


/* --------------------------------------------------------------- */
/* class Template ------------------------------------------------ */
/* --------------------------------------------------------------- */

class Template {

private:
    int			nx, ny;	// size of fft
    int			M;		// number of complex values in FFT
    int			x0, y0;	// origin of data; must add to all results
    vector<CD>	fft;	// FFT of the data

private:
    void FillFrWithNormValues(
        vector<double>	&fr,
        const PicBase	&P,
        int				i,
        int				xmin,
        int				xmax,
        int				ymin,
        int				ymax );

public:
    Template(
        const PicBase	&P,
        int				i,
        int				xmin,
        int				xmax,
        int				ymin,
        int				ymax );

    Point Match(
        const PicBase	&P,
        int				j,
        int				xmin,
        int				xmax,
        int				ymin,
        int				ymax );
};


