

#pragma once


#include	"CPoint.h"


/* --------------------------------------------------------------- */
/* SML ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Similarity alignment testing

class SML {

private:
    typedef struct {
        // rgn-pair's pts
        vector<Point>	pa;
        vector<Point>	pb;
    } CorrPts;

private:
    vector<CorrPts>	r12Pts;

private:
    double CanAlign(
        const vector<Point>	&p1,
        const vector<Point>	&p2,
        bool				print );

public:
    void AddPOINTPair(
        int			r1,
        const Point	&p1,
        int			r2,
        const Point	&p2 );

    void TestPairAlignments();
};


