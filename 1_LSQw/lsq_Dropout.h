

#pragma once


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Dropout {
public:
    long	rmax, read, pnts, kill, cutd;
private:
    void Add( const Dropout& rhs )
    {
        rmax += rhs.rmax;
        read += rhs.read;
        pnts += rhs.pnts;
        kill += rhs.kill;
        cutd += rhs.cutd;
    };
    void GatherCounts();
public:
    Dropout() : rmax(0), read(0), pnts(0), kill(0), cutd(0) {};
    void Scan();
};


