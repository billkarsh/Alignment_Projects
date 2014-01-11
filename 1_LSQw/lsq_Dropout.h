

#pragma once


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Dropout {
public:
	long	pnts, kill, cutd;
private:
	void Add( const Dropout& rhs )
		{pnts += rhs.pnts; kill += rhs.kill; cutd += rhs.cutd;};
	void GatherCounts();
public:
	Dropout() : pnts(0), kill(0), cutd(0) {};
	void Scan();
};


