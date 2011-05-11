

#pragma once


/* --------------------------------------------------------------- */
/* class QueueEl ------------------------------------------------- */
/* --------------------------------------------------------------- */

class QueueEl {

public:
	int		x, y;	// we can get to (x,y)
	int		cost;	// for this cost

public:
	QueueEl( int xx, int yy, int cc )
		{x = xx; y = yy; cost = cc;};

	bool operator < ( const QueueEl &rhs ) const
		{return cost > rhs.cost;};	// lowest first
};


