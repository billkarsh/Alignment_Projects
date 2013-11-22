

#pragma once


/* --------------------------------------------------------------- */
/* class QueueEl ------------------------------------------------- */
/* --------------------------------------------------------------- */

class QueueEl {

public:
	int		x, y;	// we can get to (x,y)
	int		cost;	// for this cost

public:
	QueueEl( int x, int y, int cost )
	: x(x), y(y), cost(cost) {};

	bool operator < ( const QueueEl &rhs ) const
		{return cost > rhs.cost;};	// lowest first
};


