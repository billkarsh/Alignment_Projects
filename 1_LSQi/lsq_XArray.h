

#pragma once


#include	<vector>
using namespace std;

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class XArray {
public:
	int						NE;
	vector<vector<double> >	X;
public:
	void Load( const char *path );
};


