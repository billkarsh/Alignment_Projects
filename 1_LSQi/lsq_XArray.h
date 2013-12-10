

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
	void Size( int ne );
	void Load( const char *path );
	void Updt();
	void Save();
public:
	static bool PriorIsAffine( const char *path );
};


