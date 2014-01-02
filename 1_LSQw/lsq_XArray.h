

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
	static bool PriorIsAffine( const char *path );
public:
	void Size( int ne );
	void Load( const char *path );
	void UpdtFS();
	void Save();
private:
	bool Send( int zlo, int zhi, int XorU, int toLorR );
	bool Recv( int zlo, int zhi, int XorU, int fmLorR );
public:
	bool Updt();
};


