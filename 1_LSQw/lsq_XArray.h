

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
	void Resize( int ne );
	void Load( const char *path );
	void Save();
	void UpdtFS();
private:
	bool Send( int zlo, int zhi, int XorF, int toLorR );
	bool Recv( int zlo, int zhi, int XorF, int fmLorR );
public:
	bool Updt();
};


