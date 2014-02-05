

#pragma once


#include	"lsq_XArray.h"

#include	<map>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Split {
public:
	const XArray&			X;
	vector<vector<int> >	C;
	int						splitmin;
public:
	Split( const XArray& X, int splitmin )
	: X(X), splitmin(splitmin) {};
private:
	void Resize();
	void ReportCount( const map<int,int>& m );
	void CountColors( map<int,int>& m );
	void Save();
	void BreakOut( const map<int,int>& m );
public:
	void Run();
};


