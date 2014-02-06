

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
	vector<vector<int> >	K;	// color
	int						splitmin;
public:
	Split( const XArray& X, int splitmin )
	: X(X), splitmin(splitmin) {};
private:
	void Resize();
	void ColorMontages();
	bool Send( int zlo, int zhi, int toLorR );
	bool Recv( int zlo, int zhi, int fmLorR );
	bool Updt();
	void ReportCount( const map<int,int>& m );
	void CountColors( map<int,int>& m );
	void Save();
	void BreakOut( const map<int,int>& m );
public:
	void Run();
};


