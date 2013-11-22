

#pragma once


#include	"CPoint.h"

#include	<stdio.h>

#include	<map>
#include	<set>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* class CorrPair ------------------------------------------------ */
/* --------------------------------------------------------------- */

class CorrPair {

	friend class CCorrImages;

private:
	int				i1, i2;	// image 1 and image 2
	vector<Point>	p1s;	// points in image 1
	vector<Point>	p2s;	// corresponding points in image 2

public:
	CorrPair( int i1, int i2 )
	: i1(i1), i2(i2) {};

	bool operator < ( const CorrPair &rhs ) const
		{return i1 < rhs.i1 || (i1 == rhs.i1 && i2 < rhs.i2);};
};

/* --------------------------------------------------------------- */
/* class CCorrImages --------------------------------------------- */
/* --------------------------------------------------------------- */

class CCorrImages {

private:
	map<string,int>	names;
	set<CorrPair>	corrs;

public:
	static CCorrImages *Read( const char *fname, FILE *flog );

	int Write( const char *fname );

	int Find(
		string			name1,
		string			name2,
		vector<Point>	&p1s,
		vector<Point>	&p2s );

	int Add(
		string			name1,
		string			name2,
		vector<Point>	&p1s,
		vector<Point>	&p2s );
};


