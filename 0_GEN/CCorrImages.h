

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

class CorrImages;	// forward reference


class CorrPair {

	friend class CorrImages;
	friend CorrImages *ReadImageCorrelations(
				const char *fname, FILE *flog );

private:
	int				i1, i2;	// image 1 and image 2
	vector<Point>	p1s;	// points in image 1
	vector<Point>	p2s;	// corresponding points in image 2

public:
	CorrPair( int ii1, int ii2 )
		{i1 = ii1; i2 = ii2;};

	bool operator < ( const CorrPair &rhs ) const
		{return i1 < rhs.i1 || (i1 == rhs.i1 && i2 < rhs.i2);};
};

/* --------------------------------------------------------------- */
/* class CorrImages ---------------------------------------------- */
/* --------------------------------------------------------------- */

class CorrImages {

	friend CorrImages *ReadImageCorrelations(
				const char *fname, FILE *flog );

private:
	map<string,int>	names;
	set<CorrPair>	corrs;

public:
	int WriteImageCorrelations( const char *fname );
	int FindImageCorr(
		string			name1,
		string			name2,
		vector<Point>	&p1s,
		vector<Point>	&p2s );
	int AddImageCorr(
		string			name1,
		string			name2,
		vector<Point>	&p1s,
		vector<Point>	&p2s );
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

CorrImages *ReadImageCorrelations( const char *fname, FILE *flog );


