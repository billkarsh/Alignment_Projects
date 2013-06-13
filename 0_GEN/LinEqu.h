

#pragma once


#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// A-matrix element in system A.X = B

class LHSElem {

public:
	double	val;
	int		row;
};

typedef vector<LHSElem>	LHSCol;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void AddConstraint(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nnz,
	const int		*j_nnz,
	const double	*Ai,
	double			Bi );

void WriteSolveRead(
	vector<double>			&X,
	const vector<LHSCol>	&LHS,
	const vector<double>	&RHS,
	const char				*jobtag,
	int						nproc,
	bool					uniqueNames );


