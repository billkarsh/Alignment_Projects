

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

public:
	LHSElem( double val, int row )
		{this->val=val; this->row=row;};
};

typedef vector<LHSElem>	LHSCol;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void AddConstraint_Quick(
	double			*LHS,
	double			*RHS,
	int				n,
	int				nnz,
	const int		*j_nnz,
	const double	*Ai,
	double			Bi );

void AddConstraint(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nnz,
	const int		*j_nnz,
	const double	*Ai,
	double			Bi );

void Solve_Quick(
	double	*LHS,
	double	*RHS,
	int		n );

void WriteSolveRead(
	vector<double>			&X,
	const vector<LHSCol>	&LHS,
	const vector<double>	&RHS,
	const char				*jobtag,
	int						nproc,
	bool					uniqueNames );


