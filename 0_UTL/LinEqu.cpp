

#include	"File.h"
#include	"LinEqu.h"

#include	"Maths.h"






/* --------------------------------------------------------------- */
/* AddToElem ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Accumulate summed terms into the {row,col} element of LHS.
// New rows are added to LHS as needed.
//
static void AddToElem(
	vector<LHSCol>	&LHS,
	int				row,
	int				col,
	double			val )
{
	LHSCol&	C  = LHS[col];
	int		ne = C.size();

	for( int i = 0; i < ne; ++i ) {

		if( C[i].row == row ) {
			C[i].val += val;
			return;
		}
	}

// Not found - Add

	LHSElem	e;

	e.val = val;
	e.row = row;

	C.push_back( e );
}

/* --------------------------------------------------------------- */
/* Print4x4 ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void Print4x4( double a[4][4] )
{
	printf( "[ " );

	for( int row = 0; row < 4; ++row ) {

		for( int col = 0; col < 4; ++col )
			printf( "%f ", a[row][col] );

		if( row != 3 )
			printf( ";" );
	}

	printf( " ]\n" );
}

/* --------------------------------------------------------------- */
/* SolveExplicit4x4 ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void SolveExplicit4x4(
	vector<double>			&X,
	const vector<LHSCol>	&LHS,
	const vector<double>	&RHS )
{
	double a[4][4];
	double b[4][4];

// Convert from sparse

	memset( &a[0][0], 0, 16 * sizeof(double) );

	for( int col = 0; col < 4; ++col ) {

		const LHSCol&	C  = LHS[col];
		int				ne = C.size();

		for( int i = 0; i < ne; ++i )
			a[C[i].row][col] = C[i].val;
	}

// Solve

	Invert4x4Matrix( b, a );

	for( int i = 0; i < 4; ++i ) {

		double	sum = 0.0;

		for( int j = 0; j < 4; ++j )
			sum += b[i][j] * RHS[j];

		X[i] = sum;
	}
}

/* --------------------------------------------------------------- */
/* AddConstraint ------------------------------------------------- */
/* --------------------------------------------------------------- */

// This is a tool for building up a set of normal equations,
// one constraint (equation) at a time.
//
// Often one has an overdetermined system wherein there are
// N unknowns and M > N equations for them. If the unknowns
// are represented by (Nx1) column vector X, we have A.X = B,
// with A (MxN) and B (Mx1). We can apply a linear solver by
// multiplying both sides on the left by A-transpose (NxM),
// which produces the N "normal equations" for the system.
//
// Let the original (MxN) system be represented as:
//
//	Aij . Xj = Bi,  i=[1,,M] eqns, j=[1,,N] prms.
//
// Using this notation, call this function for each (i) of
// the M equations, supplying information for the ith row:
//
// - nnz:	The number of non-zero columns (j) on this row of Aij.
// - j_nnz:	The indices (j) of those non-zero elements.
// - Ai:	The values of those non-zero elements.
// - Bi:	The rhs constant term for this equation.
//
// This function computes the current contributions to the LHS
// and RHS sums for the transpose multiplication.
//
// To help see how the matrix mechanics really works out right,
// note that the rows of the original matrix Aij are equations,
// so the columns of Aki-transpose are the same equations making
// the product Aki.Aij look like:
//
//    |  .   .   .   .      |   | ..eq1.. |       | b1 |
//    |  .   .   .   .      |   | ..eq2.. |       | b2 |
//    | eq1 eq2 eq3 eq4 ... | x | ..eq3.. | ; B = | b3 |
//    |  .   .   .   .      |   | ..eq4.. |       | b4 |
//    |  .   .   .   .      |   | ....... |       | .. |
//
// Here one can see that the terms of eq-i only multiply other
// terms of eq-i and no others, so the LHS normal matrix elements
// can be computed incrementally, one equation at a time. One can
// also see this for the RHS, where a given b-i only multiplies
// terms from eq-i.
//
void AddConstraint(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	int				nnz,
	const int		*j_nnz,
	const double	*Ai,
	double			Bi )
{
	for( int i = 0; i < nnz; ++i ) {

		int	ii = j_nnz[i];

		for( int j = 0; j < nnz; ++j ) {

			int	jj = j_nnz[j];

			AddToElem( LHS, ii, jj, Ai[i] * Ai[j] );
		}

		RHS[ii] += Ai[i] * Bi;
	}
}

/* --------------------------------------------------------------- */
/* WriteSolveRead ------------------------------------------------ */
/* --------------------------------------------------------------- */

// This is a general solver for a system of linear equations of
// the form A.X = B, where A is the LHS (NxN) matrix, and the B
// are the RHS (Nx1) constants.
//
// If N (# of unknowns) is 4 this is solved by simple matrix
// inversion using routines in this source file.
//
// Otherwise the data are written to disk file 'triples', sent
// to external program 'SuperLUSymSolve' and the results are
// read back in from the disk file 'results'.
//
void WriteSolveRead(
	vector<double>			&X,
	const vector<LHSCol>	&LHS,
	const vector<double>	&RHS,
	bool					uniqueNames )
{
	int	nvars = RHS.size();

/* --------------------------------- */
/* Handle at least this special case */
/* --------------------------------- */

	if( nvars == 4 ) {
		SolveExplicit4x4( X, LHS, RHS );
		return;
	}

/* ----------------------------------- */
/* Print equations into 'triples' file */
/* ----------------------------------- */

// Name files

	char	tname[2048], rname[2048], buf[2048];

	if( uniqueNames ) {

		int	pid = getpid();

		gethostname( buf, sizeof(buf) );

		sprintf( tname, "/tmp/triples_%s_%d", buf, pid );
		sprintf( rname, "/tmp/results_%s_%d", buf, pid );
	}
	else {
		strcpy( tname, "triples" );
		strcpy( rname, "results" );
	}

// Open file

	FILE	*f	= FileOpenOrDie( tname, "w" );
	int		nnz	= 0;	// number of non-zero terms

// Header

	for( int col = 0; col < nvars; ++col )
		nnz += LHS[col].size();

	fprintf( f, "%d %d\n", nvars, nnz );

// LHS

	for( int col = 0; col < nvars; ++col ) {

		const LHSCol&	C  = LHS[col];
		int				ne = C.size();

		for( int i = 0; i < ne; ++i ) {

			int	row = C[i].row;

			// convert to 1 based indexing

			fprintf( f, "%d %d %.16f\n",
			row + 1, col + 1, C[i].val );
		}
	}

// RHS

	for( int i = 0; i < nvars; ++i )
		fprintf( f, "%.16f\n", RHS[i] );

// Close

	fclose( f );

/* ----- */
/* Solve */
/* ----- */

	printf( "\n[[ Invoke solver ]]\n" );
	fflush( stdout );

	sprintf( buf, "SuperLUSymSolve -t -o=%s <%s", rname, tname );
	system( buf );

	fflush( stdout );
	printf( "[[ Exit solver ]]\n\n" );

/* ------------ */
/* Read results */
/* ------------ */

	f = FileOpenOrDie( rname, "r" );

	for( int i = 0; i < nvars; ++i )
		fscanf( f, "%lf", &X[i] );

	fclose( f );
}


