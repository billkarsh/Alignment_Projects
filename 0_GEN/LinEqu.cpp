

#include	"Disk.h"
#include	"File.h"
#include	"LinEqu.h"

#include	"Maths.h"
#include	"Memory.h"

#include	<string.h>
#include	<unistd.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	SWAP( a, b )											\
    {double temp = (a); (a) = (b); (b) = temp;}

/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	kMaxDirectN	256






/* --------------------------------------------------------------- */
/* Zero_Quick ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void Zero_Quick(
    double			*LHS,
    double			*RHS,
    int				n )
{
    memset( LHS, 0, n*n * sizeof(double) );
    memset( RHS, 0,   n * sizeof(double) );
}

/* --------------------------------------------------------------- */
/* AddConstraint_Quick ------------------------------------------- */
/* --------------------------------------------------------------- */

// Fastest way to build normal equations one constraint at a time.
//
// Assumes small LHS matrix in packed format:
// rows end-to-end in contiguous 1-D array.
//
void AddConstraint_Quick(
    double			*LHS,
    double			*RHS,
    int				n,
    int				nnz,
    const int		*j_nnz,
    const double	*Ai,
    double			Bi )
{
    for( int i = 0; i < nnz; ++i ) {

        int	ii = j_nnz[i],
            ni = n * ii;

        for( int j = 0; j < nnz; ++j )
            LHS[ni + j_nnz[j]] += Ai[i] * Ai[j];

        RHS[ii] += Ai[i] * Bi;
    }
}

/* --------------------------------------------------------------- */
/* AddConstraint_QuickWt ----------------------------------------- */
/* --------------------------------------------------------------- */

// Similar to AddConstraint_Quick(),
// but adds equation weighting factor.
//
void AddConstraint_QuickWt(
    double			*LHS,
    double			*RHS,
    int				n,
    int				nnz,
    const int		*j_nnz,
    const double	*Ai,
    double			Bi,
    double			Wt )
{
    Wt *= Wt;

    for( int i = 0; i < nnz; ++i ) {

        int	ii = j_nnz[i],
            ni = n * ii;

        for( int j = 0; j < nnz; ++j )
            LHS[ni + j_nnz[j]] += Wt * Ai[i] * Ai[j];

        RHS[ii] += Wt * Ai[i] * Bi;
    }
}

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

    C.push_back( LHSElem( val, row ) );
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
/* MATludcmp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* MATludcmp ---------------------------------------------------------
 *
 * LU decomposition.
 *
 * An adaptation of ludcmp from "Numerical Recipes in C."
 *
 * On entry:
 *
 * [A] is packed n x n matrix.
 * [indx] is dim-n workspace.
 * d is scalar workspace.
 * [vv] is dim-n workspace.
 *
 * On exit:
 *
 * [A] is replaced by row-permuted LU decomposition.
 * [indx] (dimension n) records perms.
 * d is {+1,-1} for {even,odd} # of perms.
 *
 * Return true if [A] nondegenerate.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

static int MATludcmp(
    double*		A,
    int*		indx,
    double*		d,
    double*		vv,
    int			n )
{
    const double	TINY = 1e-20;

    *d = 1.0;	// zero interchanges at start

// vv gets inv of biggest elem in each row

    for( int i = 0; i < n; ++i ) {

        double	big = 0.0;

        for( int j = 0; j < n; ++j ) {

            double	t;

            if( (t = fabs( A[n*i+j] )) > big )
                big = t;
        }

        if( !big )
            return false;

        vv[i] = 1.0 / big;
    }

// Crout's method: loop over cols

    for( int j = 0; j < n; ++j ) {

        for( int i = 0; i < j; ++i ) {

            double	sum = A[n*i+j];

            for( int k = 0; k < i; ++k )
                sum -= A[n*i+k] * A[n*k+j];

            A[n*i+j] = sum;
        }

        double	big = 0.0;
        int		imx;

        for( int i = j; i < n; ++i ) {

            double	sum = A[n*i+j];

            for( int k = 0; k < j; ++k )
                sum -= A[n*i+k] * A[n*k+j];

            A[n*i+j] = sum;

            double	t;

            if( (t = vv[i] * fabs( sum )) >= big ) {
                big = t;
                imx = i;
            }
        }

        // interchange rows?

        if( j != imx ) {

            for( int k = 0; k < n; ++k ) {

                double	t	= A[n*imx+k];
                A[n*imx+k]	= A[n*j+k];
                A[n*j+k]	= t;
            }

            *d		= -*d;
            vv[imx]	= vv[j];
        }

        indx[j] = imx;

        if( !A[n*j+j] )
            A[n*j+j] = TINY;

        if( j != n - 1 ) {

            double	t = 1.0 / A[n*j+j];

            for( int i = j + 1; i < n; ++i )
                A[n*i+j] *= t;
        }
    }

    return true;
}

/* --------------------------------------------------------------- */
/* MATlubksb ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* MATlubksb ---------------------------------------------------------
 *
 * Solve Ax=B using LU decomp from MATludcmp.
 *
 * An adaptation of lubksb from "Numerical Recipes in C."
 *
 * On entry:
 *
 * [B] is dim-n RHS col vector.
 * [A] is square n x n matrix from MATludcmp.
 * [indx] is dim-n result from MATludcmp.
 * d is scalar workspace.
 * [vv] is dim-n workspace.
 *
 * On exit:
 *
 * [B] is replaced by solution vector x.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

static void MATlubksb(
    double*			B,
    const double*	A,
    const int*		indx,
    int				n )
{
    int	ii = -1;

    for( int i = 0; i < n; ++i ) {

        int		ip  = indx[i];
        double	sum = B[ip];

        B[ip] = B[i];

        if( ii >= 0 ) {

            for( int j = ii; j < i; ++j )
                sum -= A[n*i+j] * B[j];
        }
        else if( sum )
            ii = i;

        B[i] = sum;
    }

    for( int i = n - 1; i >= 0; --i ) {

        double	sum = B[i];

        for( int j = i + 1; j < n; ++j )
            sum -= A[n*i+j] * B[j];

        B[i] = sum / A[n*i+i];
    }
}

/* --------------------------------------------------------------- */
/* SolveDirectLU ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SolveDirectLU(
    vector<double>			&X,
    const vector<LHSCol>	&LHS,
    const vector<double>	&RHS,
    int						n )
{
    double	A[n*n], vv[n], d;
    int		indx[n];

    memset( A, 0, n*n * sizeof(double) );

    for( int col = 0; col < n; ++col ) {

        const LHSCol&	C  = LHS[col];
        int				ne = C.size();

        for( int i = 0; i < ne; ++i )
            A[n*C[i].row+col] = C[i].val;
    }

    X = RHS;

    MATludcmp( A, indx, &d, vv, n );
    MATlubksb( &X[0], A, indx, n );
}

/* --------------------------------------------------------------- */
/* MATGaussJ ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* MATGaussJ ---------------------------------------------------------
 *
 * Solve linear equations [A].[X] = [B] by Gauss-Jordan elimination
 * with full pivoting.
 *
 * An adaptation of gaussj from "Numerical Recipes in C."
 *
 * Important differences from gaussj:
 *
 *	- Strictly zero-based array and matrix addressing.
 *	- Rearranged parameter list.
 *	- Return value indicates degeneracy.
 *	- Optimized pointer handling.
 *
 * On entry:
 *
 * [A] is square n x n matrix.
 * [B] is n x m, that is, m columns of length n.
 *
 * On exit:
 *
 * [A] is replaced by its inverse.
 * [B] is the set of solution column vectors.
 *
 * Return true if [A] nondegenerate.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

static int MATGaussJ(
    vector<vector<double> >	&A,
    vector<vector<double> >	&B,
    int						n,
    int						m )
{
    double		*p1, *p2;
    double		t, piv;
    vector<int>	indxc( n ),
                indxr( n ),
                ipiv( n, 0 );
    int			i, j, k, icol, irow, ok = true;

    for( i = 0; i < n; ++i ) {

        piv = 0.0;

        for( j = 0; j < n; ++j ) {

            if( ipiv[j] != 1 ) {

                p1 = &A[j][0];

                for( k = 0; k < n; ++k ) {

                    if( !ipiv[k] ) {

                        t = fabs( p1[k] );

                        if( t >= piv ) {
                            piv		= t;
                            irow	= j;
                            icol	= k;
                        }
                    }
                    else if( ipiv[k] > 1 ) {
                        ok = false;
                        goto exit;
                    }
                }
            }
        }

        ++ipiv[icol];

        /* ---------------------------------- */
        /* Swap rows to put pivot on diagonal */
        /* ---------------------------------- */

        if( irow != icol ) {

            p1 = &A[irow][0];
            p2 = &A[icol][0];

            for( j = 0; j < n; ++j )
                SWAP( p1[j], p2[j] );

            p1 = &B[irow][0];
            p2 = &B[icol][0];

            for( j = 0; j < m; ++j )
                SWAP( p1[j], p2[j] );
        }

        indxr[i]	= irow;
        indxc[i]	= icol;
        p1			= &A[icol][0];
        p2			= &B[icol][0];

        if( p1[icol] == 0.0 ) {
            ok = false;
            goto exit;
        }

        /* ------------------- */
        /* Normalize pivot row */
        /* ------------------- */

        piv			= 1.0 / p1[icol];
        p1[icol]	= 1.0;

        for( j = 0; j < n; ++j )
            p1[j] *= piv;

        for( j = 0; j < m; ++j )
            p2[j] *= piv;

        /* ------------ */
        /* Ellimination */
        /* ------------ */

        for( j = 0; j < n; ++j ) {

            if( j != icol ) {

                double	*p3	= &A[j][0];
                double	t	= p3[icol];

                p3[icol] = 0.0;

                for( k = 0; k < n; ++k )
                    p3[k] -= p1[k] * t;

                p3 = &B[j][0];

                for( k = 0; k < m; ++k )
                    p3[k] -= p2[k] * t;
            }
        }
    }

/* ------------------------------- */
/* Unswap columns in reverse order */
/* ------------------------------- */

    for( i = n - 1; i >= 0; --i ) {

        if( (irow = indxr[i]) != (icol = indxc[i]) ) {

            for( j = 0; j < n; ++j ) {

                double	*p1 = &A[j][0];

                SWAP( p1[irow], p1[icol] );
            }
        }
    }

exit:
    return ok;
}

/* --------------------------------------------------------------- */
/* SolveDirectGJ ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SolveDirectGJ(
    vector<double>			&X,
    const vector<LHSCol>	&LHS,
    const vector<double>	&RHS,
    int						n )
{
    vector<vector<double> >	a, b;

    a.resize( n );
    b.resize( n );

    for( int i = 0; i < n; ++i ) {
        a[i].resize( n, 0 );
        b[i].resize( 1, RHS[i] );
    }

    for( int col = 0; col < n; ++col ) {

        const LHSCol&	C  = LHS[col];
        int				ne = C.size();

        for( int i = 0; i < ne; ++i )
            a[C[i].row][col] = C[i].val;
    }

    MATGaussJ( a, b, n, 1 );

    for( int i = 0; i < n; ++i )
        X[i] = b[i][0];
}

/* --------------------------------------------------------------- */
/* Solve_Quick --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fastest in-place solver for small packed matrices.
//
// LHS and RHS are sent directly to NR routines, so are replaced
// by LU decomp and result vector directly. Caller should build
// RHS directly in output X array to minimize data copying.
//
bool Solve_Quick(
    double	*LHS,
    double	*RHS,
    int		n )
{
    double	vv[n], d;
    int		indx[n];

    if( MATludcmp( LHS, indx, &d, vv, n ) ) {
        MATlubksb( RHS, LHS, indx, n );
        return true;
    }

    return false;
}

/* --------------------------------------------------------------- */
/* WriteSolveRead ------------------------------------------------ */
/* --------------------------------------------------------------- */

// This is a general solver for a system of linear equations of
// the form A.X = B, where A is the LHS (NxN) matrix, and the B
// are the RHS (Nx1) constants. This is a blocking call.
//
// If N (# of unknowns) is <= kMaxDirectN this is solved by
// simple matrix inversion using routines in this source file.
//
// Otherwise the data are written to disk file 'triples', sent
// to external program 'SuperLUSymSolve' and the results are
// read back in from the disk file 'results'.
//
// If nproc == 1 we use single threaded SuperLUSymSolve; the
// system call will block until completion. Otherwise, we call
// SuperLUSymSolveMPI and wait for it to create semaphore file
// 'slu_signal'.
//
void WriteSolveRead(
    vector<double>			&X,
    const vector<LHSCol>	&LHS,
    const vector<double>	&RHS,
    const char				*jobtag,
    int						nproc,
    bool					uniqueNames )
{
    int	nvars = RHS.size();

/* ------ */
/* Direct */
/* ------ */

    if( nvars <= kMaxDirectN ) {
        SolveDirectLU( X, LHS, RHS, nvars );
        return;
    }

/* ----------------------------------- */
/* Print equations into 'triples' file */
/* ----------------------------------- */

// Name files

    char	iname[2048], oname[2048], buf[2048];

    if( uniqueNames ) {

        int	pid = getpid();

        gethostname( buf, sizeof(buf) );

        sprintf( iname, "triples_%s_%d", buf, pid );
        sprintf( oname, "results_%s_%d", buf, pid );
    }
    else {
        strcpy( iname, "triples" );
        strcpy( oname, "results" );
    }

// Delete any previous results file

    sprintf( buf, "rm -f %s", oname );
    system( buf );

// Open triples file

    FILE	*f	= FileOpenOrDie( iname, "w" );
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

    printf( "\n[[ Invoke solver: %s ]]\n", jobtag );
    VMStats( stdout );
    fflush( stdout );

    if( nproc == 1 ) {

        sprintf( buf, "SuperLUSymSolve -t -o=%s <%s", oname, iname );
        system( buf );
    }
    else {

        // remove signal file
        system( "rm -f slu_signal" );

        // submit mpi job
        sprintf( buf,
        "qsub -N %s -cwd -V -b y -pe impi %d"
        " 'mpirun -np %d"
        " SuperLUSymSolveMPI"
        " -r=1 -c=%d -t -j=%s -o=%s -i=%s'",
        jobtag, nproc, nproc, nproc, jobtag, oname, iname );
        system( buf );

        // await completion
        for(;;) {

            if( DskExists( "slu_signal" ) )
                break;

            sleep( 2 );
        }
    }

    VMStats( stdout );
    printf( "[[ Exit solver ]]\n\n" );
    fflush( stdout );

/* ------------ */
/* Read results */
/* ------------ */

    f = FileOpenOrDie( oname, "r" );

    for( int i = 0; i < nvars; ++i )
        fscanf( f, "%lf\n", &X[i] );

    fclose( f );
}


