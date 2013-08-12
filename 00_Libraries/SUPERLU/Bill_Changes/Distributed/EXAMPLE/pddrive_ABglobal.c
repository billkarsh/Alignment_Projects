

/*! @file 
 * \brief Driver program for pdgssvx_ABglobal example
 *
 * <pre>
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 * </pre>
 */

#include <math.h>
#include "superlu_ddefs.h"

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * The driver program pddrive_ABglobal.
 *
 * This example illustrates how to use pdgssvx_ABglobal with the full
 * (default) options to solve a linear system.
 * 
 * Five basic steps are required:
 *   1. Initialize the MPI environment and the SuperLU process grid
 *   2. Set up the input matrix and the right-hand side
 *   3. Set the options argument
 *   4. Call pdgssvx_ABglobal
 *   5. Release the process grid and terminate the MPI environment
 *
 * On an IBM SP, the program may be run by typing
 *    poe pddrive_ABglobal -r <proc rows> -c <proc columns> <input_file> -procs <p>
 * </pre>
 */

int main(int argc, char *argv[])
{
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    gridinfo_t grid;
    double   *berr;
    double   *a, *b, *xtrue;
    int_t    *asub, *xa;
    int_t    m, n, nnz;
    int_t    nprow, npcol;
    int      iam, info, ldb, ldx, nrhs;
    char     trans[1];
    char     *jobtag, *outfile, *infile, ftriple = 0;
    extern int cpp_defs();


    nprow = 1;  /* Default process rows.      */
    npcol = 1;  /* Default process columns.   */
    nrhs = 1;   /* Number of right-hand side. */

    /* ------------------------------------------------------------
       INITIALIZE MPI ENVIRONMENT. 
       ------------------------------------------------------------*/
    MPI_Init( &argc, &argv );

    /* Parse command line argv[]. */

	int	i;
	for( i = 1; i < argc; ++i ) {
	
		if( argv[i][1] == 'r' ) {

			nprow = atoi( argv[i] + 3 );
		}
		else if( argv[i][1] == 'c' ) {

			npcol = atoi( argv[i] + 3 );
		}
		else if( argv[i][1] == 't' ) {

			ftriple = 1;
		}
		else if( argv[i][1] == 'j' ) {

			jobtag = argv[i] + 3;
		}
		else if( argv[i][1] == 'o' ) {

			outfile = argv[i] + 3;
		}
		else if( argv[i][1] == 'i' ) {

			infile = argv[i] + 3;
		}
	}

    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID. 
       ------------------------------------------------------------*/
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
    iam = grid.iam;

    /* Bail out if I do not belong in the grid. */
    if ( iam >= nprow * npcol )
	goto out;

	char	slog[128];
	sprintf( slog, "slu_%s_%d.txt", jobtag, iam );
	freopen( slog, "a", stdout );
	printf( "\n\n----- New Job -----\n" );
	printf( "in/out/t/r/c = %s/%s/%d/%d/%d\n", infile, outfile, ftriple, nprow, npcol );
	VMStats( stdout );

#if ( VAMPIR>=1 )
    VT_traceoff();
#endif

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter main()");
#endif

    /* ------------------------------------------------------------
       PROCESS 0 READS THE MATRIX A, AND THEN BROADCASTS IT TO ALL
       THE OTHER PROCESSES.
       ------------------------------------------------------------*/

    if ( !iam ) {
	/* Print the CPP definitions. */
	cpp_defs();
	
	/* Read the matrix stored on disk in *cpp */
	FILE	*fp = fopen( infile, "r" );

	if( !fp )
		ABORT("File does not exist");

	if( ftriple )
		dreadtriple(fp, &m, &n, &nnz, &a, &asub, &xa);
	else
		dreadhb_dist(iam, fp, &m, &n, &nnz, &a, &asub, &xa);
	
	printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);

	/* And read rhs */
    if (!(b=doubleMalloc_dist(m*nrhs))) ABORT("Malloc fails for b[]");
	int	k;
	for( k = 0; k < m; ++k )
		fscanf( fp, "%lf", b + k );

	/* Broadcast input data to the other PEs. */
	MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid.comm );
	MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );
	MPI_Bcast( b,    m,   MPI_DOUBLE, 0, grid.comm );
   } else {
	/* Receive matrix A from PE 0. */
	MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );

	/* Allocate storage for compressed column representation. */
	dallocateA_dist(n, nnz, &a, &asub, &xa);
	if (!(b=doubleMalloc_dist(m*nrhs))) ABORT("Malloc fails for b[]");

	MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid.comm );
	MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );
	MPI_Bcast( b,    m,   MPI_DOUBLE, 0, grid.comm );
    }

    /* Create compressed column matrix for A. */
    dCreate_CompCol_Matrix_dist(&A, m, n, nnz, a, asub, xa,
				SLU_NC, SLU_D, SLU_GE);

    /* Generate the exact solution and compute the right-hand side. */
    if (!(xtrue=doubleMalloc_dist(n*nrhs))) ABORT("Malloc fails for xtrue[]");
    *trans = 'N';
    ldx = n;
    ldb = m;
    dGenXtrue_dist(n, nrhs, xtrue, ldx);
//  dFillRHS_dist(trans, nrhs, xtrue, ldx, &A, b, ldb);


    if ( !(berr = doubleMalloc_dist(nrhs)) )
	ABORT("Malloc fails for berr[].");

    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM.
       ------------------------------------------------------------*/

    /* Set the default input options:
        options.Fact = DOFACT;
        options.Equil = YES;
        options.ColPerm = METIS_AT_PLUS_A;
        options.RowPerm = LargeDiag;
        options.ReplaceTinyPivot = YES;
        options.Trans = NOTRANS;
        options.IterRefine = DOUBLE;
        options.SolveInitialized = NO;
        options.RefineInitialized = NO;
        options.PrintStat = YES;
     */
    set_default_options_dist(&options);

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
		     &LUstruct, berr, &stat, &info);

	VMStats( stdout );

  /* Check the accuracy of the solution. */
    if ( !iam ) {

		if( info == 0 ) {

			dinf_norm_error_dist(n, nrhs, b, ldb, xtrue, ldx, &grid);

			FILE	*f = fopen( outfile, "w" );

			if( f ) {

				int	k;
				for( k = 0; k < n; ++k ) {
	
					if( k && !(k%6) )
						fprintf( f, "\n" );
	
					fprintf( f, "%.16e ", b[k] );
				}
	
				fprintf( f, "\n" );
				fclose( f );
			}
			else
				printf( "Can't open output file '%s'\n", outfile );
		}
		else
			printf("pdgssvx_ABglobal() returns INFO= %d\n", info);

		fclose( fopen( "slu_signal", "w" ) );
    }

   PStatPrint(&options, &stat, &grid);        /* Print the statistics. */

    /* ------------------------------------------------------------
       DEALLOCATE STORAGE.
       ------------------------------------------------------------*/
    PStatFree(&stat);
    Destroy_CompCol_Matrix_dist(&A);
    Destroy_LU(n, &grid, &LUstruct);
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);
    SUPERLU_FREE(b);
    SUPERLU_FREE(xtrue);
    SUPERLU_FREE(berr);

    /* ------------------------------------------------------------
       RELEASE THE SUPERLU PROCESS GRID.
       ------------------------------------------------------------*/
out:
    superlu_gridexit(&grid);

    /* ------------------------------------------------------------
       TERMINATES THE MPI EXECUTION ENVIRONMENT.
       ------------------------------------------------------------*/
    MPI_Finalize();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit main()");
#endif

}


int cpp_defs()
{
    printf(".. CPP definitions:\n");
#if ( PRNTlevel>=1 )
    printf("\tPRNTlevel = %d\n", PRNTlevel);
#endif
#if ( DEBUGlevel>=1 )
    printf("\tDEBUGlevel = %d\n", DEBUGlevel);
#endif
#if ( PROFlevel>=1 )
    printf("\tPROFlevel = %d\n", PROFlevel);
#endif
#if ( StaticPivot>=1 )
    printf("\tStaticPivot = %d\n", StaticPivot);
#endif
    printf("....\n");
    return 0;
}
