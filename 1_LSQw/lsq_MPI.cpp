

#include	<mpi.h>		// put me ahead of all includes

#include	"lsq_MPI.h"

#include	"Cmdline.h"

#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

int	wkid = 0,	// my worker id (main=0)
    nwks = 1;	// total number workers






/* --------------------------------------------------------------- */
/* MPIInit ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set wkid, and set nwks which determines if MPI is used,
// or not in the case that nwks == 1.
//
void MPIInit( int& argc, char**& argv )
{
// Look for the -nwks=j argument

    for( int i = 1; i < argc; ++i ) {
        if( GetArg( &nwks, "-nwks=%d", argv[i] ) )
            break;
    }

// Get wkid

    int	thrdtype = -1;

    if( nwks > 1 ) {

        // Start MPI

        MPI_Init_thread( &argc, &argv,
            MPI_THREAD_FUNNELED, &thrdtype );

        MPI_Comm_rank( MPI_COMM_WORLD, &wkid );
    }

// Create log file 'lsqw_i'

    char slog[32];
    sprintf( slog, "lsqw_%d.txt", wkid );
    freopen( slog, "w", stdout );

    printf( "---- Read params ----\n" );

// MPI sanity checks

    if( nwks > 1 ) {

        // Verify worker count

        int	size;

        MPI_Comm_size( MPI_COMM_WORLD, &size );

        if( size != nwks ) {

            printf(
            "MPI: Bad worker count: %d, expected %d\n.",
            size, nwks );

            MPIExit();
            exit( 42 );
        }

        // Verify desired threading type.
        // Funneled means that the application can use
        // multiple threads but MPI communicates only
        // with the main thread.

        if( thrdtype != MPI_THREAD_FUNNELED ) {

            printf(
            "MPI: Not funneled thread type."
            " Try linking with -mt_mpi.\n" );

            MPIExit();
            exit( 42 );
        }
    }

    printf( "Worker %d / %d\n", wkid, nwks );
}

/* --------------------------------------------------------------- */
/* MPIExit ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void MPIExit()
{
    if( nwks > 1 )
        MPI_Finalize();
}

/* --------------------------------------------------------------- */
/* MPIWaitForOthers ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Caller is blocked here until all workers have also called it.
//
void MPIWaitForOthers()
{
    if( nwks > 1 )
        MPI_Barrier( MPI_COMM_WORLD );
}

/* --------------------------------------------------------------- */
/* MPISend ------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool MPISend( void* buf, int bytes, int wdst, int tag )
{
    return MPI_SUCCESS ==
    MPI_Ssend( buf, bytes, MPI_CHAR, wdst, tag, MPI_COMM_WORLD );
}

/* --------------------------------------------------------------- */
/* MPIRecv ------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool MPIRecv( void* buf, int bytes, int wsrc, int tag )
{
    MPI_Status	st;

    return MPI_SUCCESS ==
    MPI_Recv( buf, bytes, MPI_CHAR, wsrc, tag, MPI_COMM_WORLD, &st );
}


