

#include	"EZThreads.h"

#include	<limits.h>

#include	<vector>
using namespace std;






/* --------------------------------------------------------------- */
/* EZThreads ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Run nthr instances of given EZThreadproc
// and await completion of all threads.
//
// stksize_factor:
// <= 0:	use default stack size (~2MB)
//    1:	use 1 x PTHREAD_STACK_MIN (~16KB)
//    n:	use n x PTHREAD_STACK_MIN.
//
// msgname:
// Name of proc for use in error message.
//
// Return true if launches successful.
//
bool EZThreads(
    EZThreadproc	proc,
    int				nthr,
    int				stksize_factor,
    const char		*msgname,
    FILE			*flog )
{
    vector<pthread_t>	vthr( nthr );

// Start coworkers [1..nthr)

    if( nthr > 1 ) {

        pthread_attr_t	attr;
        pthread_attr_init( &attr );

        pthread_attr_setdetachstate( &attr,
            PTHREAD_CREATE_JOINABLE );

        if( stksize_factor > 0 ) {

            pthread_attr_setstacksize( &attr,
                stksize_factor * PTHREAD_STACK_MIN );
        }

        int	err;

        for( int i = 1; i < nthr; ++i ) {

            err = pthread_create(
                    &vthr[i], &attr,
                    proc, reinterpret_cast<void*>(i) );

            if( err ) {

                fprintf( flog,
                "Error [%d] starting '%s' thread, index [%d].\n",
                err, msgname, i );

                for( int j = 1; j < i; ++j )
                    pthread_cancel( vthr[j] );
            }
        }

        pthread_attr_destroy( &attr );

        if( err )
            return false;
    }

// Run worker 0 locally

    proc( 0 );

// Wait for coworkers

    if( nthr > 1 ) {

        for( int i = 1; i < nthr; ++i ) {
            pthread_join( vthr[i], NULL );
            pthread_detach( vthr[i] );
        }
    }

    return true;
}


