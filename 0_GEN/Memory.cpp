

#include	"Memory.h"
#include	"File.h"

#include	<string.h>
#include	<sys/resource.h>
#include	<unistd.h>






/* --------------------------------------------------------------- */
/* VMStats ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Note: Google "/proc/PID/status" to find many other examples
// of available kernel info.
//
void VMStats( FILE *flog )
{
    char			host[128];
    struct rusage	usage;
    int				pid = getpid();

    gethostname( host, sizeof(host) );
    getrusage( RUSAGE_SELF, &usage );

    fprintf( flog, "\n---- Memory ----\n" );
    fprintf( flog, "Host name: %s\n", host );
    fprintf( flog, "PID:       %d\n", pid );
    fprintf( flog, "User time:   %5ld seconds.\n",
        usage.ru_utime.tv_sec );
    fprintf( flog, "System time: %5ld seconds.\n",
        usage.ru_stime.tv_sec );

    CLineScan	LS;

    LS.bufsize	= 1024;
    LS.line		= (char*)malloc( LS.bufsize );

    sprintf( LS.line, "/proc/%d/status", pid );

    FILE	*f = fopen( LS.line, "r" );

    if( f ) {

        while( LS.Get( f ) > 0 ) {

            if( strstr( LS.line, "Vm" ) || strstr( LS.line, "ctxt" ) )
                fprintf( flog, "%s", LS.line );
        }

        fclose( f );
    }

    fprintf( flog, "\n" );
}


