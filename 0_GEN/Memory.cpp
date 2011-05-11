

#include	"Memory.h"

#include	<stdlib.h>
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
	struct rusage	usage;

	getrusage( RUSAGE_SELF, &usage );

	fprintf( flog, "\n---- Memory ----\n" );
	fprintf( flog, "User time:   %5d seconds.\n", usage.ru_utime );
	fprintf( flog, "System time: %5d seconds.\n", usage.ru_stime );

	char	line[1024];

	sprintf( line, "/proc/%d/status", getpid() );

	FILE	*f = fopen( line, "r" );

	if( f ) {

		while( fgets( line, sizeof(line), f ) ) {

			if( strstr( line, "Vm" ) || strstr( line, "ctxt" ) )
				fprintf( flog, "%s", line );
		}

		fclose( f );
	}
}


