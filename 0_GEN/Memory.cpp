

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
	struct rusage	usage;

	getrusage( RUSAGE_SELF, &usage );

	fprintf( flog, "\n---- Memory ----\n" );
	fprintf( flog, "User time:   %5d seconds.\n", usage.ru_utime );
	fprintf( flog, "System time: %5d seconds.\n", usage.ru_stime );

	CLineScan	LS;

	LS.bufsize	= 1024;
	LS.line		= (char*)malloc( LS.bufsize );

	sprintf( LS.line, "/proc/%d/status", getpid() );

	FILE	*f = fopen( LS.line, "r" );

	if( f ) {

		while( LS.Get( f ) > 0 ) {

			if( strstr( LS.line, "Vm" ) || strstr( LS.line, "ctxt" ) )
				fprintf( flog, "%s", LS.line );
		}

		fclose( f );
	}
}


