#include	"superlu_ddefs.h"
#include	<sys/resource.h>


void VMStats( FILE *flog )
{
	char			host[128];
	struct rusage	usage;

	gethostname( host, sizeof(host) );
	getrusage( RUSAGE_SELF, &usage );

	fprintf( flog, "\n---- Memory ----\n" );
	fprintf( flog, "Host name: %s\n", host );
	fprintf( flog, "PID:       %d\n", getpid() );
	fprintf( flog, "User time:   %5d seconds.\n", usage.ru_utime );
	fprintf( flog, "System time: %5d seconds.\n", usage.ru_stime );

	size_t	bufsize	= 1024;
	char	*line	= (char*)malloc( bufsize );

	sprintf( line, "/proc/%d/status", getpid() );

	FILE	*f = fopen( line, "r" );

	if( f ) {

		while( getline( &line, &bufsize, f ) > 0 ) {

			if( strstr( line, "Vm" ) || strstr( line, "ctxt" ) )
				fprintf( flog, "%s", line );
		}

		fclose( f );
	}

	free( line );

	fprintf( flog, "\n" );
	fflush( flog );
}


