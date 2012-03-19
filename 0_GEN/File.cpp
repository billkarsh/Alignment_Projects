

#include	"File.h"

#include	<errno.h>






/* --------------------------------------------------------------- */
/* FileOpenOrDie ------------------------------------------------- */
/* --------------------------------------------------------------- */

FILE *FileOpenOrDie(
	const char	*name,
	const char	*rw,
	FILE		*flog )
{
	FILE	*f = fopen( name, rw );

	if( !f ) {

		if( !flog )
			flog = stdout;

		fprintf( flog,
		"Can't open file [%s] for op [%s] errno [%d].\n",
		name, rw, errno );

		exit( 42 );
	}

	return f;
}


