

#include	"File.h"

#include	<cstdlib>
using namespace std;






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
		"Can't open file [%s] for op [%s].\n", name, rw );
		exit( 42 );
	}

	return f;
}


