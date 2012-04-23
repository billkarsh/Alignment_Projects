

#include	"File.h"

#include	<errno.h>
#include	<string.h>






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

/* --------------------------------------------------------------- */
/* FileNamePtr --------------------------------------------------- */
/* --------------------------------------------------------------- */

const char* FileNamePtr( const char *path )
{
	const char	*name;

	if( name = strrchr( path, '/' ) )
		++name;
	else
		name = path;

	return name;
}

/* --------------------------------------------------------------- */
/* FileDotPtr ---------------------------------------------------- */
/* --------------------------------------------------------------- */

const char* FileDotPtr( const char *path )
{
	const char	*dot;

	if( !(dot = strrchr( path, '.' )) )
		dot = path + strlen( path );

	return dot;
}

/* --------------------------------------------------------------- */
/* FileCloneNamePart --------------------------------------------- */
/* --------------------------------------------------------------- */

char* FileCloneNamePart( const char *path )
{
	char		buf[2048];
	const char	*name	= FileNamePtr( path ),
				*dot	= FileDotPtr( name );

	sprintf( buf, "%.*s", dot - name, name );

	return strdup( buf );
}


