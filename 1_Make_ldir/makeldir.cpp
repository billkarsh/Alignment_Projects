

#include	"File.h"

#include	<stdio.h>
#include	<string.h>






// Deduce Z-ordering of images from given TrackEM2.xml file.
//
// While one might guess that true Z corresponds to the
// name tags "A1_", "B1_", etc. These image labels merely
// give the slice acquisition order. It is the order of
// these slice labels in the TrackEM2 file that gives
// true Z order.


int main( int argc, char **argv )
{
	char	line[2048], curKey[32];
	FILE	*f1	= FileOpenOrDie( argv[1], "r" );
	FILE	*f2	= FileOpenOrDie( "ldir", "w" );
	int		Z	= 0;

	curKey[0] = 0;

	for(;;) {

		if( !fgets( line, sizeof(line), f1 ) )
			goto close;

		if( !strstr( line, "file_path=" ) )
			continue;

		char *s = strrchr( line, '/' ) + 1;

		for( int i = 0; i < 10; ++i ) {

			if( !s[i] )
				goto close;

			if( s[i] == '_' ) {
				s[i] = 0;
				goto gotkey;
			}
		}

		goto close;

gotkey:
		if( !strcmp( curKey, s ) )
			continue;

		strcpy( curKey, s );
		fprintf( f2, "DIR %d %s_\n", Z++, s );
	}

close:
	if( f2 )
		fclose( f2 );

	if( f1 )
		fclose( f1 );

	return 0;
}


