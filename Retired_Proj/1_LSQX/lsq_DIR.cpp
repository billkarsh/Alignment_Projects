

#include	"lsq_DIR.h"

#include	"File.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* ReadDIRFile --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read file who's entries associate:
//
//	[layer-#] <-> [imgname-tag].
//
// Example entry:
//
//	DIR 0 A1_
//
void DIR::ReadDIRFile( const char *dirfile, FILE *FOUT )
{
    printf( "---- Read DIR ----\n" );

    if( !dirfile ) {
        printf( "No DIR file - assuming layer zero only.\n" );
        return;
    }

// Read file

    FILE		*f = FileOpenOrDie( dirfile, "r" );
    CLineScan	LS;

    for(;;) {

        if( LS.Get( f ) <= 0 )
            break;

        if( !strncmp( LS.line, "DIR", 3 ) ) {

            fprintf( FOUT, LS.line );

            char	buf[1024];
            int		z;

            sscanf( LS.line, "DIR %d %[^\r\n]", &z, buf );
            dirTbl[string(buf)] = z;
        }
        else {
            printf( "Bad line '%s' in DIR file.\n", LS.line );
            exit( 42 );
        }
    }

    fclose( f );
}

/* --------------------------------------------------------------- */
/* ZFromName ----------------------------------------------------- */
/* --------------------------------------------------------------- */

int DIR::ZFromName( const char *name ) const
{
    map<string,int>::const_iterator	it = dirTbl.find( string(name) );

    if( it != dirTbl.end() )
        return it->second;

    return 0;
}


