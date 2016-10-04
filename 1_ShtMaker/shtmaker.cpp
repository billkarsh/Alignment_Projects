

#include	<stdio.h>






// Read file stdin and emit C-code to print that content
// line by line.
//
int main(int argc, char **argv)
{
    char	*line;
    size_t	size = 0;

    while( getline( &line, &size, stdin ) > 0 ) {

        printf( "\tfprintf( f, \"" );

        char	*s = line, c;

        while( c = *s++ ) {

            switch( c ) {

                case '\n':
                    printf( "\\n\" );\n" );
                break;

                case '\t':
                    printf( "\\t" );
                break;

                case '"':
                    printf( "\\\"" );
                break;

                case '\\':
                    printf( "\\\\" );
                break;

                case '%':
                    printf( "%%%%" );
                break;

                default:
                    printf( "%c", c );
                break;
            }
        }
    }

    return 0;
}


