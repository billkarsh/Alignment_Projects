

#include	"ImageIO.h"

#include	<stdlib.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main(int argc, char **argv)
{
    if( argc < 2 ) {
        printf("Usage: qa <list of fm files>\n");
        exit( 42 );
    }

    for(int i=1; i<argc; i++) {

        uint32	w, h;
        uint8*	a = Raster8FromTif( argv[i], w, h );
        uint32	np = w*h;
        uint32	ps[256];
        uint32	j;

        memset( ps, 0, 256 * sizeof(uint32) );

        for( j = 0; j < np; ++j )
            ++ps[a[j]];

        int		npa = 0;  // number of patches
        uint32	nnz = 0;  // pixels in all of them

        for( j = 0; j < 256; ++j ) {

            if( ps[j] ) {
                ++npa;
                nnz += ps[j];
            }
        }

        printf( "file %s, %4d patches, %6.1f percent covered\n",
        argv[i], npa, double(nnz)/double(np)*100.0);

        RasterFree( a );
    }

    return 0;
}


