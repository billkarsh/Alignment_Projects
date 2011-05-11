

#include	"ImageIO.h"
#include	"Draw.h"






static void DrawImage()
{
	int		w = 1376, h = 1040, n = w * h;
	uint8*	r = (uint8*)RasterAlloc( n );

	memset( r, 0, n );

	for( int y = 9; y <= h-9; y += 44 ) {
		for( int x = 9; x <= w-9; x += 44 ) {
			for( int i = -8; i <= 8; ++i ) {
				for( int j = -8; j <= 8; ++j ) {
					r[x+j + w*(y+i)] = 255;
				}
			}
		}
	}

	Raster8ToTif8( "Lines.tif", r, w, h );

	RasterFree( r );
}


int main( int argc, char **argv )
{
	DrawImage();

	return 0;
}


