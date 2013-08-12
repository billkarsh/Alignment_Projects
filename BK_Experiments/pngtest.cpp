

#include	"ImageIO.h"

#include	<stdlib.h>






int main( int argc, char **argv )
{
	const char*	sin = "/groups/apig/tomo/EX2/TIF1/A1_2_2.tif";
	uint16*		raster16;
	uint8*		raster8;
	uint32*		raster32;
	uint32		w=0, h=0;

	raster8 = Raster8FromTif( sin, w, h );
	raster16 = (uint16*)malloc( w*h*sizeof(uint16) );
	raster32 = (uint32*)malloc( w*h*sizeof(uint32) );
	for( int i = 0; i < w*h; ++i ) {
		raster16[i] = raster8[i];
		raster32[i] = raster8[i];
	}

	Raster8ToPng8( "png8out.png", raster8, w, h );
	Raster16ToPng16( "png16out.png", raster16, w, h );
	Raster32ToPngRGBA( "pngRGBAout.png", raster32, w, h );
	free( raster8 );
	free( raster16 );
	free( raster32 );

	w=0, h=0;
	raster16 = Raster16FromPng( "png8out.png", w, h, stdout );
	raster8 = (uint8*)malloc( w*h*sizeof(uint8) );
	for( int i = 0; i < w*h; ++i )
		raster8[i] = raster16[i];
	Raster8ToPng8( "VIEW8.png", raster8, w, h );
	free( raster8 );
	free( raster16 );

	w=0, h=0;
	raster16 = Raster16FromPng( "png16out.png", w, h, stdout );
	raster8 = (uint8*)malloc( w*h*sizeof(uint8) );
	for( int i = 0; i < w*h; ++i )
		raster8[i] = raster16[i];
	Raster8ToPng8( "VIEW16.png", raster8, w, h );
	free( raster8 );
	free( raster16 );

	w=0, h=0;
	raster16 = Raster16FromPng( "pngRGBAout.png", w, h, stdout );
	raster8 = (uint8*)malloc( w*h*sizeof(uint8) );
	for( int i = 0; i < w*h; ++i )
		raster8[i] = raster16[i];
	Raster8ToPng8( "VIEW32.png", raster8, w, h );
	free( raster8 );
	free( raster16 );

	return 0;
}


