

#include	"CCropLoad.h"
#include	"ImageIO.h"
#include	"File.h"






/* --------------------------------------------------------------- */
/* ReadIDB ------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CCropLoad::ReadIDB( const string &idb, FILE* flog )
{
	this->flog = flog;

// Exists?

	if( !(isfile = !idb.empty()) )
		return true;

	char		path[2048];
	sprintf( path, "%s/crop.txt", idb.c_str() );

	FILE		*f = fopen( path, "r" );
	CLineScan	LS;
	uint8		is[4] = {0,0,0,0};
	bool		ok = false;

// Exists?

	if( !(isfile = (f != NULL)) )
		return true;

// Header

	//if( LS.Get( f ) <= 0 ) {

	//	fprintf( flog,
	//	"Crop::Read: Empty file [%s].\n", path );
	//	goto exit;
	//}

// Entries

	while( LS.Get( f ) > 0 ) {

		CBox	b;
		int		cam;

		if( 5 != sscanf( LS.line, "%d%d%d%d%d",
				&cam, &b.x0, &b.y0, &b.dx, &b.dy ) ) {

			fprintf( flog,
			"Crop::Read: Bad line [%s].\n", LS.line );
			goto exit;
		}

		if( cam < 0 || cam > 3 ) {

			fprintf( flog,
			"Crop::Read: Bad cam index [%s].\n", LS.line );
			goto exit;
		}

		B[cam] = b;
		is[cam] = 1;
	}

// Got all four?

	if( 4 != is[0] + is[1] + is[2] + is[3] ) {

		fprintf( flog,
		"Crop::Read: Missing entry [%s].\n", path );
		goto exit;
	}

	ok = true;

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* Raster8 ------------------------------------------------------- */
/* --------------------------------------------------------------- */

uint8* CCropLoad::Raster8(
	const char*	name,
	int			cam,
	uint32		&w,
	uint32		&h,
	FILE*		flog,
	bool		transpose )
{
	uint8*	ras = Raster8FromAny( name, w, h, flog, transpose );

	if( ras && isfile ) {

		int		dx		= B[cam].dx,
				dy		= B[cam].dy;
		uint8	*dst	= ras,
				*src	= ras + B[cam].x0 + w*B[cam].y0;

		for( int y = 0; y < dy; ++y, dst += dx, src += w )
			memmove( dst, src, dx );
	}

	return ras;
}


