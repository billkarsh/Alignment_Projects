

#include	"CCropMask.h"
#include	"Disk.h"
#include	"File.h"






/* --------------------------------------------------------------- */
/* IsFile -------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CCropMask::IsFile( const string &idb )
{
	if( isfile < 0 ) {

		char	path[2048];
		sprintf( path, "%s/crop.txt", idb.c_str() );
		isfile = DskExists( path );
	}

	return isfile;
}

/* --------------------------------------------------------------- */
/* ReadIDB ------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CCropMask::ReadIDB( const string &idb, FILE* flog )
{
	this->flog = flog;

// Exists?

	if( !(isfile = !idb.empty()) )
		return true;

	char	path[2048];
	sprintf( path, "%s/crop.txt", idb.c_str() );

	FILE		*f = fopen( path, "r" );
	CLineScan	LS;
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

		int	cam, x0, y0, dx, dy;

		if( 5 != sscanf( LS.line, "%d%d%d%d%d",
				&cam, &x0, &y0, &dx, &dy ) ) {

			fprintf( flog,
			"Crop::Read: Bad line [%s].\n", LS.line );
			goto exit;
		}

		if( cam < 0 || cam > 3 ) {

			fprintf( flog,
			"Crop::Read: Bad cam index [%s].\n", LS.line );
			goto exit;
		}

		B[cam].L	= x0;
		B[cam].R	= x0 + dx;
		B[cam].B	= y0;
		B[cam].T	= y0 + dy;
	}

	ok = true;

exit:
	if( f )
		fclose( f );

	return ok;
}


