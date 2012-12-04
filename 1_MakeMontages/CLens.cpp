

#include	"CLens.h"
#include	"GenDefs.h"
#include	"File.h"






/* --------------------------------------------------------------- */
/* ReadFile ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CLens::ReadFile( const char *path, FILE* flog )
{
	this->flog = flog;

	FILE		*f = fopen( path, "r" );
	CLineScan	LS;
	uint8		is[4] = {0,0,0,0};
	bool		ok = false;

// Header

	if( LS.Get( f ) <= 0 ) {

		fprintf( flog,
		"Lens::Read: Empty file [%s].\n", path );
		goto exit;
	}

// Entries

	while( LS.Get( f ) > 0 ) {

		TForm	T;
		int		cam;

		if( 7 != sscanf( LS.line,
				"%d%lf%lf%lf%lf%lf%lf",
				&cam,
				&T.t[0], &T.t[3], &T.t[1],
				&T.t[4], &T.t[2], &T.t[5] ) ) {

			fprintf( flog,
			"Lens::Read: Bad line [%s].\n", LS.line );
			goto exit;
		}

		if( cam < 0 || cam > 3 ) {

			fprintf( flog,
			"Lens::Read: Bad cam index [%s].\n", LS.line );
			goto exit;
		}

		Tf[cam] = T;
		is[cam] = 1;
	}

// Got all four?

	if( 4 != is[0] + is[1] + is[2] + is[3] ) {

		fprintf( flog,
		"Lens::Read: Missing entry [%s].\n", path );
		goto exit;
	}

// Make inverses

	for( int i = 0; i < 4; ++i ) {

		Tf[i].SetXY( 0, 0 );
		InvertTrans( Ti[i], Tf[i] );
	}

	ok = true;

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* Tdfm ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLens::Tdfm( TForm &T, int a, int b )
{
	if( a != b )
		MultiplyTrans( T, Ti[b], Tf[a] );
	else
		T.NUSetOne();
}

/* --------------------------------------------------------------- */
/* CamID --------------------------------------------------------- */
/* --------------------------------------------------------------- */

int CLens::CamID( const char *name )
{
	const char	*s = FileNamePtr( name );

	if( !(s = strstr( s, "_cam" )) ) {

		fprintf( flog,
		"Lens::PrintArg: Missing cam-id [%s].\n", name );
		exit( 42 );
	}

	return atoi( s + 4 );
}

/* --------------------------------------------------------------- */
/* PrintArg ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CLens::PrintArg(
	char		*buf,
	const char	*aname,
	const char	*bname )
{
	TForm	T;

	Tdfm( T, CamID( aname ), CamID( bname ) );

	sprintf( buf, " -Tdfm=%g,%g,%g,%g,%g,%g",
		T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5] );
}


