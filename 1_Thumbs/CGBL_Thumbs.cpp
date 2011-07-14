

#include	"CGBL_Thumbs.h"

#include	"Cmdline.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_Thumbs	GBL;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Object management --------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_Thumbs::CGBL_Thumbs()
{
	arg.Transpose		= false;
	arg.NoFolds			= false;
	arg.SingleFold		= false;

	A.layer	= 0;
	A.tile	= 0;

	B.layer	= 0;
	B.tile	= 0;
}

/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CGBL_Thumbs::SetCmdLine( int argc, char* argv[] )
{
// Parse args

	vector<char*>	noa;	// non-option arguments

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			noa.push_back( argv[i] );
		else if( IsArg( "-tr", argv[i] ) )
			arg.Transpose = true;
		else if( IsArg( "-nf", argv[i] ) )
			arg.NoFolds = true;
		else if( IsArg( "-sf", argv[i] ) )
			arg.SingleFold = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			return false;
		}
	}

// Decode labels in noa[0]

	if( !noa.size() ||
		(4 != sscanf( noa[0],
				"%d/%d@%d/%d",
				&A.layer, &A.tile,
				&B.layer, &B.tile )) ) {

		printf( "main: Usage: thumbs <za/ia@zb/ib>.\n" );
		return false;
	}

// Rename stdout using image labels

	OpenPairLog( A.layer, A.tile, B.layer, B.tile );

	printf( "\n---- Thumbnail matching ----\n" );

// Record start time

	time_t	t0 = time( NULL );
	printf( "main: Start: %s\n", ctime(&t0) );

// Get parameters

	if( !ReadThmParams( thm, A.layer, stdout ) )
		return false;

// Fetch Til2Img entries

	if( !ReadTil2Img( A.t2i, A.layer, A.tile, stdout ) ||
		!ReadTil2Img( B.t2i, B.layer, B.tile, stdout ) ) {

		return false;
	}

	PrintTil2Img( stdout, 'A', A.t2i );
	PrintTil2Img( stdout, 'B', B.t2i );

// Extract file name as useful label

	A.file = ExtractFilename( A.t2i.path );
	B.file = ExtractFilename( B.t2i.path );

	return true;
}


