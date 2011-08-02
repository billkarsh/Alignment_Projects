

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
	arg.SCALE			= 999.0;
	arg.XSCALE			= 999.0;
	arg.YSCALE			= 999.0;
	arg.SKEW			= 999.0;
	arg.CTR				= 999.0;
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
		else if( GetArg( &arg.SCALE, "-SCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &arg.XSCALE, "-XSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &arg.YSCALE, "-YSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &arg.SKEW, "-SKEW=%lf", argv[i] ) )
			;
		else if( GetArg( &arg.CTR, "-CTR=%lf", argv[i] ) )
			;
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

// Get default parameters

	if( !ReadThmParams( thm, A.layer, stdout ) )
		return false;

// Commandline overrides

	printf( "\n---- Command-line overrides ----\n" );

	if( A.layer != B.layer ) {

		if( arg.SCALE != 999.0 && arg.SCALE != thm.SCALE ) {
			thm.SCALE = arg.SCALE;
			printf( "SCALE=%g\n", arg.SCALE );
		}

		if( arg.XSCALE != 999.0 && arg.XSCALE != thm.XSCALE ) {
			thm.XSCALE = arg.XSCALE;
			printf( "XSCALE=%g\n", arg.XSCALE );
		}

		if( arg.YSCALE != 999.0 && arg.YSCALE != thm.YSCALE ) {
			thm.YSCALE = arg.YSCALE;
			printf( "YSCALE=%g\n", arg.YSCALE );
		}

		if( arg.SKEW != 999.0 && arg.SKEW != thm.SKEW ) {
			thm.SKEW = arg.SKEW;
			printf( "SKEW=%g\n", arg.SKEW );
		}
	}

	if( arg.CTR != 999.0 )
		printf( "CTR=%g\n", arg.CTR );

	printf( "\n" );

// Fetch Til2Img entries

	printf( "\n---- Input images ----\n" );

	if( !ReadTil2Img( A.t2i, A.layer, A.tile, stdout ) ||
		!ReadTil2Img( B.t2i, B.layer, B.tile, stdout ) ) {

		return false;
	}

	PrintTil2Img( stdout, 'A', A.t2i );
	PrintTil2Img( stdout, 'B', B.t2i );

	printf( "\n" );

// Extract file name as useful label

	A.file = ExtractFilename( A.t2i.path );
	B.file = ExtractFilename( B.t2i.path );

	return true;
}


