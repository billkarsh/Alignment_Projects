

#include	"CGBL_dmesh.h"

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

CGBL_dmesh	GBL;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Object management --------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_dmesh::CGBL_dmesh()
{
	arg.SCALE			= 999.0;
	arg.XSCALE			= 999.0;
	arg.YSCALE			= 999.0;
	arg.SKEW			= 999.0;
	arg.CTR				= 999.0;
	arg.Transpose		= false;
	arg.WithinSection	= false;
	arg.Verbose			= false;
	arg.NoFolds			= false;
	arg.SingleFold		= false;
	arg.Heatmap			= false;

	A.layer	= 0;
	A.tile	= 0;

	B.layer	= 0;
	B.tile	= 0;
}

/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CGBL_dmesh::SetCmdLine( int argc, char* argv[] )
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
		else if( IsArg( "-ws", argv[i] ) )
			arg.WithinSection = true;
		else if( IsArg( "-nf", argv[i] ) )
			arg.NoFolds = true;
		else if( IsArg( "-sf", argv[i] ) )
			arg.SingleFold = true;
		else if( IsArg( "-v", argv[i] ) )
			arg.Verbose = true;
		else if( IsArg( "-heatmap", argv[i] ) )
			arg.Heatmap = true;
		else if( !strncmp( argv[i], "-TRA=", 5 ) ) {

			TForm	a;

			if( 6 != sscanf(
				argv[i] + 5, "%lf,%lf,%lf,%lf,%lf,%lf",
				&a.t[0], &a.t[1], &a.t[2],
				&a.t[3], &a.t[4], &a.t[5] ) ) {

				printf(
				"main: WARNING: Bad format in -TRA [%s].\n",
				argv[i] );
			}
			else
				Tusr.push_back( a );
		}
		else if( !strncmp( argv[i], "-EXY=", 5 ) ) {

			Point	a;

			if( 2 != sscanf(
				argv[i] + 5, "%lf,%lf",
				&a.x, &a.y ) ) {

				printf(
				"main: WARNING: Bad format in -EXY [%s].\n",
				argv[i] );
			}
			else
				XYusr.push_back( a );
		}
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

		printf( "main: Usage: ptest <za/ia@zb/ib>.\n" );
		return false;
	}

// Rename stdout using image labels

	OpenPairLog( A.layer, A.tile, B.layer, B.tile );

	printf( "\n---- dmesh start ----\n" );

// Record start time

	time_t	t0 = time( NULL );
	printf( "main: Start: %s\n", ctime(&t0) );

// Get default parameters

	if( !ReadThmParams( thm, A.layer, stdout ) )
		return false;

	if( !ReadMeshParams( msh, A.layer, stdout ) )
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


