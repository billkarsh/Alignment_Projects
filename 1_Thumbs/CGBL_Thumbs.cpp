

#include	"CGBL_Thumbs.h"

#include	"Cmdline.h"
#include	"File.h"
#include	"CAffineLens.h"
#include	"Debug.h"


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
	_arg.SCALE			= 999.0;
	_arg.XSCALE			= 999.0;
	_arg.YSCALE			= 999.0;
	_arg.SKEW			= 999.0;
	_arg.ima			= NULL;
	_arg.imb			= NULL;
	_arg.FLD			= 0;
	_arg.MODE			= 0;

	arg.CTR				= 999.0;
	arg.fma				= NULL;
	arg.fmb				= NULL;
	arg.Transpose		= false;
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

	vector<double>	vD;
	char			*key;

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			key = argv[i];
		else if( GetArgList( vD, "-Tdfm=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tdfm.push_back( TAffine( &vD[0] ) );
			else {
				printf(
				"main: WARNING: Bad format in -Tdfm [%s].\n",
				argv[i] );
			}
		}
		else if( GetArgList( vD, "-Tab=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tab.push_back( TAffine( &vD[0] ) );
			else {
				printf(
				"main: WARNING: Bad format in -Tab [%s].\n",
				argv[i] );
			}
		}
		else if( GetArg( &_arg.SCALE, "-SCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.XSCALE, "-XSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.YSCALE, "-YSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.SKEW, "-SKEW=%lf", argv[i] ) )
			;
		else if( GetArgStr( _arg.ima, "-ima=", argv[i] ) )
			;
		else if( GetArgStr( _arg.imb, "-imb=", argv[i] ) )
			;
		else if( GetArg( &_arg.FLD, "-FLD=%c", argv[i] ) )
			;
		else if( GetArg( &_arg.MODE, "-MODE=%c", argv[i] ) )
			;
		else if( GetArg( &arg.CTR, "-CTR=%lf", argv[i] ) )
			;
		else if( GetArgStr( arg.fma, "-fma=", argv[i] ) )
			;
		else if( GetArgStr( arg.fmb, "-fmb=", argv[i] ) )
			;
		else if( IsArg( "-tr", argv[i] ) )
			arg.Transpose = true;
		else if( IsArg( "-nf", argv[i] ) )
			_arg.FLD = 'N';
		else if( IsArg( "-sf", argv[i] ) )
			arg.SingleFold = true;
		else if( IsArg( "-dbgcor", argv[i] ) )
			dbgCor = true;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			return false;
		}
	}

// Decode labels in key

	if( !key ||
		(4 != sscanf( key, "%d/%d@%d/%d",
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

	if( !ReadMatchParams( mch, A.layer, B.layer ) )
		return false;

// Which file params to use according to (same,cross) layer

	double	cSCALE=1, cXSCALE=1, cYSCALE=1, cSKEW=0;
	int		cDfmFromTab;

	ctx.FLD = mch.FLD;

	if( A.layer == B.layer ) {

		cDfmFromTab	= mch.TAB2DFM_SL;

		//ctx.Tdfm = identity (default)
		ctx.XYCONF	= mch.XYCONF_SL;
		ctx.NBMXHT	= mch.NBMXHT_SL;
		ctx.HFANGDN	= mch.HFANGDN_SL;
		ctx.HFANGPR	= mch.HFANGPR_SL;
		ctx.RTRSH	= mch.RTRSH_SL;
		ctx.OLAP2D	= mch.OLAP2D_SL;
		ctx.MODE	= mch.MODE_SL;
		ctx.THMDEC	= mch.THMDEC_SL;
		ctx.OLAP1D	= mch.OLAP1D_SL;
		ctx.LIMXY	= mch.LIMXY_SL;
	}
	else {

		cSCALE	= mch.SCALE;
		cXSCALE	= mch.XSCALE;
		cYSCALE	= mch.YSCALE;
		cSKEW	= mch.SKEW;

		ctx.Tdfm.ComposeDfm( cSCALE, cXSCALE, cYSCALE, 0, cSKEW );

		cDfmFromTab	= mch.TAB2DFM_XL;

		ctx.XYCONF	= mch.XYCONF_XL;
		ctx.NBMXHT	= mch.NBMXHT_XL;
		ctx.HFANGDN	= mch.HFANGDN_XL;
		ctx.HFANGPR	= mch.HFANGPR_XL;
		ctx.RTRSH	= mch.RTRSH_XL;
		ctx.OLAP2D	= mch.OLAP2D_XL;
		ctx.MODE	= mch.MODE_XL;
		ctx.THMDEC	= mch.THMDEC_XL;
		ctx.OLAP1D	= mch.OLAP1D_XL;
		ctx.LIMXY	= mch.LIMXY_XL;
	}

// Fetch Til2Img entries

	printf( "\n---- Input images ----\n" );

	IDBReadImgParams( idb );

	if( !IDBTil2Img( A.t2i, idb, A.layer, A.tile, _arg.ima ) ||
		!IDBTil2Img( B.t2i, idb, B.layer, B.tile, _arg.imb ) ) {

		return false;
	}

	PrintTil2Img( stdout, 'A', A.t2i );
	PrintTil2Img( stdout, 'B', B.t2i );

	printf( "\n" );

// Extract file name as useful label

	A.file = FileCloneNamePart( A.t2i.path.c_str() );
	B.file = FileCloneNamePart( B.t2i.path.c_str() );

// Commandline overrides

	printf( "\n---- Command-line overrides ----\n" );

	if( _arg.Tab.size() ) {

		Tab = _arg.Tab[0];

		// remove lens parts of Tab coming from cross_thisblock

		if( mch.PXLENS && A.layer != B.layer ) {

			CAffineLens	LN;

			if( !LN.ReadIDB( idb ) )
				return false;

			LN.UpdateTFormRHS( Tab, A.t2i.path.c_str(), true );
			LN.UpdateTFormLHS( Tab, B.t2i.path.c_str(), false );
		}

		Tab.TPrint( stdout, "Tab= " );
	}
	else
		Tab.FromAToB( A.t2i.T, B.t2i.T );

	int	altTdfm = false;

	if( _arg.Tdfm.size() ) {

		ctx.Tdfm	= _arg.Tdfm[0];
		altTdfm		= true;
	}
	else {

		if( _arg.SCALE != 999.0 ) {
			cSCALE	= _arg.SCALE;
			altTdfm	= true;
			printf( "SCALE=%g\n", _arg.SCALE );
		}

		if( _arg.XSCALE != 999.0 ) {
			cXSCALE	= _arg.XSCALE;
			altTdfm	= true;
			printf( "XSCALE=%g\n", _arg.XSCALE );
		}

		if( _arg.YSCALE != 999.0 ) {
			cYSCALE	= _arg.YSCALE;
			altTdfm	= true;
			printf( "YSCALE=%g\n", _arg.YSCALE );
		}

		if( _arg.SKEW != 999.0 ) {
			cSKEW	= _arg.SKEW;
			altTdfm	= true;
			printf( "SKEW=%g\n", _arg.SKEW );
		}

		if( altTdfm )
			ctx.Tdfm.ComposeDfm( cSCALE, cXSCALE, cYSCALE, 0, cSKEW );
	}

	if( !altTdfm && cDfmFromTab ) {

		TAffine	R;
		R.NUSetRot( -Tab.GetRadians() );

		ctx.Tdfm = Tab;
		ctx.Tdfm.SetXY( 0, 0 );
		ctx.Tdfm = R * ctx.Tdfm;
	}

	ctx.Tdfm.TPrint( stdout, "Tdfm=" );

	if( _arg.FLD ) {
		ctx.FLD = _arg.FLD;
		printf( "FLD=%c\n", _arg.FLD );
	}

	if( ctx.FLD = 'X' ) {
		ctx.FLD = (GBL.A.layer == GBL.B.layer ? 'N' : 'Y');
		printf( "FLD=%c (was X)\n", ctx.FLD );
	}

	if( _arg.MODE ) {
		ctx.MODE = _arg.MODE;
		printf( "MODE=%c\n", _arg.MODE );
	}

	if( ctx.MODE == 'Z' ) {
		ctx.MODE = 'C';
		arg.CTR = 0.0;
		printf( "MODE=C (was Z)\n", arg.CTR );
	}
	else if( ctx.MODE == 'M' ) {
		ctx.MODE = 'N';
		arg.CTR = 0.0;
		printf( "MODE=N (was M)\n", arg.CTR );
	}

	if( arg.CTR != 999.0 )
		printf( "CTR=%g\n", arg.CTR );

	printf( "\n" );

	return true;
}


