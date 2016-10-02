

#include	"CGBL_dmesh.h"

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

CGBL_dmesh	GBL;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Object management --------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_dmesh::CGBL_dmesh()
{
	_arg.SCALE			= 999.0;
	_arg.XSCALE			= 999.0;
	_arg.YSCALE			= 999.0;
	_arg.SKEW			= 999.0;
	_arg.matchparams	= NULL;
	_arg.ima			= NULL;
	_arg.imb			= NULL;
	_arg.FLD			= 0;
	_arg.MODE			= 0;

	arg.CTR				= 999.0;
	arg.fma				= NULL;
	arg.fmb				= NULL;
	arg.comp_png		= NULL;
	arg.registered_png	= NULL;
	arg.Transpose		= false;
	arg.WithinSection	= false;
	arg.SingleFold		= false;
	arg.JSON			= false;
	arg.Verbose			= false;
	arg.Heatmap			= false;

	A.z		= 0;
	A.id	= 0;

	B.z		= 0;
	B.id	= 0;
}

/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CGBL_dmesh::SetCmdLine( int argc, char* argv[] )
{
// Parse args

	const char		*key = NULL;
	vector<double>	vD;

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			key = argv[i];
		else if( GetArgStr( _arg.matchparams, "-prm=", argv[i] ) )
			;
		else if( GetArgList( vD, "-Tdfm=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tdfm.push_back( TAffine( &vD[0] ) );
			else {
				fprintf( stderr,
				"main: WARNING: Bad format in -Tdfm [%s].\n",
				argv[i] );
			}
		}
		else if( GetArgList( vD, "-Tab=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tab.push_back( TAffine( &vD[0] ) );
			else {
				fprintf( stderr,
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
		else if( GetArgStr( arg.comp_png, "-comp_png=", argv[i] ) )
			;
		else if( GetArgStr( arg.registered_png, "-registered_png=", argv[i] ) )
			;
		else if( IsArg( "-tr", argv[i] ) )
			arg.Transpose = true;
		else if( IsArg( "-ws", argv[i] ) )
			arg.WithinSection = true;
		else if( IsArg( "-nf", argv[i] ) )
			_arg.FLD = 'N';
		else if( IsArg( "-sf", argv[i] ) )
			arg.SingleFold = true;
		else if( IsArg( "-json", argv[i] ) )
			arg.JSON = true;
		else if( IsArg( "-v", argv[i] ) )
			arg.Verbose = true;
		else if( IsArg( "-heatmap", argv[i] ) )
			arg.Heatmap = true;
		else if( IsArg( "-dbgcor", argv[i] ) )
			dbgCor = true;
		else if( GetArgList( vD, "-Tmsh=", argv[i] ) ) {

			if( 6 == vD.size() )
				Tmsh.push_back( TAffine( &vD[0] ) );
			else {
				fprintf( stderr,
				"main: WARNING: Bad format in -Tmsh [%s].\n",
				argv[i] );
			}
		}
		else if( GetArgList( vD, "-XYexp=", argv[i] ) ) {

			if( 2 == vD.size() )
				XYexp.push_back( Point( vD[0], vD[1] ) );
			else {
				fprintf( stderr,
				"main: WARNING: Bad format in -XYexp [%s].\n",
				argv[i] );
			}
		}
		else {
			fprintf( stderr,
			"Did not understand option '%s'.\n", argv[i] );
			return false;
		}
	}

// Decode labels in key

	if( !key ||
		(4 != sscanf( key, "%d.%d^%d.%d",
				&A.z, &A.id,
				&B.z, &B.id )) ) {

		fprintf( stderr, "main: Usage: ptest <za.ia^zb.ib>.\n" );
		return false;
	}

// Rename stdout using image labels

	fprintf( stderr, "\n---- dmesh start ----\n" );

// Record start time

	time_t	t0 = time( NULL );
	fprintf( stderr, "main: Start: %s\n", ctime(&t0) );

// Get default parameters

	if( !ReadMatchParams( mch, A.z, B.z, _arg.matchparams, stderr ) )
		return false;

// Which file params to use according to (same,cross) layer

	double	cSCALE=1, cXSCALE=1, cYSCALE=1, cSKEW=0;
	int		cDfmFromTab;

	ctx.FLD = mch.FLD;

	if( A.z == B.z ) {

		cDfmFromTab	= mch.TAB2DFM_SL;

		//ctx.Tdfm = identity (default)
		ctx.XYCONF	= mch.XYCONF_SL;
		ctx.NBMXHT	= mch.NBMXHT_SL;
		ctx.HFANGDN	= mch.HFANGDN_SL;
		ctx.HFANGPR	= mch.HFANGPR_SL;
		ctx.RTRSH	= mch.RTRSH_SL;
		ctx.RIT		= mch.RIT_SL;
		ctx.RFA		= mch.RFA_SL;
		ctx.RFT		= mch.RFT_SL;
		ctx.OLAP2D	= mch.OLAP2D_SL;
		ctx.MODE	= mch.MODE_SL;
		ctx.THMDEC	= mch.THMDEC_SL;
		ctx.OLAP1D	= mch.OLAP1D_SL;
		ctx.LIMXY	= mch.LIMXY_SL;
		ctx.OPT		= mch.OPT_SL;
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
		ctx.RIT		= mch.RIT_XL;
		ctx.RFA		= mch.RFA_XL;
		ctx.RFT		= mch.RFT_XL;
		ctx.OLAP2D	= mch.OLAP2D_XL;
		ctx.MODE	= mch.MODE_XL;
		ctx.THMDEC	= mch.THMDEC_XL;
		ctx.OLAP1D	= mch.OLAP1D_XL;
		ctx.LIMXY	= mch.LIMXY_XL;
		ctx.OPT		= true;
	}

// Fetch Til2Img entries

	fprintf( stderr, "\n---- Input images ----\n" );

	IDBFromTemp( idb, "../../" );

	if( !IDBT2IGet1( A.t2i, idb, A.z, A.id, _arg.ima ) ||
		!IDBT2IGet1( B.t2i, idb, B.z, B.id, _arg.imb ) ) {

		return false;
	}

	PrintTil2Img( stderr, 'A', A.t2i );
	PrintTil2Img( stderr, 'B', B.t2i );

	fprintf( stderr, "\n" );

// Commandline overrides

	fprintf( stderr, "\n---- Command-line overrides ----\n" );

	if( _arg.Tab.size() ) {

		Tab = _arg.Tab[0];

		// remove lens parts of Tab coming from cross_thisblock

		if( mch.PXLENS && A.z != B.z ) {

			CAffineLens	LN;

			if( !LN.ReadIDB( idb ) )
				return false;

			LN.UpdateTFormRHS( Tab, A.t2i.cam, true );
			LN.UpdateTFormLHS( Tab, B.t2i.cam, false );
		}

		Tab.TPrint( stderr, "Tab= " );
	}
	else if( _arg.Ta.size() || _arg.Tb.size() ) {

		TAffine	Ta, Tb;

		if( _arg.Ta.size() )
			Ta = _arg.Ta[0];

		if( _arg.Tb.size() )
			Tb = _arg.Tb[0];

		Tab.FromAToB( Ta, Tb );

		Ta.TPrint( stderr,  "Ta=  " );
		Tb.TPrint( stderr,  "Tb=  " );
		Tab.TPrint( stderr, "Tab= " );
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
			fprintf( stderr, "SCALE=%g\n", _arg.SCALE );
		}

		if( _arg.XSCALE != 999.0 ) {
			cXSCALE	= _arg.XSCALE;
			altTdfm	= true;
			fprintf( stderr, "XSCALE=%g\n", _arg.XSCALE );
		}

		if( _arg.YSCALE != 999.0 ) {
			cYSCALE	= _arg.YSCALE;
			altTdfm	= true;
			fprintf( stderr, "YSCALE=%g\n", _arg.YSCALE );
		}

		if( _arg.SKEW != 999.0 ) {
			cSKEW	= _arg.SKEW;
			altTdfm	= true;
			fprintf( stderr, "SKEW=%g\n", _arg.SKEW );
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

	ctx.Tdfm.TPrint( stderr, "Tdfm=" );

	if( _arg.FLD ) {
		ctx.FLD = _arg.FLD;
		fprintf( stderr, "FLD=%c\n", _arg.FLD );
	}

	if( ctx.FLD == 'X' ) {
		ctx.FLD = (GBL.A.z == GBL.B.z ? 'N' : 'Y');
		fprintf( stderr, "FLD=%c (was X)\n", ctx.FLD );
	}

	if( _arg.MODE ) {
		ctx.MODE = _arg.MODE;
		fprintf( stderr, "MODE=%c\n", _arg.MODE );
	}

	if( ctx.MODE == 'Z' ) {
		ctx.MODE = 'C';
		arg.CTR = 0.0;
		fprintf( stderr, "MODE=C (was Z)\n" );
	}
	else if( ctx.MODE == 'M' ) {
		ctx.MODE = 'N';
		arg.CTR = 0.0;
		fprintf( stderr, "MODE=N (was M)\n" );
	}

	if( arg.CTR != 999.0 )
		fprintf( stderr, "CTR=%g\n", arg.CTR );

	fprintf( stderr, "\n" );

	return true;
}


