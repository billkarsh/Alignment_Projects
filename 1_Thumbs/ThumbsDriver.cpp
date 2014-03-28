

#include	"CGBL_Thumbs.h"
#include	"Thumbs.h"

#include	"ImageIO.h"
#include	"Timer.h"
#include	"Memory.h"


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	clock_t		t0 = StartTiming();

/* ------------------ */
/* Parse command line */
/* ------------------ */

	if( !GBL.SetCmdLine( argc, argv ) )
		return 42;

/* ---------- */
/* Get images */
/* ---------- */

	PixPair	px;

	if( !px.Load(
			GBL.A, GBL.B, GBL.idb,
			GBL.mch.PXLENS, GBL.mch.PXRESMSK, GBL.mch.PXBRO,
			GBL.mch.PXDOG, GBL.mch.PXDOG_R1, GBL.mch.PXDOG_R2,
			stdout, GBL.arg.Transpose ) ) {

		goto exit;
	}

/* ------------------- */
/* Scaling adjustments */
/* ------------------- */

	GBL.ctx.OLAP1D /=  px.scl;
	GBL.ctx.OLAP2D /= (px.scl * px.scl);

/* ---------------------------------- */
/* Use foldmasks according to context */
/* ---------------------------------- */

	if( GBL.ctx.FLD == 'N' )
		Thumbs_NoCR( px, stdout );
	else {

		/* --------------- */
		/* Load fold masks */
		/* --------------- */

		clock_t				t1 = StartTiming();
		uint8				*fold_mask_a, *fold_mask_b;
		vector<ConnRegion>	Acr, Bcr;

		printf( "\n---- Foldmaps ----\n" );

		// Note that the foldmasks are always at full resolution.

		{
			CCropMask	CM, *pCM = &CM;

			if( !CM.ReadIDB( GBL.idb ) )
				pCM = NULL;

			fold_mask_a =
				GetFoldMask(
					GBL.idb, GBL.A, GBL.arg.fma,
					px.resmska, pCM,
					px.wf, px.hf,
					false, GBL.arg.Transpose,
					GBL.arg.SingleFold );

			fold_mask_b =
				GetFoldMask(
					GBL.idb, GBL.B, GBL.arg.fmb,
					px.resmskb, pCM,
					px.wf, px.hf,
					false, GBL.arg.Transpose,
					GBL.arg.SingleFold );
		}

		/* ----------------------- */
		/* Convert to Conn regions */
		/* ----------------------- */

		// Note that the connected region lists are always
		// at the reduced resolution, if this is used.

		printf( "\n---- Connected regions ----\n" );

		ConnRgnsFromFoldMask( Acr, fold_mask_a,
			px.wf, px.hf, px.scl,
			GBL.ctx.OLAP2D, stdout );

		ConnRgnsFromFoldMask( Bcr, fold_mask_b,
			px.wf, px.hf, px.scl,
			GBL.ctx.OLAP2D, stdout );

		StopTiming( stdout, "Conn regions", t1 );

		/* --------------------- */
		/* Run the pair mappings */
		/* --------------------- */

		for( int i = 0; i < Acr.size(); ++i ) {

			for( int j = 0; j < Bcr.size(); ++j ) {

				printf( "\n---- Begin A-%d to B-%d ----\n",
				Acr[i].id, Bcr[j].id );

				Thumbs( px, Acr[i], Bcr[j], stdout );
			}
		}

		RasterFree( fold_mask_a );
		RasterFree( fold_mask_b );
	}

/* ---- */
/* Done */
/* ---- */

exit:
	StopTiming( stdout, "Total", t0 );
	VMStats( stdout );

	return 0;
}


