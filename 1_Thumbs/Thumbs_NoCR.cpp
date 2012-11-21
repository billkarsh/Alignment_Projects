

#include	"CGBL_Thumbs.h"
#include	"Thumbs.h"






/* --------------------------------------------------------------- */
/* Thumbs_NoCR --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Seek transform A->B composed as: Trigid * Tdfm * Tpretwk.
//
// The output is an entry in the ThmPair_lyrA_@_lyrB.txt file.
//
// Discussion
// ----------
// See notes in ApproximateMatch.cpp.
//
bool Thumbs_NoCR( const PixPair &px, FILE* flog )
{
/* ------------------- */
/* Handle bypass modes */
/* ------------------- */

	if( GBL.ctx.MODE == 'E' || GBL.ctx.MODE == 'F' )
		return true;

/* ------------------- */
/* Handle search modes */
/* ------------------- */

	CThmUtil	U( GBL.A, 1, GBL.B, 1, px,
					GBL.Tab, GBL.ctx.OLAP2D, flog );

	U.SetParams(
		GBL.ctx.HFANGDN, GBL.ctx.HFANGPR, GBL.ctx.RTRSH,
		GBL.ctx.OLAP1D, GBL.ctx.MODE, GBL.ctx.LIMXY );

	CThmScan	S;
	CorRec		best;

	S.Initialize( flog, best );
	S.SetTdfm( GBL.ctx.Tdfm );
	S.SetRThresh( GBL.ctx.RTRSH );
	S.SetNbMaxHt( GBL.ctx.NBMXHT );
	S.SetSweepConstXY( true );
	S.SetSweepPretweak( GBL.mch.PRETWEAK );
	S.SetUseCorrR( false );
	S.SetDisc( 0, 0, -1, -1 );

/* ----------------------- */
/* Create image thumbnails */
/* ----------------------- */

	OlapRec	olp;
	ThmRec	thm;
	int		nPriorAngles = U.SetStartingAngle( GBL.arg.CTR );

	U.Crop_NoCR( olp, GBL.ctx.XYCONF );

	if( !U.MakeThumbs( thm, olp, GBL.mch.THMDEC ) )
		return false;

/* ------ */
/* Search */
/* ------ */

	if( GBL.ctx.MODE == 'N' ) {

		if( !U.Disc( best, S, thm, olp, GBL.mch.PRETWEAK ) )
			return false;
	}
	else if( !U.Sweep( best, S, thm, nPriorAngles ) )
		return false;

/* ------ */
/* Finish */
/* ------ */

	return U.Finish( best, S, thm, olp, GBL.mch.TWEAKS );
}


