

#include	"CGBL_dmesh.h"
#include	"ApproximateMatch.h"






/* --------------------------------------------------------------- */
/* ApproximateMatch_NoCR ----------------------------------------- */
/* --------------------------------------------------------------- */

// Seek transform A->B composed as: Trgd * Tdfm * Tpretwk.
//
// Return true and add one entry to guesses if successful.
//
// Also, output an entry in the ThmPair_lyrA^lyrB.txt file.
//
// Discussion
// ----------
// See notes in ApproximateMatch.cpp.
//
bool ApproximateMatch_NoCR(
    vector<TAffine>	&guesses,
    const PixPair	&px,
    FILE*			flog )
{
    CThmUtil	U( GBL.A, 1, GBL.B, 1, px,
                    GBL.Tab, GBL.ctx.OLAP2D, flog );

    U.SetParams(
        GBL.ctx.HFANGDN, GBL.ctx.HFANGPR, GBL.ctx.RTRSH,
        GBL.ctx.OLAP1D, GBL.ctx.MODE, GBL.ctx.LIMXY,
        GBL.mch.WTHMPR );

/* ------------------- */
/* Handle bypass modes */
/* ------------------- */

    if( GBL.ctx.MODE == 'E' )
        return U.Echo( guesses );

    if( GBL.ctx.MODE == 'F' )
        return U.FromLog( guesses );

/* ------------------- */
/* Handle search modes */
/* ------------------- */

    CThmScan	S;
    CorRec		best;

    S.Initialize( flog, best );
    S.SetTdfm( GBL.ctx.Tdfm );
    S.SetRThresh( GBL.ctx.RTRSH );
    S.SetNbMaxHt( GBL.ctx.NBMXHT );
    S.SetSweepConstXY( true );
    S.SetSweepPretweak( GBL.mch.PRETWEAK );
    S.SetUseCorrR( true );
    S.SetDisc( 0, 0, -1, -1 );

/* ----------------------- */
/* Create image thumbnails */
/* ----------------------- */

    OlapRec	olp;
    ThmRec	thm;
    int		nPriorAngles = U.SetStartingAngle(
                                GBL.ctx.Tdfm, GBL.arg.CTR );

    U.Crop_NoCR( olp, GBL.ctx.XYCONF );

    if( !U.MakeThumbs( thm, olp, GBL.ctx.THMDEC ) )
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

    if( U.Finish( best, S, thm, olp, GBL.mch.TWEAKS ) ) {

        guesses.push_back( best.T );
        return true;
    }

    return false;
}


