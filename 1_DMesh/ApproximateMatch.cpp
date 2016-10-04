

#include	"CGBL_dmesh.h"
#include	"ApproximateMatch.h"


/*
    Notes
    -----
    matchparams.txt cases out most parameters to allow separate
    treatment for same and cross layer.

    Here we seek a composite transform: Trgd * Tdfm * Tpretwk.

    File params PRETWEAK and TWEAKS determine whether Tpretwk
    is different from unity. PRETWEAK enables small amounts of
    distortion to be automatically added to boost signal iff
    no solution > RTRSH is found otherwise. TWEAKS is applied
    after a solution is found, to optimize the starting point
    for subsequent mesh phase.

    Tdfm is hierarchically set by separating it from the initial
    Tab, but can be overridden by file params {SCALE, XSCALE,
    YSCALE, SKEW}. These are overridden if provided on command
    line. All that is overridden if command line -Tdfm given.

    The primary job of this function, though, is to seek a
    refinement of the initial tform Tab. Several strategies
    (modes) are available to accommodate user's assessment
    of Tab initial accuracy.

    MODE=N (small neighborhood search) is used when the Tab are
    already roughly known from a prior alignment. Here, file param
    LIMXY is used as the search radius within the correlation image,
    and the simple maximum R peak criterion is used. If -CTR= is
    given on the command line then that angle is used instead of
    that from Tab. This forces greater angle consistency in a block.

    MODE=M (N + zero angle) is really the same as N mode with
    -CTR=0, but the command line -CTR=0 can be omitted for
    convenience. This is appropriate for montaging to force
    greater angle consistency.

    In all other modes where an angle sweep is used, correlation
    peak determination uses F mode. This convolves the correlation
    image with a peak enhancing (LOG + I) filter, and peak hunting
    requires that no F-value be higher than NBMXHT * F-peak within
    a specified guard zone of a putative F-peak.

    Sweeping in MODE=Y (dynamic angle hunting) is used when little
    is known about the angle to begin with. As best angles are found
    within a local block of tiles they are tabulated so that a median
    of these can be taken as the starting guess for subsequent tile
    pairs. At the start of a block when there are too few entries
    to get a consensus, a denovo guess is made from Tab or given
    -CTR= angle. The range of angles explored is cased out as
    HFANGDN for the denovo phase and HFANGPR for the prior angles
    phase, which can be set much narrower than the former.

    MODE=C (center angle search) is used when the approximate
    angle is already known but still benefits from a refining
    sweep centered upon that estimate. Here, the central angle
    can either be implicit in the Tab of the idb, overridden
    by the command line Tab, or given explicitly as -CTR=.

    MODE=Z (zero angle search) is really the same as C mode, but
    the angle used is identically zero and need not be specified
    as a parameter. This is appropriate for denovo montaging when
    the stage coordinates are of unknown quality. If the stage
    coordinates are trustworthy then N or M mode is preferred.

    If the inital Tab is believed close, then the overlapping areas
    of the images can be determined from Tab. The FFT operations
    are sped up linearly vs the image area. Moreover, correlation
    peak hunting can be foiled by noise, so making smaller areas
    increases robustness. We crop the images using Tab, but guided
    by file param XYCONF which is a confidence value in range [0,1].
    Zero signifies no confidence in Tab and the images are not
    cropped. If the value is one, the overlap rectangles derived
    from Tab are used without modification, while intermediate
    confidence values add some margin to the rectangles to avoid
    excluding the peak.

    All correlation work is done at reduced scale THMDEC (2^N).
    This increases speed and can squeeze out fine detail which
    is often less reliably correlated than larger features. The
    final result is double-checked at full scale and if consistent
    with THMDEC scale, full scale coordinates are used.

    In angle sweep modes, file param LIMXY (if non-zero) is used
    to constrain solutions to a disc of that size. It is often
    better to reject an uncertain result to avoid introducing
    stresses in later pipeline stages.
/*






/* --------------------------------------------------------------- */
/* ApproximateMatch ---------------------------------------------- */
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
bool ApproximateMatch(
    vector<TAffine>		&guesses,
    const PixPair		&px,
    const ConnRegion	&acr,
    const ConnRegion	&bcr,
    FILE*				flog )
{
    CThmUtil	U( GBL.A, acr.id, GBL.B, bcr.id, px,
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

    if( !U.Crop( olp, acr, bcr, GBL.ctx.XYCONF ) )
        return false;

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


