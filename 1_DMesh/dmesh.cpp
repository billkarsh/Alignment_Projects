

#include	"CGBL_dmesh.h"
#include	"dmesh.h"
#include	"ApproximateMatch.h"
#include	"RegionToRegionMap.h"

#include	"Disk.h"
#include	"Inspect.h"
#include	"LinEqu.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	FITAFF		0
#define	FITHMG		0
#define	FITDRAW		0
#define	FITTAB		0

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CStatus {
public:
    int	argn,
        brgn,
        thmok,
        ntri;
public:
    CStatus( int argn, int brgn )
    : argn(argn), brgn(brgn), thmok(0), ntri(0) {};
};


class Match {
public:
    double	weight;
    Point	pa, pb;
    int		ra, rb;	// subregions
public:
    Match(
        int				ra,
        int				rb,
        const Point&	pa,
        const Point&	pb )
    : weight(1.0), pa(pa), pb(pb), ra(ra), rb(rb) {};
};


class Matches {
private:
    vector<Match>	vM;
public:
    void Tabulate(
        const PixPair			&px,
        const vector<CStatus>	&vstat,
        const ffmap				&maps,
        const uint8*			fold_mask_a,
        const uint8*			fold_mask_b,
        int						wf,
        int						hf,
        FILE					*flog );
    void WriteAs_CPOINT2();
    void WriteAs_JSON();
    static void WriteEmpty_JSON();
private:
    static bool GetMutex( CMutex &M, const char *tag, const char **sud );
    static void JSON_head();
    static void JSON_tail();
    void JSON_ps();
    void JSON_qs();
    void JSON_ws();
    int index( int &i0, int &iLim, int ra, int rb );
    void FitAffine(
        const PixPair	&px,
        int				argn,
        int				brgn,
        FILE			*flog );
    void FitHmgphy(
        const PixPair	&px,
        int				argn,
        int				brgn,
        FILE			*flog );
};

/* --------------------------------------------------------------- */
/* Class Matches ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Catalogue good feature matches, grouped by (argn,brgn).
//
void Matches::Tabulate(
    const PixPair			&px,
    const vector<CStatus>	&vstat,
    const ffmap				&maps,
    const uint8*			fold_mask_a,
    const uint8*			fold_mask_b,
    int						wf,
    int						hf,
    FILE					*flog )
{
    fprintf( flog,
    "\n---- Tabulating point matches ----\n" );

    int	nstat = vstat.size(), tr0 = 0;

// For each region pair (vstat entry)...

    for( int istat = 0; istat < nstat; tr0 += vstat[istat++].ntri ) {

        int	trlim = vstat[istat].ntri;

        if( !trlim )
            continue;

        trlim += tr0;

        // Good points for this region pair

        vector<Match>	vGood;

        // For each A-triangle...

        for( int itr = tr0; itr < trlim; ++itr ) {

            Point	pa = maps.centers[itr],
                    pb = pa;
            int		ra, rb, ix, iy;

            // Lookup for A-point

            ix = int(pa.x);
            iy = int(pa.y);

            if( ix >= 0 && ix < wf && iy >= 0 && iy < hf ) {

                ra = fold_mask_a[ix + wf*iy];

                if( ra < 0 ) {

                    fprintf( flog,
                    "Tabulate: Rgn %d -> %d Tri %d:"
                    " A-centroid has bad mask value: mask=%d @ (%d, %d).\n",
                    vstat[istat].argn, vstat[istat].brgn, itr,
                    ra, ix, iy );

                    continue;
                }

                if( ra > 0 && ra != vstat[istat].argn ) {

                    fprintf( flog,
                    "Tabulate: Rgn %d -> %d Tri %d:"
                    " A-centroid not in stat block: mask=%d.\n",
                    vstat[istat].argn, vstat[istat].brgn, itr,
                    ra );

                    continue;
                }
            }
            else {

                fprintf( flog,
                "Tabulate: Rgn %d -> %d Tri %d:"
                " A-centroid out of A-image bounds: (%d, %d).\n",
                vstat[istat].argn, vstat[istat].brgn, itr,
                ix, iy );

                continue;
            }

            // Lookup for B-point

            maps.transforms[itr].Transform( pb );

            ix = int(pb.x);
            iy = int(pb.y);

            if( ix >= 0 && ix < wf && iy >= 0 && iy < hf ) {

                rb = fold_mask_b[ix + wf*iy];

                if( rb < 0 ) {

                    fprintf( flog,
                    "Tabulate: Rgn %d -> %d Tri %d:"
                    " B-centroid has bad mask value: mask=%d @ (%d, %d).\n",
                    vstat[istat].argn, vstat[istat].brgn, itr,
                    rb, ix, iy );

                    continue;
                }

                if( rb > 0 && rb != vstat[istat].brgn ) {

                    fprintf( flog,
                    "Tabulate: Rgn %d -> %d Tri %d:"
                    " B-centroid not in stat block: mask=%d.\n",
                    vstat[istat].argn, vstat[istat].brgn, itr,
                    rb );

                    continue;
                }
            }
            else {

                fprintf( flog,
                "Tabulate: Rgn %d -> %d Tri %d:"
                " B-centroid out of B-image bounds: (%d, %d).\n",
                vstat[istat].argn, vstat[istat].brgn, itr,
                ix, iy );

                continue;
            }

            // It's good; set it aside

            vGood.push_back( Match( ra, rb, pa, pb ) );
        }

        // Keep good points only if sufficiently many...
        // Although the mesh triangles have already passed sanity checks,
        // a bad mapping at this stage should still raise doubts about
        // whether the optimizer worked. We'll keep the good ones only
        // if a clear majority of matches look good.

        int	nG = vGood.size();

        if( nG >= 0.80 * vstat[istat].ntri )
            vM.insert( vM.end(), vGood.begin(), vGood.end() );
        else
            nG = 0;

        fprintf( flog,
        "Tabulate: Rgn %d -> %d:"
        " Keeping %d of %d point-pairs.\n",
        vstat[istat].argn, vstat[istat].brgn,
        nG, vstat[istat].ntri );

        if( !nG )
            continue;

        // Model the transforms obtained from point pairs

#if FITAFF
        FitAffine( px,
            vstat[istat].argn,
            vstat[istat].brgn, flog );
#endif

#if FITHMG
        FitHmgphy( px,
            vstat[istat].argn,
            vstat[istat].brgn, flog );
#endif
    }
}


void Matches::WriteAs_CPOINT2()
{
    int	np = vM.size();

    if( !np )
        return;

    const char	*sud;
    CMutex		M;

    if( GetMutex( M, "P", &sud ) ) {

        for( int i = 0; i < np; ++i ) {

            const Match	&m = vM[i];

            printf(
            "CPOINT2"
            " %d.%d-%d %f %f"
            " %d.%d-%d %f %f\n",
            GBL.A.z, GBL.A.id, m.ra, m.pa.x, m.pa.y,
            GBL.B.z, GBL.B.id, m.rb, m.pb.x, m.pb.y );
        }

        fflush( stdout );
    }

    M.Release();
}


// Pattern:
//
// {
//     "matches": {
//         "p": [[a-x-values],[a-y-values]],
//         "q": [[b-x-values],[b-y-values]],
//         "w": [weights]
//     }
// }
//
// Note: value lists are comma-separated.
//
void Matches::WriteAs_JSON()
{
    CMutex	M;

    if( GetMutex( M, "P", NULL ) ) {

        JSON_head();

        if( vM.size() ) {
            JSON_ps();
            JSON_qs();
            JSON_ws();
        }

        JSON_tail();
    }

    M.Release();
}


void Matches::WriteEmpty_JSON()
{
    CMutex	M;

    if( GetMutex( M, "P", NULL ) ) {

        JSON_head();
        JSON_tail();
    }

    M.Release();
}


bool Matches::GetMutex( CMutex &M, const char *tag, const char **sud )
{
    const char	*_sud;
    char		name[256];

    if( GBL.A.z < GBL.B.z )
        _sud = "up";
    else if( GBL.A.z == GBL.B.z )
        _sud = "same";
    else
        _sud = "down";

    sprintf( name, "%s_%d_%s", _sud, GBL.A.z, tag );

    if( sud )
        *sud = _sud;

    return M.Get( name );
}


void Matches::JSON_head()
{
    printf( "{\n" );
    printf( "    \"matches\": {\n" );
}


void Matches::JSON_tail()
{
    printf( "    }\n" );
    printf( "}\n" );

    fflush( stdout );
}


void Matches::JSON_ps()
{
    int	n = vM.size();

    printf( "        \"p\": [[" );

        printf( "%.4f", vM[0].pa.x );
        for( int i = 1; i < n; ++i )
            printf( ",%.4f", vM[i].pa.x );

    printf( "],[" );

        printf( "%.4f", vM[0].pa.y );
        for( int i = 1; i < n; ++i )
            printf( ",%.4f", vM[i].pa.y );

    printf( "]],\n" );
}


void Matches::JSON_qs()
{
    int	n = vM.size();

    printf( "        \"q\": [[" );

        printf( "%.4f", vM[0].pb.x );
        for( int i = 1; i < n; ++i )
            printf( ",%.4f", vM[i].pb.x );

    printf( "],[" );

        printf( "%.4f", vM[0].pb.y );
        for( int i = 1; i < n; ++i )
            printf( ",%.4f", vM[i].pb.y );

    printf( "]],\n" );
}


void Matches::JSON_ws()
{
    int	n = vM.size();

    printf( "        \"w\": [" );

        printf( "%.4g", vM[0].weight );
        for( int i = 1; i < n; ++i )
            printf( ",%.4g", vM[i].weight );

    printf( "]\n" );
}


// Get count and index range [i0,iLim) of vM entries matching (ra,rb).
//
int Matches::index( int &i0, int &iLim, int ra, int rb )
{
    i0		= -1;
    iLim	= -1;

    int	n = vM.size();

// seek start

    for( int i = 0; i < n; ++i ) {

        const Match	&m = vM[i];

        if( m.ra == ra && m.rb == rb ) {
            i0		= i;
            iLim	= i + 1;
            goto seek_lim;
        }
    }

    return 0;

// seek lim

seek_lim:
    while( iLim < n ) {

        const Match	&m = vM[iLim];

        if( m.ra == ra && m.rb == rb )
            ++iLim;
        else
            break;
    }

exit:
    return iLim - i0;
}


void Matches::FitAffine(
    const PixPair	&px,
    int				argn,
    int				brgn,
    FILE			*flog )
{
    int	i0, iLim, np = index( i0, iLim, argn, brgn );

    if( np < 3 ) {
        fprintf( flog,
        "Pipe: Too few points to fit affine [%d].\n", np );
        return;
    }

// Create system of normal equations

    double	RHS[6];
    double	LHS[6*6];
    int		i1[3] = { 0, 1, 2 },
            i2[3] = { 3, 4, 5 };

    Zero_Quick( LHS, RHS, 6 );

    for( int i = i0; i < iLim; ++i ) {

        const Point&	A = vM[i].pa;
        const Point&	B = vM[i].pb;

        double	v[3] = { A.x, A.y, 1.0 };

        AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
        AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
    }

// Solve

    fprintf( flog,
    "Pipe: Aff solver returns: %d\n", Solve_Quick( LHS, RHS, 6 ) );

    TAffine	T( &RHS[0] );

// Report

    T.TPrint( flog, "Pipe: FitAffine: " );

#if FITTAB
    {
        const char	*sud;
        CMutex		M;

        if( GetMutex( M, "A", &sud ) ) {

            char	name[256];
            sprintf( name, "aff.%s", sud );
            FILE *f = fopen( name, "a" );

            if( f ) {

                // Write entry

                fprintf( f,
                "AFFINE"
                " %d.%d-%d %d.%d-%d"
                " %f %f %f %f %f %f\n",
                GBL.A.z, GBL.A.id, argn,
                GBL.B.z, GBL.B.id, brgn,
                RHS[0], RHS[1], RHS[2],
                RHS[3], RHS[4], RHS[5] );

                fflush( f );
                fclose( f );
            }
        }

        M.Release();
    }
#endif

// RMS error

    double	E = 0;

    for( int i = i0; i < iLim; ++i ) {

        Point	a = vM[i].pa;

        T.Transform( a );

        double	err = a.DistSqr( vM[i].pb );

        E += err;
    }

    E = sqrt( E / np );

    fprintf( flog, "Pipe: FitAffineRMSerr: %g\n", E );

// Paint

#if FITDRAW
    YellowView( px, T, flog );
#endif
}


void Matches::FitHmgphy(
    const PixPair	&px,
    int				argn,
    int				brgn,
    FILE			*flog )
{
    int	i0, iLim, np = index( i0, iLim, argn, brgn );

    if( np < 4 ) {
        fprintf( flog,
        "Pipe: Too few points to fit homography [%d].\n", np );
        return;
    }

// Create system of normal equations

    double	RHS[8];
    double	LHS[8*8];
    int		i1[5] = { 0, 1, 2, 6, 7 },
            i2[5] = { 3, 4, 5, 6, 7 };

    Zero_Quick( LHS, RHS, 8 );

    for( int i = i0; i < iLim; ++i ) {

        const Point&	A = vM[i].pa;
        const Point&	B = vM[i].pb;

        double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

        AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

        v[3] = -A.x*B.y;
        v[4] = -A.y*B.y;

        AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
    }

// Solve

    fprintf( flog,
    "Pipe: Hmg solver returns: %d\n", Solve_Quick( LHS, RHS, 8 ) );

    THmgphy	T( &RHS[0] );

// Report

    T.TPrint( flog, "Pipe: FitHmgphy: " );

#if FITTAB
    {
        const char	*sud;
        CMutex		M;

        if( GetMutex( M, "H", &sud ) ) {

            char	name[256];
            sprintf( name, "hmg.%s", sud );
            FILE *f = fopen( name, "a" );

            if( f ) {

                // Write entry

                fprintf( f,
                "HMGPHY"
                " %d.%d-%d %d.%d-%d"
                " %f %f %f %f %f %f %.12g %.12g\n",
                GBL.A.z, GBL.A.id, argn,
                GBL.B.z, GBL.B.id, brgn,
                RHS[0], RHS[1], RHS[2],
                RHS[3], RHS[4], RHS[5],
                RHS[6], RHS[7] );

                fflush( f );
                fclose( f );
            }
        }

        M.Release();
    }
#endif

// RMS error

    double	E = 0;

    for( int i = i0; i < iLim; ++i ) {

        Point	a = vM[i].pa;

        T.Transform( a );

        double	err = a.DistSqr( vM[i].pb );

        E += err;
    }

    E = sqrt( E / np );

    fprintf( flog, "Pipe: FitHmgphyRMSerr: %g\n", E );

// Paint

#if FITDRAW
    YellowView( px, T, flog );
#endif
}

/* --------------------------------------------------------------- */
/* RoughMatch ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool RoughMatch(
    vector<TAffine>		&guesses,
    const PixPair		&px,
    CCropMask			&CM,
    const ConnRegion	&acr,
    const ConnRegion	&bcr,
    FILE*				flog )
{
    if( guesses.size() > 0 )
        return true;

    if( GBL.ctx.FLD == 'N'
        && !GBL.mch.PXRESMSK
        && !CM.IsFile( GBL.idb ) ) {

        // Call NoCR at most once. Possible states
        // are {0=never called, 1=failed, 2=success}.

        static vector<TAffine>	T;
        static int				state = 0;

        int	calledthistime = false;

        if( !state ) {
            state = 1 + ApproximateMatch_NoCR( T, px, flog );
            calledthistime = true;
        }

        if( state == 2 ) {

            if( !calledthistime ) {
                fprintf( flog, "\n---- Thumbnail matching ----\n" );
                T[0].TPrint( flog, "Reuse Approx: Best transform " );
            }

            guesses.push_back( T[0] );
            return true;
        }

        if( !calledthistime ) {
            fprintf( flog, "\n---- Thumbnail matching ----\n" );
            fprintf( flog, "FAIL: Approx: Case already failed.\n" );
        }

        return false;
    }
    else
        return ApproximateMatch( guesses, px, acr, bcr, flog );
}

/* --------------------------------------------------------------- */
/* UpscaleCoords ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpscaleCoords( ffmap &maps, int scale )
{
    int	nT;

    if( scale > 1 && (nT = maps.transforms.size()) ) {

        for( int i = 0; i < nT; ++i ) {

            maps.transforms[i].MulXY( scale );

            maps.centers[i].x *= scale;
            maps.centers[i].y *= scale;
        }
    }
}

/* --------------------------------------------------------------- */
/* PipelineDeformableMap ----------------------------------------- */
/* --------------------------------------------------------------- */

// The main routine the pipeline should call.
//
// Ntrans		- returned tform count
// tr_array		- array of these values
// map_mask		- <10=no map, else tform[pix-10]
// px			- a and b image pixels
// fold_mask_a	- 0=fold, 1,2,3...=region #
// fold_mask_b	- 0=fold, 1,2,3...=region #
// flog			- detailed output
//
void PipelineDeformableMap(
    int				&Ntrans,
    double*			&tr_array,
    uint16*			map_mask,
    const PixPair	&px,
    const uint8*	fold_mask_a,
    const uint8*	fold_mask_b,
    FILE*			flog )
{
    int	wf = px.wf, hf = px.hf;

/* --------------- */
/* Default results */
/* --------------- */

    Ntrans		= 0;
    tr_array	= NULL;

    memset( map_mask, 0, wf * hf * sizeof(uint16) );

/* --------------------------------- */
/* Create the connected region lists */
/* --------------------------------- */

// Note that the connected region lists are always
// at the reduced resolution, if this is used.

    fprintf( flog, "\n---- Connected regions ----\n" );

    vector<ConnRegion>	Acr, Bcr;
    CCropMask			CM;

    if( GBL.ctx.FLD == 'N'
        && !GBL.mch.PXRESMSK
        && !CM.IsFile( GBL.idb ) ) {

        fprintf( flog, "Forcing single connected region.\n" );

        ConnRgnForce1( Acr, px.ws, px.hs );
        Bcr = Acr;
    }
    else {

        ConnRgnsFromFoldMask( Acr, fold_mask_a,
            wf, hf, px.scl, uint32(0.9 * GBL.mch.MMA), flog );

        ConnRgnsFromFoldMask( Bcr, fold_mask_b,
            wf, hf, px.scl, uint32(0.9 * GBL.mch.MMA), flog );
    }

/* ----------------------------------------- */
/* Find mappings for each region-region pair */
/* ----------------------------------------- */

    vector<CStatus>	vstat;
    ffmap			maps;  // transforms and centers
    FILE			*ftri	= NULL;

    //ftri = fopen( "Triangles.txt", "w" );

    for( int i = 0; i < Acr.size(); ++i ) {

        for( int j = 0; j < Bcr.size(); ++j ) {

            fprintf( flog, "\n---- Begin A-%d to B-%d ----\n",
            Acr[i].id, Bcr[j].id );

            // start list with user's transform arguments

            CStatus			stat( Acr[i].id, Bcr[j].id );
            vector<TAffine>	guesses = GBL.Tmsh;

            if( RoughMatch( guesses,
                    px, CM, Acr[i], Bcr[j], flog ) ) {

                stat.thmok = true;

                // Try to get detailed mesh solution from each
                // guess {user + all returned from RoughMatch}.
                // The first to be successful (diff count > 0)
                // causes break.

                for( int k = 0; k < guesses.size(); ++k ) {

                    int	count = maps.transforms.size();

                    // Downscale coordinates
                    guesses[k].MulXY( 1.0 / px.scl );

                    RegionToRegionMap( maps, map_mask,
                        px, Acr[i], Bcr[j],
                        guesses[k], flog, ftri );

                    count = maps.transforms.size() - count;

                    stat.ntri = count;

                    if( count )
                        break;
                }
            }

            vstat.push_back( stat );
        }
    }

    //if( ftri )
    //	fclose( ftri );

/* ------------ */
/* Report total */
/* ------------ */

    Ntrans = maps.transforms.size();

    fprintf( flog,
    "Pipe: Returning %d total transforms.\n", Ntrans );

/* ---------------- */
/* Report by region */
/* ---------------- */

    fprintf( flog,
    "\n---- Summary Region-region results ----\n" );

    fprintf( flog,
    "Pipe: Table Headers {Az t r Bz t r thmok ntri}\n" );

    int	nstat = vstat.size();

    for( int i = 0; i < nstat; ++i ) {

        const CStatus& S = vstat[i];

        fprintf( flog,
        "FOUND: %5d %4d %3d  %5d %4d %3d %3d %3d\n",
        GBL.A.z, GBL.A.id, S.argn,
        GBL.B.z, GBL.B.id, S.brgn,
        S.thmok, S.ntri );
    }

/* ---------------- */
/* Feature matches  */
/* ---------------- */

    if( Ntrans ) {

        // Report results at full image size

        UpscaleCoords( maps, px.scl );

        Matches	AB;

        AB.Tabulate( px, vstat, maps,
            fold_mask_a, fold_mask_b, wf, hf, flog );

        if( GBL.arg.JSON )
            AB.WriteAs_JSON();
        else
            AB.WriteAs_CPOINT2();

        tr_array = (double*)malloc( Ntrans * 6 * sizeof(double) );

        for( int i = 0; i < Ntrans; ++i ) {

            maps.transforms[i].ToMatlab();
            maps.transforms[i].CopyOut( tr_array + i*6 );
        }

        fprintf( flog, "\n" );
    }
    else {

        if( GBL.arg.JSON )
            Matches::WriteEmpty_JSON();

        memset( map_mask, 0, wf * hf * sizeof(uint16) );
    }
}


