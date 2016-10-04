

#include	"CGBL_dmesh.h"
#include	"RegionToRegionMap.h"
#include	"CreateMesh.h"
#include	"ImproveMesh.h"

#include	"Maths.h"
#include	"Metrics.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

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
/* class CMatch -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Used for testing how well triangles match.
//
class CMatch {
public:
    vector<Point>	pts;	// points in image a
    vector<double>	a;		// value in image a (source)
    vector<double>	b;		// value in image b (target)
};

/* --------------------------------------------------------------- */
/* Metric -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static double Metric(
    const vector<Point>		&pts,
    const vector<double>	&av,
    const vector<double>	&bv,
    const char				*msg,
    FILE*					flog )
{
    if( GBL.mch.EMM ) {
        return	EarthMoversMetric( pts, av, bv,
                    GBL.mch.WDI, msg, flog );
    }
    else {
        return	FourierMatch( pts, av, bv, 25,
                    GBL.mch.WDI, msg, flog );
    }
}

/* --------------------------------------------------------------- */
/* SelfConsistent ------------------------------------------------ */
/* --------------------------------------------------------------- */

// All checking done at scaled (downsampled) sizes.
//
// Return non-zero if ok.
//
static int SelfConsistent( const vector<TAffine> &vT, FILE* flog )
{
    int	nT = vT.size();

    if( nT <= 1 )
        return nT;

// Compute value spread for each transform element

    double			ang_min = PI;
    double			ang_max = -PI;
    vector<double>	pmax( 6, -1e30 );
    vector<double>	pmin( 6,  1e30 );
    vector<double>	pspan( 6 );
    bool			extreme = false;

// Collect maxima and minima

    for( int i = 0; i < nT; ++i ) {

        double	ang =  vT[i].GetRadians();

        // handle atan() branch cuts

        if( i > 0 ) {

            if( ang > ang_max + PI )
                ang = ang - 2.0*PI;
            else if( ang < ang_min - PI )
                ang = ang + 2.0*PI;
        }

        ang_min = fmin( ang_min, ang );
        ang_max = fmax( ang_max, ang );

        for( int j = 0; j < 6; ++j ) {
            pmax[j] = max( pmax[j], vT[i].t[j] );
            pmin[j] = min( pmin[j], vT[i].t[j] );
        }
    }

    fprintf( flog,
    "\nAngle span: min, max, delta = %f %f %f\n",
    ang_min, ang_max, ang_max - ang_min );

// Calculate spans

    for( int j = 0; j < 6; ++j ) {

        pspan[j] = pmax[j] - pmin[j];

        if( j == 2 || j == 5 )	// translation
            extreme |= pspan[j] > GBL.mch.LDC;
        else					// rotation
            extreme |= pspan[j] > GBL.mch.LDR;
    }

// Any extreme spans?

    if( ang_max - ang_min > GBL.mch.LDA || extreme ) {

        fprintf( flog,
        "FAIL: Triangles too different:"
        " angles %f %f %f; LDA %f\n",
        ang_min, ang_max, ang_max-ang_min, GBL.mch.LDA );

        fprintf( flog,
        "FAIL: Triangles too different:"
        " elements %f %f %f %f %f %f; LDR %f LDC %f\n",
        pspan[0], pspan[1], pspan[2],
        pspan[3], pspan[4], pspan[5],
        GBL.mch.LDR, GBL.mch.LDC );

        nT = 0;
    }

    return nT;
}

/* --------------------------------------------------------------- */
/* ThmbConsistent ------------------------------------------------ */
/* --------------------------------------------------------------- */

static int ThmbConsistent(
    const vector<TAffine>	&vT,
    const TAffine			&T0,
    FILE*					flog )
{
// Average the final trans

    TAffine	A	= vT[0];
    int		nT	= vT.size();

    if( nT > 1 ) {

        for( int i = 1; i < nT; ++i ) {

            for( int j = 0; j < 6; ++j )
                A.t[j] += vT[i].t[j];
        }

        for( int j = 0; j < 6; ++j )
            A.t[j] /= nT;
    }

// Differences

    const int	f = 2;
    double		LDR = f * GBL.mch.LDR,
                LDC = f * GBL.mch.LDC;
    double		D[6];

    for( int i = 0; i < 6; ++i )
        D[i] = fabs( A.t[i] - T0.t[i] );

    fprintf( flog,
    "Mesh-Thm:%7.4f %7.4f %8.2f   %7.4f %7.4f %8.2f"
    " %dx(LDR,LDC)=(%.2f,%.2f) -- ",
    D[0], D[1], D[2], D[3], D[4], D[5],
    f, LDR, LDC );

// Check

    for( int i = 0; i < 6; ++i ) {

        if( i == 2 || i == 5 ) {
            if( D[i] > LDC )
                goto reject;
        }
        else if( D[i] > LDR )
            goto reject;
    }

    fprintf( flog, "OK.\n" );
    return true;

reject:
    fprintf( flog, "Reject.\n" );
    return false;
}

/* --------------------------------------------------------------- */
/* UserConsistent ------------------------------------------------ */
/* --------------------------------------------------------------- */

static int UserConsistent( const vector<TAffine> &vT, FILE* flog )
{
// Translations in bounds?

    int	nT = vT.size(), nXY = GBL.XYexp.size();

    if( nT && nXY ) {

        bool	allok = true;

        for( int i = 0; i < nT; ++i ) {

            bool	iok = false;

            // iok true if ANY bound satisfied

            for( int j = 0; j < nXY; ++j ) {

                Point	p( vT[i].t[2], vT[i].t[5] );

                iok |= p.Dist( GBL.XYexp[j] ) <= GBL.mch.DXY;
            }

            if( !iok ) {

                fprintf( flog,
                "Transform translation (%f %f)"
                " not on allowed list, tolerance %f\n",
                vT[i].t[2], vT[i].t[5], GBL.mch.DXY );
            }

            allok &= iok;
        }

        if( !allok ) {

            fprintf( flog,
            "FAIL: Transform not on allowed list.\n" );

            nT = 0;
        }
    }

    return nT;
}

/* --------------------------------------------------------------- */
/* ReportCenters ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReportCenters(
    const vector<TAffine>	&vT,
    const vector<Point>		&centers,
    int						napts,
    int						a_id,
    int						b_id,
    FILE					*flog )
{
    int	nc = centers.size();

    fprintf( flog, "\nFinal remapping of %d points.\n", napts );

    for( int i = 0; i < nc; ++i ) {

        const Point&	ca = centers[i];
        Point			cb = ca;

        vT[i].Transform( cb );

        fprintf( flog, "Center %f %f: ", ca.x, ca.y );
        vT[i].TPrint( flog );

        fprintf( flog,
        "Mapping region %d xy= %f %f to region %d xy= %f %f\n",
        a_id, ca.x, ca.y, b_id, cb.x, cb.y );
    }

    fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* WriteTrackEMTriangles ----------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTrackEMTriangles(
    vector<TAffine>			&vT,
    const vector<triangle>	&tri,
    const vector<vertex>	&ctl,
    FILE					*ftri )
{
    if( !ftri )
        return;

    int		ntri = tri.size();

    for( int j = 0; j < ntri; ++j ) {

        Point	p[3];

        for( int k = 0; k < 3; ++k ) {

            const vertex&	v = ctl[tri[j].v[k]];

            p[k] = Point( v.x, v.y );
        }

        fprintf( ftri,
        "%f %f %f %f %f %f       ",
        p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y );

        for( int k = 0; k < 3; ++k )
            vT[j].Transform( p[k] );

        fprintf( ftri,
        "%f %f %f %f %f %f\n",
        p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y );
    }
}

/* --------------------------------------------------------------- */
/* RegionToRegionMap --------------------------------------------- */
/* --------------------------------------------------------------- */

// Starting with approximate transform, find detailed correspondence.
// Parameter (maps) gets the list of found transforms and the values
// in image (ids) are indices (+10) into that list.
//
// If ftri != NULL, write found triangles there.
//
// Note
// ----
// All work and reporting is in scaled (downsampled) coords.
//
void RegionToRegionMap(
    ffmap				&maps,
    uint16				*ids,
    const PixPair		&px,
    const ConnRegion	&acr,
    const ConnRegion	&bcr,
    TAffine				tr_guess,
    FILE				*flog,
    FILE				*ftri )
{
    TAffine	T0		= tr_guess;
    int		w		= px.ws,
            h		= px.hs,
            sc		= px.scl,
            npix	= w * h,
            fullw	= px.wf,
            napts	= acr.pts.size();

    fprintf( flog, "\n---- Starting detailed region mapping ----\n" );

/* ---------------- */
/* Condition images */
/* ---------------- */

// Intended raster usage:
//
// - {av_vfy, acr.pts}	to QA results.
// - {bv_vfy, w, h}		to QA results.
//
// - {av_msh, ap_msh}	to build and refine the mesh.
// - {bv_msh, w, h}		to build and refine the mesh.
//
// We also create a bitmap IsIn[] with ones only inside the
// B connected region.

    const vector<double>&	av_vfy = *px.avs_vfy;
    const vector<double>&	bv_vfy = *px.bvs_vfy;
    const vector<double>&	av_aln = *px.avs_aln;
    vector<char>			IsIn( npix, 0 );

// Fill IsIn

    {
        int		nbpts = bcr.pts.size();

        for( int i = 0; i < nbpts; ++i )
            IsIn[int(bcr.pts[i].x) + w * int(bcr.pts[i].y)] = 1;
    }

// {av_msh, ap_msh} are subset of A mapping to IsIn

    vector<Point>	ap_msh;
    vector<double>	av_msh;

    for( int i = 0; i < napts; ++i ) {

        const Point&	ap = acr.pts[i];
        Point			bp = ap;

        tr_guess.Transform( bp );

        int ix = int(bp.x);
        int iy = int(bp.y);

        if( ix >= 0 && ix < w &&
            iy >= 0 && iy < h &&
            IsIn[ix+w*iy] ) {

            ap_msh.push_back( ap );
            av_msh.push_back( av_aln[int(ap.x) + w*int(ap.y)] );
        }
    }

// Region large enough?

    fprintf( flog, "Roughly %ld pixels map to B.\n",
    ap_msh.size() );

    if( ap_msh.size() < GBL.mch.MMA ) {

        fprintf( flog,
        "FAIL: Region too small - %ld pixels, MMA %d\n",
        ap_msh.size(), GBL.mch.MMA );

        return;
    }

// Normalize pixel values for ImproveMesh. The av_msh remain
// constant and can be normalized once now. ImproveMesh will
// select and normalize subsets of the bv_msh in GradDescStep.

    vector<double>	bv_msh = *px.bvs_aln;

    Normalize( av_msh );

/* ------------------------ */
/* Get common ap_msh bounds */
/* ------------------------ */

    IBox	B;

    MeshGetBounds( B, ap_msh, flog );

    if( B.R - B.L + 1 < 10 ||
        B.T - B.B + 1 < 10 ) {

        fprintf( flog,
        "FAIL: Region too slim - %d pixels, min %d\n",
        min( B.R - B.L, B.T - B.B ) + 1, 10 );

        return;
    }

/* -------------------------------------------- */
/* First create a single-triangle "affine" mesh */
/* -------------------------------------------- */

// We first try a single affine transformation for one large
// triangle on the region. If we can't get a good signal from
// that then a mesh of smaller triangles is unlikely to work.

    fprintf( flog, "\n---- Building mesh - affine ----\n" );

    vector<triangle>	tri;
    vector<vertex>		ctl;
    int					ntri = 1;

    MeshMakeSingleTri( tri, ctl, B, flog );

/* --------------------------- */
/* Run the mesh improver on it */
/* --------------------------- */

    vector<TAffine>	transforms;
    vector<Point>	centers;
    double			corr;

    corr = ImproveMesh(
            transforms, centers,
            tri, ctl,
            ap_msh, av_msh,
            bv_msh, w, h,
            tr_guess,
            GBL.ctx.RFA,
            flog, "affine" );

    if( !transforms.size() )
        return;
    else {

        TAffine	inv;
        inv.InverseOf( transforms[0] );

        transforms[0].TPrint(
            flog, "Best affine transform: " );

        inv.TPrint(
            flog, "    Inverse transform: " );

        tr_guess = transforms[0];
    }

/* ------------------------ */
/* Check affine mesh result */
/* ------------------------ */

    fprintf( flog, "\n---- QC - affine ----\n" );

    {

    // collect mapped points

        CMatch	mat;

        for( int k = 0; k < napts; ++k ) {

            const Point&	ap = acr.pts[k];
            Point			bp = ap;

            tr_guess.Transform( bp );

            if( bp.x >= 0 && bp.x < w-1 &&
                bp.y >= 0 && bp.y < h-1 &&
                IsIn[int(bp.x) + w*int(bp.y)] ) {

                double va = av_vfy[int(ap.x) + w*int(ap.y)];
                double vb = InterpolatePixel( bp.x, bp.y, bv_vfy, w );

                mat.pts.push_back( ap );
                mat.a.push_back( va );
                mat.b.push_back( vb );
            }
        }

        Normalize( mat.a );
        Normalize( mat.b );

    // get cross-corr, Fourier metric, yellow

        double	fm, y;
        double	sum	= 0.0;
        int		N	= mat.a.size();

        for( int i = 0; i < N; ++i )
            sum += mat.a[i] * mat.b[i];

        fm = Metric( mat.pts, mat.a, mat.b, "AFF", flog );

        y = PercentYellow( mat.a, mat.b, flog );

        fprintf( flog,
        "Affine Triangle: %d points, corr %f, fm %f, yellow %.2f\n",
        N, sum/N, fm, y );

    // Assess result quality

        if( GBL.mch.EMM ) {

            if( fm > 2.0*GBL.mch.EMT ) {

                fprintf( flog,
                "FAIL: Kicking out - Initial EMM too big (%f),"
                " 2*EMT %f\n", fm, 2.0*GBL.mch.EMT );

                return;
            }
        }
        else {	// use older Fourier Metric

            if( fm < GBL.mch.IFM ) {

                fprintf( flog,
                "FAIL: Kicking out - Initial FM too low (%f),"
                " IFM %f\n", fm, GBL.mch.IFM );

                for( int wvlen = 5; wvlen <= 40; wvlen += 5 ) {

                    double	ffm;

                    ffm = FourierMatch(
                            mat.pts, mat.a, mat.b, wvlen,
                            wvlen==5 || GBL.mch.WDI,
                            "AF2", flog );

                    fprintf( flog,
                    " wavelength %d, metric %f\n", wvlen, ffm );
                }

                return;
            }
        }
    }

    if( GBL.mch.ONE ) {
        fprintf( flog,
        "Cmd-line directive to use one affine only.\n" );
        goto quality_control;
    }

/* ----------------------------- */
/* OK, now try a deformable mesh */
/* ----------------------------- */

    fprintf( flog, "\n---- Building mesh - deformable ----\n" );

    if( MeshCreateX( tri, ctl, ap_msh, B, flog ) ) {

        fprintf( flog,
        "FAIL: Deformable triangular mesh failed - Small overlap?"
        " %ld pixels, MMA %d\n",
        ap_msh.size(), GBL.mch.MMA );

        return;
    }

    ntri = tri.size();

    if( !ntri ) {

        fprintf( flog,
        "FAIL: Deformable triangular mesh failed - No triangles."
        " %ld pixels, MMA %d\n",
        ap_msh.size(), GBL.mch.MMA );

        return;
    }

/* ----------------------- */
/* Improve deformable mesh */
/* ----------------------- */

    corr = ImproveMesh(
                transforms, centers,
                tri, ctl,
                ap_msh, av_msh,
                bv_msh, w, h,
                tr_guess,
                GBL.ctx.RFT,
                flog, "deformable mesh" );

    if( !transforms.size() )
        return;

/* --------------------------------- */
/* Generate quality control measures */
/* --------------------------------- */

quality_control:
    fprintf( flog, "\n---- QC - deformable ----\n" );

    vector<CMatch>	matches( ntri );	// each triangle
    CMatch			allp;				// over all triangles

    fprintf( flog, "Remapping %d points.\n", napts );

    for( int k = 0; k < napts; ++k ) {

        const Point&	ap = acr.pts[k];
        Point			bp = ap;
        int				t  = BestTriangle( tri, ctl, ap );

        transforms[t].Transform( bp );

        if( bp.x >= 0 && bp.x < w-1 &&
            bp.y >= 0 && bp.y < h-1 &&
            IsIn[int(bp.x) + w*int(bp.y)] ) {

            double va = av_vfy[int(ap.x) + w*int(ap.y)];
            double vb = InterpolatePixel( bp.x, bp.y, bv_vfy, w );

            matches[t].pts.push_back( ap );
            matches[t].a.push_back( va );
            matches[t].b.push_back( vb );

            allp.pts.push_back( ap );
            allp.a.push_back( va );
            allp.b.push_back( vb );
        }
    }

    Normalize( allp.a );
    Normalize( allp.b );

// Get Fourier metric over all triangles

    double	dfm;

    dfm = Metric( allp.pts, allp.a, allp.b, "DEF", flog );

    fprintf( flog, "All points, deformable, dfm %f\n\n", dfm );

// Sum over triangles: cross-corr, Fourier metric, yellow...

    double	weighted_sum	= 0.0;
    double	weighted_yellow	= 0.0;
    int		sum_pts			= 0;

    for( int k = 0; k < ntri; ++k ) {

        double	fm, y;
        double	sum	= 0.0;
        int		N	= matches[k].a.size();

        if( N ) {

            if( !Normalize( matches[k].a ) ||
                !Normalize( matches[k].b ) ) {

                fprintf( flog,
                "Triangle %d: Likely empty; skip metrics.\n\n", k );
                continue;
            }

            for( int i = 0; i < N; ++i )
                sum += matches[k].a[i] * matches[k].b[i];

            fm = Metric( matches[k].pts,
                    matches[k].a, matches[k].b, "TRI", flog );

            y = PercentYellow( matches[k].a, matches[k].b, flog );

            fprintf( flog,
            "Triangle %d: %d points, corr %f, fm %f, yellow %.2f\n\n",
            k, N, sum/N, fm, y );

            weighted_sum	+= fm * N;
            weighted_yellow	+= y * N;
            sum_pts			+= N;
        }
    }

// ...and get weighted average scores

    double	score	= weighted_sum / sum_pts;
    double	yell	= weighted_yellow / sum_pts;

/* --------------------- */
/* Assess result quality */
/* --------------------- */

    fprintf( flog, "\n---- Final reports ----\n" );

    char	buf[128];
    string	reason;

    if( yell < GBL.mch.FYL ) {
        sprintf( buf, " [Yellow=%g < FYL=%g]", yell, GBL.mch.FYL );
        reason += buf;
    }

    if( GBL.mch.EMM ) {

        fprintf( flog,
        "STAT: Overall %d points, corr %.4f,"
        " EMM %.4f, weighted EMM %.4f,"
        " cor+dfm %.4f, weighted yellow %6.4f\n",
        sum_pts, corr, dfm, score, corr+dfm, yell );

        if( score > GBL.mch.EMT ) {
            sprintf( buf,
            " [Weighted EMM=%g > EMT=%g]", score, GBL.mch.EMT );
            reason += buf;
        }

        if( dfm > GBL.mch.EMT*1.10 ) {
            sprintf( buf,
            " [EMM=%g > %g = (EMT=%g * 1.10)]",
            dfm, GBL.mch.EMT*1.10, GBL.mch.EMT );
            reason += buf;
        }

        //if( score > dfm*1.1 && score > GBL.mch.EMT/2.0 ) {
        //	reason +=	" [not great EMM, and Weighted EMM"
        //				" is higher by more than 10%]";
        //}
    }
    else {	// use older Fourier metric

        fprintf( flog,
        "STAT: Overall %d points, corr %.4f,"
        " dfm %.4f, weighted Fourier %.4f,"
        " cor+dfm %.4f, weighted yellow %6.4f\n",
        sum_pts, corr, dfm, score, corr+dfm, yell );

        if( dfm < GBL.mch.FFM )
            reason += " [Final Fourier too low]";

        if( sum_pts < 400000 && score < GBL.mch.FFM )
            reason += " [Small region and weighted metric < dfm]";

        if( score < 0.75*dfm )
            reason += " [Weighted less than 0.75 of dfm]";

        if( corr+dfm < GBL.mch.CPD )
            reason += " [Sum of corr+dfm too low]";
    }

// Failed if any reason cited

    if( reason.length() ) {

        fprintf( flog,
        "FAIL: Overall rejected:%s.\n", reason.c_str() );

        return;
    }

/* ------------------ */
/* Consistency checks */
/* ------------------ */

    if( !SelfConsistent( transforms, flog ) ||
        !ThmbConsistent( transforms, T0, flog ) ||
        !UserConsistent( transforms, flog ) ) {

        return;
    }

/* ---------------------- */
/* Success...make reports */
/* ---------------------- */

// Report centers in log

    ReportCenters( transforms, centers, napts, acr.id, bcr.id, flog );

// Print triangles {A} and Tr{A}, for trakEM

    WriteTrackEMTriangles( transforms, tri, ctl, ftri );

/* ------------------------------------ */
/* Collect results into data structures */
/* ------------------------------------ */

// Append maps entries

    int	next_id	= maps.transforms.size() + 10;

    for( int k = 0; k < ntri; ++k ) {

        maps.transforms.push_back( transforms[k] );
        maps.centers.push_back( centers[k] );
    }

// Paint ids with points that map a -> b

    for( int k = 0; k < napts; ++k ) {

        const Point&	ap = acr.pts[k];
        Point			bp = ap;
        int				t  = BestTriangle( tri, ctl, ap );

        transforms[t].Transform( bp );

        if( bp.x >= 0 && bp.x < w &&
            bp.y >= 0 && bp.y < h &&
            IsIn[int(bp.x) + w*int(bp.y)] ) {

            int ix = int(ap.x);
            int iy = int(ap.y);

            for( int x = 0; x < sc; ++x ) {

                for( int y = 0; y < sc; ++y )
                    ids[sc*ix+x + fullw*(sc*iy+y)] = next_id + t;
            }
        }
    }
}


