

#include	"lsq_MTrans.h"
#include	"lsq_MAffine.h"

#include	"TrakEM2_UTL.h"
#include	"File.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SetPointPairs(
    vector<LHSCol>	&LHS,
    vector<double>	&RHS,
    double			sc )
{
    int	nc	= vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        double	fz =
        (vRgn[C.r1].z == vRgn[C.r2].z ? same_strength : 1);

        double	x1 = C.p1.x * fz / sc,
                y1 = C.p1.y * fz / sc,
                x2 = C.p2.x * fz / sc,
                y2 = C.p2.y * fz / sc;
        int		j  = vRgn[C.r1].itr * NX,
                k  = vRgn[C.r2].itr * NX;

        // A1(p1) - A2(p2) = 0

        double	v[6]  = { x1,  y1,  fz, -x2, -y2, -fz};
        int		i1[6] = {  j, j+1, j+2,   k, k+1, k+2};
        int		i2[6] = {j+3, j+4, j+5, k+3, k+4, k+5};

        AddConstraint( LHS, RHS, 6, i1, v, 0.0 );
        AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
    }
}

/* --------------------------------------------------------------- */
/* SetIdentityTForm ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Explicitly set some TForm to Identity.
// @@@ Does it matter which one we use?
//
void MAffine::SetIdentityTForm(
    vector<LHSCol>	&LHS,
    vector<double>	&RHS,
    int				itr )
{
    double	stiff	= 1.0;

    double	one	= stiff;
    int		j	= itr * NX;

    AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
    AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
    AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
    AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
    AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
    AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;

// Don't do this with unite layer

    unite_layer = -1;

// Report which tile we set

    int	nr = vRgn.size();

    for( int k = 0; k < nr; ++k ) {

        if( vRgn[k].itr == itr ) {

            printf( "Ref region z=%d, id=%d\n",
            vRgn[k].z, vRgn[k].id );
            break;
        }
    }
}

/* --------------------------------------------------------------- */
/* SetUniteLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set one layer-full of TForms to those from a previous
// solution output file gArgs.unt_file.
//
void MAffine::SetUniteLayer(
    vector<LHSCol>	&LHS,
    vector<double>	&RHS,
    double			sc )
{
/* --------- */
/* Once only */
/* --------- */

    if( unite_layer < 0 )
        return;

/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

    map<MZIDR,TAffine>	M;

    LoadTAffineTbl_RngZ( M, unite_layer, unite_layer, unt_file );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

    double	stiff	= 0.001;

    int	nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	R = vRgn[i];

        if( R.z != unite_layer || R.itr < 0 )
            continue;

        map<MZIDR,TAffine>::iterator	it;

        it = M.find( MZIDR( R.z, R.id, R.rgn ) );

        if( it == M.end() )
            continue;

        double	one	= stiff,
                *t	= it->second.t;
        int		j	= R.itr * NX;

        AddConstraint( LHS, RHS, 1, &j, &one, one*t[0] );		j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[1] );		j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[2] / sc );	j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[3] );		j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[4] );		j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[5] / sc );	j++;
    }

/* --------- */
/* Once only */
/* --------- */

    unite_layer = -1;
}

/* --------------------------------------------------------------- */
/* SolveWithSquareness ------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SolveWithSquareness(
    vector<double>	&X,
    vector<LHSCol>	&LHS,
    vector<double>	&RHS,
    int				nTr )
{
/* -------------------------- */
/* Add squareness constraints */
/* -------------------------- */

    double	stiff = square_strength;

    for( int i = 0; i < nTr; ++i ) {

        int	j = i * NX;

        // equal cosines
        {
            double	V[2] = {stiff, -stiff};
            int		I[2] = {j, j+4};

            AddConstraint( LHS, RHS, 2, I, V, 0.0 );
        }

        // opposite sines
        {
            double	V[2] = {stiff, stiff};
            int		I[2] = {j+1, j+3};

            AddConstraint( LHS, RHS, 2, I, V, 0.0 );
        }
    }

/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

    WriteSolveRead( X, LHS, RHS, "A-Square", nproc, false );
    PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* SolveWithUnitMag ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Effectively, we want to constrain the cosines and sines
// so that c^2 + s^2 = 1. We can't make constraints that are
// non-linear in the variables X[], but we can construct an
// approximation using the {c,s = X[]} of the previous fit:
// c*x + s*y = 1. To reduce sensitivity to the sizes of the
// previous fit c,s, we normalize them by m = sqrt(c^2 + s^2).
//
void MAffine::SolveWithUnitMag(
    vector<double>	&X,
    vector<LHSCol>	&LHS,
    vector<double>	&RHS,
    int				nTR )
{
    double	stiff = scale_strength;

    for( int i = 0; i < nTR; ++i ) {

        int		j = i * NX;
        double	c = X[j];
        double	s = X[j+3];
        double	m = sqrt( c*c + s*s );

        // c*x/m + s*y/m = 1

        double	V[2] = {c * stiff, s * stiff};
        int		I[2] = {j, j+3};

        AddConstraint( LHS, RHS, 2, I, V, m * stiff );
    }

    WriteSolveRead( X, LHS, RHS, "A-Unimag", nproc, false );
    printf( "\t\t\t\t" );
    PrintMagnitude( X );
}

/* --------------------------------------------------------------- */
/* RescaleAll ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::RescaleAll(
    vector<double>	&X,
    double			sc )
{
    int	nr	= vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        X[itr+2] *= sc;
        X[itr+5] *= sc;
    }
}

/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::RotateAll(
    vector<double>	&X,
    double			degcw )
{
    TAffine	T, R;
    int		nr	= vRgn.size();

    R.SetCWRot( degcw, Point(0,0) );

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        TAffine	t( &X[itr] );

        T = R * t;
        T.CopyOut( &X[itr] );
    }
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::NewOriginAll(
    vector<double>	&X,
    double			xorg,
    double			yorg )
{
    int	nr	= vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        X[itr+2] -= xorg;
        X[itr+5] -= yorg;
    }
}

/* --------------------------------------------------------------- */
/* DevFromTrans -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to see how much translation terms of each affine
// have moved from the trans-only starting values. We just list
// all tiles with dev > XXX, but do nothing with that for now.
//
void MAffine::DevFromTrans(
    const vector<double>	&T,
    const vector<double>	&X )
{
    int	nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[i];

        if( I.itr < 0 )
            continue;

        const double	*J = &T[I.itr * 2],
                        *K = &X[I.itr * NX];
        double			dx = J[0] - K[2],
                        dy = J[1] - K[5];

        if( (dx = sqrt( dx*dx + dy*dy  )) > 200 ) {
            printf( "Dev: %d.%d-%d dr= %d\n",
            I.z, I.id, I.rgn, int(dx) );
        }
    }
}

/* --------------------------------------------------------------- */
/* DevFromPrior -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to see how much translation terms of each affine
// have moved from their prior starting values. We just list
// all tiles with dev > XXX, but do nothing with that for now.
//
void MAffine::DevFromPrior(
    const vector<double>	&A,
    const vector<double>	&X )
{
    int	nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[i];

        if( I.itr < 0 )
            continue;

        const double	*J = &A[I.itr * NX],
                        *K = &X[I.itr * NX];
        double			dx = J[2] - K[2],
                        dy = J[5] - K[5];

        if( (dx = sqrt( dx*dx + dy*dy  )) > 200 ) {
            printf( "Dev: %d.%d-%d dr= %d\n",
            I.z, I.id, I.rgn, int(dx) );
        }
    }
}

/* --------------------------------------------------------------- */
/* LoadAffTable -------------------------------------------------- */
/* --------------------------------------------------------------- */

// z0 and nz resp. get first layer and num layers in scaffold.
//
void MAffine::LoadAffTable(
    vector<double>	&X,
    int				&z0,
    int				&nz,
    int				nTr )
{
    X.resize( nTr * NX );
    z0	= -1;
    nz	= 0;

// Load table

    printf( "Aff: Loading existing table.\n" );

    map<MZIDR,TAffine>	M;
    set<int>			Z;

    LoadTAffineTbl_AllZ( M, Z, priorafftbl );

// Data range

    if( !Z.size() ) {
        printf( "Aff: No layers in scaffold.\n" );
        exit( 42 );
    }

    z0 = *Z.begin();
    nz = *Z.rbegin();

    printf( "Aff: Z range [%d %d].\n", z0, nz );
    nz = nz - z0 + 1;

// Fill into X

    printf( "Aff: Mapping prior solutions.\n" );

    int	nr = vRgn.size(), nmapped = 0;

    for( int i = 0; i < nr; ++i ) {

        const RGN&	R = vRgn[i];

        if( R.itr < 0 )
            continue;

        map<MZIDR,TAffine>::iterator	it;

        it = M.find( MZIDR( R.z, R.id, R.rgn ) );

        if( it != M.end() ) {

            memcpy( &X[R.itr*NX], it->second.t, NX*sizeof(double) );
            ++nmapped;
        }
        else {
            // mark as no solution
            X[R.itr * NX] = 999.0;
            printf( "No prior for %d.%d-%d\n", R.z, R.id, R.rgn );
        }
    }

// Done

    printf( "Aff: Mapped %d affines.\n\n", nmapped );
    fflush( stdout );
}

/* --------------------------------------------------------------- */
/* UntwistScaffold ----------------------------------------------- */
/* --------------------------------------------------------------- */

class RgdSums {
public:
    double	Xa, Ya, Xb, Yb, XaXb, YaYb, XaYb, YaXb;
    int		za, zb, N;
};


// The standard way to obtain an external scaffold for use with
// the -prior option is to use XMLGetTF on HiRes.xml; the result
// of aligning low-res strips from very good montages. However,
// better angle calculations can be made later, after the down
// correspondence points are formed. Adjustments are calculated
// here as rigid transforms T{theta, kx, ky) from layer a to b.
//
// Let c=cos(theta), s=sin(theta), Sum over all point-pairs:
//
// E = Sum[Xb - cXa + sYa - kx]^2 + [Yb - sXa - cYa - ky]^2
//
// The params u = {theta,kx,ky) are determined by dE/du = 0.
// If we use notation [arg] => Sum[argi] over point-pairs,
//
// kx = ([Xb] - c[Xa] + s[Ya]) / N
// ky = ([Yb] - s[Xa] - c[Ya]) / N
//
//				[YaXb] - [XaYb] + ([Xa][Yb] - [Ya][Xb]) / N
// tan(theta) = --------------------------------------------
//				([Xa][Xb] + [Ya][Yb]) / N - [XaXb] - [YaYb]
//
void MAffine::UntwistScaffold(
    vector<double>	&X,
    int				z0,
    int				nz )
{
/* ------------------------------------------ */
/* Size sum vector for all layers in scaffold */
/* ------------------------------------------ */

    if( nz < 2 )
        return;

    vector<RgdSums>	vS( nz );

    memset( &vS[0], 0, nz*sizeof(RgdSums) );

/* ---------------------------- */
/* Do sums over all point pairs */
/* ---------------------------- */

    int	nc = vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        const RGN	&Ra = vRgn[C.r1],
                    &Rb = vRgn[C.r2];

        if( Ra.z != Rb.z + 1 )
            continue;

        if( X[Ra.itr * NX] == 999.0 )
            continue;

        if( X[Rb.itr * NX] == 999.0 )
            continue;

        Point	pa = C.p1,
                pb = C.p2;

        L2GPoint( pa, X, Ra.itr );
        L2GPoint( pb, X, Rb.itr );

        RgdSums	&S = vS[Ra.z - z0];

        S.Xa += pa.x;
        S.Ya += pa.y;
        S.Xb += pb.x;
        S.Yb += pb.y;

        S.XaXb += pa.x * pb.x;
        S.YaYb += pa.y * pb.y;
        S.XaYb += pa.x * pb.y;
        S.YaXb += pa.y * pb.x;

        S.za = Ra.z;
        S.zb = Rb.z;
        ++S.N;
    }

/* --------------------------------------- */
/* Map cumulative rigids by affected layer */
/* --------------------------------------- */

    map<int,TAffine>	M;
    TAffine				Tprev;	// identity tform to start

    for( int i = 1; i < nz; ++i ) {

        RgdSums	&S = vS[i];

        if( S.N < 2 )
            continue;

        double
        theta = atan(
            (S.YaXb - S.XaYb + (S.Xa*S.Yb - S.Ya*S.Xb)/S.N) /
            ((S.Xa*S.Xb + S.Ya*S.Yb)/S.N - S.XaXb - S.YaYb)
        ),
        c  = cos( theta ),
        s  = sin( theta ),
        kx = (S.Xb - c*S.Xa + s*S.Ya) / S.N,
        ky = (S.Yb - s*S.Xa - c*S.Ya) / S.N;

        TAffine	T( c, -s, kx, s, c, ky );

        // propagate up the stack
        T = T * Tprev;
        Tprev = T;

        M[S.za] = T;
    }

/* ----- */
/* Apply */
/* ----- */

    TAffine	Tcache;
    int		zcache = -1, nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	R = vRgn[(*zs)[i].i];

        if( R.itr < 0 )
            continue;

        if( X[R.itr * NX] == 999.0 )
            continue;

        TAffine	T( &X[R.itr * NX] );

        if( R.z != zcache ) {

            zcache = R.z;

            map<int,TAffine>::iterator	it = M.find( zcache );

            if( it != M.end() )
                Tcache = it->second;
            else
                Tcache.NUSetOne();
        }

        T = Tcache * T;
        T.CopyOut( &X[R.itr * NX] );
    }
}

/* --------------------------------------------------------------- */
/* AffineFromFile ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::AffineFromFile( vector<double> &X, int nTr )
{
    double	sc		= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Aff: %d unknowns; %ld constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// Apply previous final results at highest level, once only

    SetUniteLayer( LHS, RHS, sc );

// Standard starting point

    SetPointPairs( LHS, RHS, sc );

// Load the Affines A

    vector<double>	A;
    int				z0, nz;

    LoadAffTable( A, z0, nz, nTr );
    UntwistScaffold( A, z0, nz );

// Relatively weighted: A(pi) = A(pj)

    int	nc = vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        // A(p1) = A(p2)
        if( A[vRgn[C.r2].itr * NX] != 999.0 ) {

            int		j  = vRgn[C.r1].itr * NX;
            double	x1 = C.p1.x * scaf_strength / sc,
                    y1 = C.p1.y * scaf_strength / sc,
                    x2,
                    y2;
            Point	g2 = C.p2;

            L2GPoint( g2, A, vRgn[C.r2].itr );
            x2 = g2.x * scaf_strength / sc;
            y2 = g2.y * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }

        // A(p2) = T(p1)
        if( A[vRgn[C.r1].itr * NX] != 999.0 ) {

            int		j  = vRgn[C.r2].itr * NX;
            double	x1 = C.p2.x * scaf_strength / sc,
                    y1 = C.p2.y * scaf_strength / sc,
                    x2,
                    y2;
            Point	g2 = C.p1;

            L2GPoint( g2, A, vRgn[C.r1].itr );
            x2 = g2.x * scaf_strength / sc;
            y2 = g2.y * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }
    }

// Solve

    //SolveWithSquareness( X, LHS, RHS, nTr );
    //SolveWithUnitMag( X, LHS, RHS, nTr );

    WriteSolveRead( X, LHS, RHS, "A-FrmFil", nproc, false );
    PrintMagnitude( X );

    RescaleAll( X, sc );

//	DevFromPrior( A, X );

    fflush( stdout );
}

/* --------------------------------------------------------------- */
/* AffineFromFile2 ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::AffineFromFile2( vector<double> &X, int nTr )
{
    double	sc		= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Aff: %d unknowns; %ld constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// Apply previous final results at highest level, once only

    SetUniteLayer( LHS, RHS, sc );

// Standard starting point

    SetPointPairs( LHS, RHS, sc );

// Get the Affines A

    vector<double>	A;
    int				z0 = (*zs).begin()->z,
                    nz = (*zs).rbegin()->z - z0 + 1;

    AffineFromFile( A, nTr );
    UntwistScaffold( A, z0, nz );

// Relatively weighted: A(pi) = A(pj)

    int	nc = vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        // A(p1) = A(p2)
        if( A[vRgn[C.r2].itr * NX] != 999.0 ) {

            int		j  = vRgn[C.r1].itr * NX;
            double	x1 = C.p1.x * scaf_strength / sc,
                    y1 = C.p1.y * scaf_strength / sc,
                    x2,
                    y2;
            Point	g2 = C.p2;

            L2GPoint( g2, A, vRgn[C.r2].itr );
            x2 = g2.x * scaf_strength / sc;
            y2 = g2.y * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }

        // A(p2) = T(p1)
        if( A[vRgn[C.r1].itr * NX] != 999.0 ) {

            int		j  = vRgn[C.r2].itr * NX;
            double	x1 = C.p2.x * scaf_strength / sc,
                    y1 = C.p2.y * scaf_strength / sc,
                    x2,
                    y2;
            Point	g2 = C.p1;

            L2GPoint( g2, A, vRgn[C.r1].itr );
            x2 = g2.x * scaf_strength / sc;
            y2 = g2.y * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }
    }

// Solve

    //SolveWithSquareness( X, LHS, RHS, nTr );
    //SolveWithUnitMag( X, LHS, RHS, nTr );

    WriteSolveRead( X, LHS, RHS, "A-FrmFil2", nproc, false );
    PrintMagnitude( X );

    RescaleAll( X, sc );

//	DevFromPrior( A, X );

    fflush( stdout );
}

/* --------------------------------------------------------------- */
/* AffineFromTransWt --------------------------------------------- */
/* --------------------------------------------------------------- */

// Preferred way to get montages:
// - Relate affines to each other A1(p1) = A2(p2).
// - Add scaffold relations A(pi) = T(pj) at reduced strength.
//
void MAffine::AffineFromTransWt( vector<double> &X, int nTr )
{
    double	sc		= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Aff: %d unknowns; %ld constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// Standard starting point

    SetPointPairs( LHS, RHS, sc );

// Get the pure translations T

    MTrans			M;
    vector<double>	T;

    M.SetModelParams( gW, gH, -1, -1, -1, -1,
        nproc, -1, NULL, NULL, zs );
    M.SolveSystem( T, nTr );

// Relatively weighted A(pi) = T(pj)

    int	nc = vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        // A(p1) = T(p2)
        {
            int		j  = vRgn[C.r1].itr * NX,
                    k  = vRgn[C.r2].itr * 2;
            double	x1 = C.p1.x * scaf_strength / sc,
                    y1 = C.p1.y * scaf_strength / sc,
                    x2 = (C.p2.x + T[k  ]) * scaf_strength / sc,
                    y2 = (C.p2.y + T[k+1]) * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }

        // A(p2) = T(p1)
        {
            int		j  = vRgn[C.r2].itr * NX,
                    k  = vRgn[C.r1].itr * 2;
            double	x1 = C.p2.x * scaf_strength / sc,
                    y1 = C.p2.y * scaf_strength / sc,
                    x2 = (C.p1.x + T[k  ]) * scaf_strength / sc,
                    y2 = (C.p1.y + T[k+1]) * scaf_strength / sc;

            double	v[3]	= {  x1,  y1, scaf_strength };
            int		i1[3]	= {   j, j+1, j+2 },
                    i2[3]	= { j+3, j+4, j+5 };

            AddConstraint( LHS, RHS, 3, i1, v, x2 );
            AddConstraint( LHS, RHS, 3, i2, v, y2 );
        }
    }

// Solve

    //SolveWithSquareness( X, LHS, RHS, nTr );
    //SolveWithUnitMag( X, LHS, RHS, nTr );

    WriteSolveRead( X, LHS, RHS, "A-WtT", nproc, false );
    PrintMagnitude( X );

    RescaleAll( X, sc );

    DevFromTrans( T, X );
}

/* --------------------------------------------------------------- */
/* SolveSystemStandard ------------------------------------------- */
/* --------------------------------------------------------------- */

// Build and solve system of linear equations.
//
// Note:
// All translational variables are scaled down by 'scale' so they
// are sized similarly to the sine/cosine variables. This is only
// to stabilize solver algorithm. We undo the scaling on exit.
//
void MAffine::SolveSystemStandard( vector<double> &X, int nTr )
{
    double	scale	= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Aff: %d unknowns; %ld constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

/* ------------------ */
/* Get rough solution */
/* ------------------ */

    SetPointPairs( LHS, RHS, scale );

    if( unite_layer < 0 )
        SetIdentityTForm( LHS, RHS, nTr / 2 );
    else
        SetUniteLayer( LHS, RHS, scale );

    SolveWithSquareness( X, LHS, RHS, nTr );

/* ----------------------------------------- */
/* Use solution to add torsional constraints */
/* ----------------------------------------- */

    SolveWithUnitMag( X, LHS, RHS, nTr );

/* --------------------------- */
/* Rescale translational terms */
/* --------------------------- */

    RescaleAll( X, scale );
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::SolveSystem( vector<double> &X, int nTr )
{
#if 0

// This older method works fine for small montages, especially
// if small in the y-dimension, but there will be some banana
// curvature visible. The newer scaffold method remains more
// square at any size, so is preferred.

    SolveSystemStandard( X, nTr );

#else

    if( priorafftbl )
//		AffineFromFile( X, nTr );
        AffineFromFile2( X, nTr );
    else
        AffineFromTransWt( X, nTr );

#endif
}

/* --------------------------------------------------------------- */
/* UpdateScaffold ------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::UpdateScaffold( vector<double> &X, int nTr )
{
    int	z0, nz;

    LoadAffTable( X, z0, nz, nTr );
    UntwistScaffold( X, z0, nz );

    double			xbnd, ybnd;
    vector<double>	lrbt;

    Bounds( xbnd, ybnd, X, lrbt, 0, NULL );
    WriteTransforms( X, false, NULL );

    WriteTrakEM( xbnd, ybnd, X, 0, 0, 0, 0 );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTransforms(
    const vector<double>	&X,
    int						bstrings,
    FILE					*FOUT )
{
    printf( "---- Write transforms ----\n" );

    FILE	*ft  = FileOpenOrDie( "TAffineTable.txt", "w" ),
            *fx  = FileOpenOrDie( "magnitude_outliers.txt", "w" );
    double	smin = 100.0,
            smax = 0.0,
            smag = 0.0;
    int		nr   = vRgn.size(), nTr = 0;

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        if( I.itr < 0 )
            continue;

        int	j = I.itr * NX;

        ++nTr;

        fprintf( ft, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
        I.z, I.id, I.rgn,
        X[j  ], X[j+1], X[j+2],
        X[j+3], X[j+4], X[j+5] );

        if( FOUT ) {

            if( !bstrings ) {

                fprintf( FOUT, "TAFFINE %d.%d-%d %f %f %f %f %f %f\n",
                I.z, I.id, I.rgn,
                X[j  ], X[j+1], X[j+2],
                X[j+3], X[j+4], X[j+5] );
            }
            else {
                fprintf( FOUT, "TRANSFORM '%s::%d' %f %f %f %f %f %f\n",
                I.GetName(), I.rgn,
                X[j  ], X[j+1], X[j+2],
                X[j+3], X[j+4], X[j+5] );
            }
        }

        double	mag = Magnitude( X, I.itr );

        smag += mag;
        smin  = fmin( smin, mag );
        smax  = fmax( smax, mag );

        if( mag < 0.9 ) {
            fprintf( fx, "Low mag %f @ %d.%d-%d\n",
                mag, I.z, I.id, I.rgn );
        }
        else if( mag > 1.1 ) {
            fprintf( fx, "Hi  mag %f @ %d.%d-%d\n",
                mag, I.z, I.id, I.rgn );
        }
    }

    fclose( fx );
    fclose( ft );

    printf(
    "Average magnitude=%f, min=%f, max=%f, max/min=%f\n\n",
    smag/nTr, smin, smax, smax/smin );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteTrakEM(
    double					xmax,
    double					ymax,
    const vector<double>	&X,
    double					trim,
    int						xml_type,
    int						xml_min,
    int						xml_max )
{
    FILE	*f = FileOpenOrDie( "MultLayAff.xml", "w" );

    int	oid = 3;

    fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

    TrakEM2WriteDTD( f );

    fprintf( f, "<trakem2>\n" );

    fprintf( f,
    "\t<project\n"
    "\t\tid=\"0\"\n"
    "\t\ttitle=\"Project\"\n"
    "\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
    "\t\tn_mipmap_threads=\"8\"\n"
    "\t/>\n" );

    fprintf( f,
    "\t<t2_layer_set\n"
    "\t\toid=\"%d\"\n"
    "\t\ttransform=\"matrix(1,0,0,1,0,0)\"\n"
    "\t\ttitle=\"Top level\"\n"
    "\t\tlayer_width=\"%.2f\"\n"
    "\t\tlayer_height=\"%.2f\"\n"
    "\t>\n",
    oid++, xmax, ymax );

    int	prev	= -1;	// will be previously written layer
    int	offset	= int(2 * trim + 0.5);
    int	nr		= vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        // skip unused tiles
        if( I.itr < 0 )
            continue;

        // changed layer
        if( (*zs)[i].z != prev ) {

            if( prev != -1 )
                fprintf( f, "\t\t</t2_layer>\n" );

            fprintf( f,
            "\t\t<t2_layer\n"
            "\t\t\toid=\"%d\"\n"
            "\t\t\tthickness=\"0\"\n"
            "\t\t\tz=\"%d\"\n"
            "\t\t>\n",
            oid++, (*zs)[i].z );

            prev = (*zs)[i].z;
        }

        const char	*path;
        char		title[128];
        DisplayStrings( title, path, I );

        // fix origin : undo trimming
        int		j = I.itr * NX;
        double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
        double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

        fprintf( f,
        "\t\t\t<t2_patch\n"
        "\t\t\t\toid=\"%d\"\n"
        "\t\t\t\twidth=\"%d\"\n"
        "\t\t\t\theight=\"%d\"\n"
        "\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
        "\t\t\t\ttitle=\"%s\"\n"
        "\t\t\t\ttype=\"%d\"\n"
        "\t\t\t\tfile_path=\"%s\"\n"
        "\t\t\t\to_width=\"%d\"\n"
        "\t\t\t\to_height=\"%d\"\n",
        oid++, gW - offset, gH - offset,
        X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
        title, xml_type, path, gW - offset, gH - offset );

        if( xml_min < xml_max ) {

            fprintf( f,
            "\t\t\t\tmin=\"%d\"\n"
            "\t\t\t\tmax=\"%d\"\n"
            "\t\t\t/>\n",
            xml_min, xml_max );
        }
        else
            fprintf( f, "\t\t\t/>\n" );
    }

    if( nr > 0 )
        fprintf( f, "\t\t</t2_layer>\n" );

    fprintf( f, "\t</t2_layer_set>\n" );
    fprintf( f, "</trakem2>\n" );
    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MAffine::WriteJython(
    const vector<double>	&X,
    double					trim,
    int						Ntr )
{
    FILE	*f = FileOpenOrDie( "JythonTransforms.txt", "w" );

    fprintf( f, "transforms = {\n" );

    int	nr = vRgn.size();

    for( int i = 0, itrf = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        // skip unused tiles
        if( I.itr < 0 )
            continue;

        ++itrf;

        const char	*path;
        DisplayStrings( NULL, path, I );

        // fix origin : undo trimming
        int		j = I.itr * NX;
        double	x_orig = X[j  ]*trim + X[j+1]*trim + X[j+2];
        double	y_orig = X[j+3]*trim + X[j+4]*trim + X[j+5];

        fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f]%s\n",
            path, X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
            (itrf == Ntr ? "" : ",") );
    }

    fprintf( f, "}\n" );
    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::G2LPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    TAffine	I, T( &X[itr * NX] );
    I.InverseOf( T );
    I.Transform( p );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MAffine::L2GPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    TAffine	T( &X[itr * NX] );
    T.Transform( p );
}


void MAffine::L2GPoint(
    vector<Point>			&p,
    const vector<double>	&X,
    int						itr )
{
    TAffine	T( &X[itr * NX] );
    T.Transform( p );
}


