

#include	"lsq_MTrans.h"
#include	"lsq_MAffine.h"
#include	"lsq_MHmgphy.h"

#include	"EZThreads.h"
#include	"File.h"
#include	"Maths.h"
#include	"TrakEM2_UTL.h"
#include	"Timer.h"

#include	<math.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const double	SQRTOL		= sin( 15 * PI/180 );
const int		EDITDELAY	= 4;

/* --------------------------------------------------------------- */
/* SetUniteLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Set one layer-full of TForms to those from a previous
// solution output file gArgs.unt_file.
//
void MHmgphy::SetUniteLayer(
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

    map<MZIDR,THmgphy>	M;

    LoadTHmgphyTbl_RngZ( M, unite_layer, unite_layer, unt_file );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

    double	stiff	= 0.01;

    int	nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	R = vRgn[i];

        if( R.z != unite_layer || R.itr < 0 )
            continue;

        map<MZIDR,THmgphy>::iterator	it;

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
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[6] );		j++;
        AddConstraint( LHS, RHS, 1, &j, &one, one*t[7] );		j++;
    }

/* --------- */
/* Once only */
/* --------- */

    unite_layer = -1;
}

/* --------------------------------------------------------------- */
/* RescaleAll ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::RescaleAll(
    vector<double>	&X,
    double			sc )
{
    THmgphy	D, U;
    int		nr = vRgn.size();

    U.NUSetScl( sc );
    D.NUSetScl( 1/sc );

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        THmgphy	T( &X[itr] );

        T = U * (T * D);

        T.CopyOut( &X[itr] );
    }
}

/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::RotateAll(
    vector<double>	&X,
    double			degcw )
{
    THmgphy	T, R;
    int		nr	= vRgn.size();

    R.SetCWRot( degcw, Point(0,0) );

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        THmgphy	t( &X[itr] );

        T = R * t;
        T.CopyOut( &X[itr] );
    }
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::NewOriginAll(
    vector<double>	&X,
    double			xorg,
    double			yorg )
{
    THmgphy	M( 1,0,-xorg, 0,1,-yorg, 0,0 );
    int		nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        THmgphy	T( &X[itr] );

        T = M * T;

        T.CopyOut( &X[itr] );
    }
}

/* --------------------------------------------------------------- */
/* HmgphyFromAffine ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Deriving homographies by approximating H(pi) = A(pj) is the
// best way to create a set of simultaneous linear relations.
//
void MHmgphy::HmgphyFromAffine( vector<double> &X, int nTr )
{
    double	sc		= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Hmg: %d unknowns; %d constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// Apply previous final results at highest level, once only

    SetUniteLayer( LHS, RHS, sc );

// Get the Affines A

    MAffine			M;
    vector<double>	A;

    M.SetModelParams( gW, gH,
        same_strength,
        square_strength,
        scale_strength,
        scaf_strength,
        nproc,
        -1, NULL,
        priorafftbl, zs );

    M.SolveSystem( A, nTr );

// SetPointPairs: H(pi) = A(pj)

    double	fz	= 1.0;
    int		nc	= vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        // H(p1) = A(p2)
        {
            int		j  = vRgn[C.r1].itr * NX;
            double	x1 = C.p1.x * fz / sc,
                    y1 = C.p1.y * fz / sc,
                    x2,
                    y2;
            Point	g2 = C.p2;

            M.L2GPoint( g2, A, vRgn[C.r2].itr );
            x2 = g2.x / sc;
            y2 = g2.y / sc;

            double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
            int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
                    i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

            AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

            v[3] = -x1*y2;
            v[4] = -y1*y2;

            AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
        }

        // H(p2) = A(p1)
        {
            int		j  = vRgn[C.r2].itr * NX;
            double	x1 = C.p2.x * fz / sc,
                    y1 = C.p2.y * fz / sc,
                    x2,
                    y2;
            Point	g2 = C.p1;

            M.L2GPoint( g2, A, vRgn[C.r1].itr );
            x2 = g2.x / sc;
            y2 = g2.y / sc;

            double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
            int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
                    i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

            AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

            v[3] = -x1*y2;
            v[4] = -y1*y2;

            AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
        }
    }

// Solve

    WriteSolveRead( X, LHS, RHS, "H-FrmA", nproc, false );
    PrintMagnitude( X );

    RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* HmgphyFromHmgphy ---------------------------------------------- */
/* --------------------------------------------------------------- */

// A second pass deriving better H from previous H is worthwhile.
// More such passes gain nothing for the extra compute time.
//
void MHmgphy::HmgphyFromHmgphy( vector<double> &X, int nTr )
{
    double	sc		= 2 * max( gW, gH );
    int		nvars	= nTr * NX;

    printf( "Hmg: %d unknowns; %d constraints.\n",
        nvars, vAllC.size() );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// Apply previous final results at highest level, once only

    SetUniteLayer( LHS, RHS, sc );

// Get the Homographies A

    vector<double>	A;
    HmgphyFromAffine( A, nTr );

// SetPointPairs: H(pi) = A(pj)

    double	fz	= 1.0;
    int		nc	= vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        // H(p1) = A(p2)
        {
            int		j  = vRgn[C.r1].itr * NX;
            double	x1 = C.p1.x * fz / sc,
                    y1 = C.p1.y * fz / sc,
                    x2,
                    y2;
            Point	g2 = C.p2;

            L2GPoint( g2, A, vRgn[C.r2].itr );
            x2 = g2.x / sc;
            y2 = g2.y / sc;

            double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
            int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
                    i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

            AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

            v[3] = -x1*y2;
            v[4] = -y1*y2;

            AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
        }

        // H(p2) = A(p1)
        {
            int		j  = vRgn[C.r2].itr * NX;
            double	x1 = C.p2.x * fz / sc,
                    y1 = C.p2.y * fz / sc,
                    x2,
                    y2;
            Point	g2 = C.p1;

            L2GPoint( g2, A, vRgn[C.r1].itr );
            x2 = g2.x / sc;
            y2 = g2.y / sc;

            double	v[5]	= { x1, y1, fz, -x1*x2, -y1*x2 };
            int		i1[5]	= {   j, j+1, j+2, j+6, j+7 },
                    i2[5]	= { j+3, j+4, j+5, j+6, j+7 };

            AddConstraint( LHS, RHS, 5, i1, v, x2 * fz );

            v[3] = -x1*y2;
            v[4] = -y1*y2;

            AddConstraint( LHS, RHS, 5, i2, v, y2 * fz );
        }
    }

// Solve

    WriteSolveRead( X, LHS, RHS, "H-FrmH", nproc, false );
    PrintMagnitude( X );

    RescaleAll( X, sc );
}

/* --------------------------------------------------------------- */
/* WriteSideRatios ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Experiment to guage trapezoidism by reporting the ratio of
// each image's left vertical side over its right side.
//
// A B
// D C
//
void MHmgphy::WriteSideRatios( const vector<double> &X )
{
    FILE	*f	= FileOpenOrDie( "SideRatios.txt", "w" );
    int		nr	= vRgn.size();
    MeanStd	M[4];

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        if( I.itr < 0 )
            continue;

        int	j = I.itr * NX;

        const Til2Img	*m;
        THmgphy			T( &X[j] );
        Point			A( 0,    0 ), B( 2200,    0 ),
                        D( 0, 2200 ), C( 2200, 2200 );
        double			d;

        RGN::GetMeta( &m, NULL, I, I );

        T.Transform( A );
        T.Transform( B );
        T.Transform( C );
        T.Transform( D );

        d = D.Dist( A ) / C.Dist( B );

        B.x -= A.x;
        B.y -= A.y;
        C.x -= A.x;
        C.y -= A.y;
        D.x -= A.x;
        D.y -= A.y;
        A.x  = 0;
        A.y  = 0;

        fprintf( f,
        "%d\t%.12g\t%.12g\t%.12g"
        "\t\t%g\t%g\t\t%g\t%g\t\t%g\t%g\t\t%g\t%g\n",
        m->cam, d, X[j+6], X[j+7],
        A.x, A.y, B.x, B.y, C.x, C.y, D.x, D.y );

        M[m->cam].Element( d );
    }

    for( int i = 0; i < 4; ++i ) {

        double	ave, std;

        M[i].Stats( ave, std );

        printf( "{cam,L/R,std}: {%d,%.12g,%.12g}\n", i, ave, std );
    }

    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* Fill_myc ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static vector<vector<int> >	myc;

static void Fill_myc()
{
    myc.resize( vRgn.size() );

    int	nc = vAllC.size();

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        if( vRgn[C.r1].itr < 0 )
            continue;

        if( vRgn[C.r2].itr < 0 )
            continue;

        myc[C.r1].push_back( i );
        myc[C.r2].push_back( i );
    }
}

/* --------------------------------------------------------------- */
/* AFromTbl ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void AFromTbl( vector<double> &X, int nTr )
{
    X.resize( nTr * 6 );

    map<MZIDR,TAffine>	M;
    set<int>			Z;

//	LoadTAffineTbl_AllZ( M, Z, "../s_iter200_Scaf_A_2400/TAffineTable.txt" );
    LoadTAffineTbl_AllZ( M, Z, "../stack/TAffineTable.txt" );

    int	nr = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        RGN&	R = vRgn[i];

        if( R.itr < 0 )
            continue;

        map<MZIDR,TAffine>::iterator	it;

        it = M.find( MZIDR( R.z, R.id, R.rgn ) );

        if( it != M.end() )
            memcpy( &X[R.itr*6], it->second.t, 6*sizeof(double) );
        else
            R.itr = -1;
    }

    Fill_myc();
}

/* --------------------------------------------------------------- */
/* OnePassTH ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CThrdat {
public:
    int			r0,
                rlim;
    vector<int>	Rslo;
    vector<int>	Rkil;
};

static vector<CThrdat>	vthr;
static vector<double>	*Xout;
static vector<double>	*Xin;
static double			w;
static int				gpass;


static void* _OnePass_HFromA( void* ithr )
{
    CThrdat	&me = vthr[(int)(long)ithr];

    int	i1[5] = { 0, 1, 2, 6, 7 },
        i2[5] = { 3, 4, 5, 6, 7 },
        ndegen = 0;

    for( int i = me.r0; i < me.rlim; ++i ) {

        const RGN&	R = vRgn[i];

        if( R.itr < 0 )
            continue;

        int	nc = myc[i].size();

        if( nc < 4 )
            continue;

        double	*RHS = &(*Xout)[R.itr * 8];
        double	LHS[8*8];
        TAffine	Ta( &(*Xin)[R.itr * 6] );
        TAffine	Tb;
        int		lastb = -1;	// cache Tb

        Zero_Quick( LHS, RHS, 8 );

        for( int j = 0; j < nc; ++j ) {

            const Constraint&	C = vAllC[myc[i][j]];
            Point				A, B;

            // Mixing old and new solutions is related to
            // "successive over relaxation" methods in other
            // iterative solution schemes. Experimentally,
            // I like w = 0.9 (same layer), 0.9 (down).

            if( C.r1 == i ) {

                int	bitr = vRgn[C.r2].itr;

                if( bitr < 0 )
                    continue;

                if( C.r2 != lastb ) {
                    Tb.CopyIn( &(*Xin)[bitr * 6] );
                    lastb = C.r2;
                }
                Tb.Transform( B = C.p2 );
                Ta.Transform( A = C.p1 );
                B.x = w * B.x + (1 - w) * A.x;
                B.y = w * B.y + (1 - w) * A.y;
                A = C.p1;
            }
            else {

                int	bitr = vRgn[C.r1].itr;

                if( bitr < 0 )
                    continue;

                if( C.r1 != lastb ) {
                    Tb.CopyIn( &(*Xin)[bitr * 6] );
                    lastb = C.r1;
                }
                Tb.Transform( B = C.p1 );
                Ta.Transform( A = C.p2 );
                B.x = w * B.x + (1 - w) * A.x;
                B.y = w * B.y + (1 - w) * A.y;
                A = C.p2;
            }

            double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

            AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

            v[3] = -A.x*B.y;
            v[4] = -A.y*B.y;

            AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
        }

        ndegen += !Solve_Quick( LHS, RHS, 8 );
    }

    if( ndegen )
        printf( "HFromA degen %d\n", ndegen );

    return NULL;
}


static void HFromH_SLOnly( double *RHS, int i, int ithr )
{
    const RGN&	R = vRgn[i];

    int	i1[5] = { 0, 1, 2, 6, 7 },
        i2[5] = { 3, 4, 5, 6, 7 };

    int	nc = myc[i].size();

    double	LHS[8*8];
    THmgphy	Ta( &(*Xin)[R.itr * 8] );
    THmgphy	Tb;
    int		lastb	= -1,	// cache Tb
            nSLc	= 0;

    Zero_Quick( LHS, RHS, 8 );

    for( int j = 0; j < nc; ++j ) {

        const Constraint&	C = vAllC[myc[i][j]];
        Point				A, B;

        // Mixing old and new solutions is related to
        // "successive over relaxation" methods in other
        // iterative solution schemes. Experimentally,
        // I like w = 0.9 (same layer), 0.9 (down).

        if( C.r1 == i ) {

            if( vRgn[C.r2].z != R.z )
                continue;

            int	bitr = vRgn[C.r2].itr;

            if( bitr < 0 )
                continue;

            if( C.r2 != lastb ) {
                Tb.CopyIn( &(*Xin)[bitr * 8] );
                lastb = C.r2;
            }
            Tb.Transform( B = C.p2 );
            Ta.Transform( A = C.p1 );
            B.x = w * B.x + (1 - w) * A.x;
            B.y = w * B.y + (1 - w) * A.y;
            A = C.p1;
        }
        else {

            if( vRgn[C.r1].z != R.z )
                continue;

            int	bitr = vRgn[C.r1].itr;

            if( bitr < 0 )
                continue;

            if( C.r1 != lastb ) {
                Tb.CopyIn( &(*Xin)[bitr * 8] );
                lastb = C.r1;
            }
            Tb.Transform( B = C.p1 );
            Ta.Transform( A = C.p2 );
            B.x = w * B.x + (1 - w) * A.x;
            B.y = w * B.y + (1 - w) * A.y;
            A = C.p2;
        }

        ++nSLc;

        double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

        AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

        v[3] = -A.x*B.y;
        v[4] = -A.y*B.y;

        AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
    }

    if( nSLc < 4 ||
        !Solve_Quick( LHS, RHS, 8 ) ||
        THmgphy( RHS ).Squareness() > SQRTOL ) {

        vthr[ithr].Rkil.push_back( i );
    }
    else
        vthr[ithr].Rslo.push_back( i );
}


static void* _OnePass_HFromH( void* ithr )
{
    CThrdat	&me = vthr[(int)(long)ithr];

    int	i1[5] = { 0, 1, 2, 6, 7 },
        i2[5] = { 3, 4, 5, 6, 7 };

    for( int i = me.r0; i < me.rlim; ++i ) {

        const RGN&	R = vRgn[i];

        if( R.itr < 0 )
            continue;

        int	nc = myc[i].size();

        if( nc < 4 )
            continue;

        double	*RHS = &(*Xout)[R.itr * 8];
        double	LHS[8*8];
        THmgphy	Ta( &(*Xin)[R.itr * 8] );
        THmgphy	Tb;
        int		lastb	= -1;	// cache Tb

        Zero_Quick( LHS, RHS, 8 );

        for( int j = 0; j < nc; ++j ) {

            const Constraint&	C = vAllC[myc[i][j]];
            Point				A, B;

            // Mixing old and new solutions is related to
            // "successive over relaxation" methods in other
            // iterative solution schemes. Experimentally,
            // I like w = 0.9 (same layer), 0.9 (down).

            if( C.r1 == i ) {

                int	bitr = vRgn[C.r2].itr;

                if( bitr < 0 )
                    continue;

                if( C.r2 != lastb ) {
                    Tb.CopyIn( &(*Xin)[bitr * 8] );
                    lastb = C.r2;
                }
                Tb.Transform( B = C.p2 );
                Ta.Transform( A = C.p1 );
                B.x = w * B.x + (1 - w) * A.x;
                B.y = w * B.y + (1 - w) * A.y;
                A = C.p1;
            }
            else {

                int	bitr = vRgn[C.r1].itr;

                if( bitr < 0 )
                    continue;

                if( C.r1 != lastb ) {
                    Tb.CopyIn( &(*Xin)[bitr * 8] );
                    lastb = C.r1;
                }
                Tb.Transform( B = C.p1 );
                Ta.Transform( A = C.p2 );
                B.x = w * B.x + (1 - w) * A.x;
                B.y = w * B.y + (1 - w) * A.y;
                A = C.p2;
            }

            double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

            AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

            v[3] = -A.x*B.y;
            v[4] = -A.y*B.y;

            AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
        }

        if( gpass < EDITDELAY ) {
            Solve_Quick( LHS, RHS, 8 );
            continue;
        }

        if( !Solve_Quick( LHS, RHS, 8 ) ||
            THmgphy( RHS ).Squareness() > SQRTOL ) {

            HFromH_SLOnly( RHS, i, (int)(long)ithr );
        }
    }

    return NULL;
}


// Remove and remark bad RGNs and Constraints.
//
static void Remark()
{
    int	ncut = 0, nkil = 0;
    int	nt = vthr.size();

    // for each thread's lists...
    for( int it = 0; it < nt; ++it ) {

        vector<int>&	vslo = vthr[it].Rslo;
        vector<int>&	vkil = vthr[it].Rkil;
        int				nr = vslo.size();

        ncut += nr;

        // for the SLOnlys...
        for( int ir = 0; ir < nr; ++ir ) {

            int				islo = vslo[ir];
            vector<int>&	vmine = myc[islo];
            int				nc = vmine.size();

            // for each constraint...
            for( int ic = 0; ic < nc; ++ic ) {

                Constraint&	C = vAllC[vmine[ic]];

                if( vRgn[C.r1].z == vRgn[C.r2].z )
                    continue;

                // mark constraint
                C.used = false;

                // remove it from other myc list
                int	iother = (C.r1 == islo ? C.r2 : C.r1);
                vector<int>&	vother = myc[iother];
                int				mc = vother.size();

                for( int jc = 0; jc < mc; ++jc ) {

                    if( vother[jc] == vmine[ic] ) {

                        vother.erase( vother.begin()+jc );

                        // mark other for kill?
                        if( vother.size() < 4 )
                            vkil.push_back( iother );

                        break;
                    }
                }

                // and remove it from my list
                vmine.erase( vmine.begin()+ic );
                --ic;
                --nc;
            }

            // kill me?
            if( nc < 4 ) {
                vkil.push_back( islo );
                --ncut;
            }
        }

        // for the kills...
        for( int ir = 0; ir < vkil.size(); ++ir ) {

            int	ikil = vkil[ir];

            // already dead?
            if( vRgn[ikil].itr < 0 )
                continue;

            vector<int>&	vmine = myc[ikil];
            int				nc = vmine.size();

            // for each constraint...
            for( int ic = 0; ic < nc; ++ic ) {

                Constraint&	C = vAllC[vmine[ic]];

                // mark constraint
                C.used = false;

                // remove it from other myc list
                int	iother = (C.r1 == ikil ? C.r2 : C.r1);
                vector<int>&	vother = myc[iother];
                int				mc = vother.size();

                for( int jc = 0; jc < mc; ++jc ) {

                    if( vother[jc] == vmine[ic] ) {

                        vother.erase( vother.begin()+jc );

                        // mark other for kill?
                        if( vother.size() < 4 )
                            vkil.push_back( iother );

                        break;
                    }
                }
            }

            // kill my list
            vmine.clear();

            // mark me dead
            vRgn[ikil].itr = -1;
            ++nkil;
        }
    }

    if( ncut || nkil )
        printf( "Tiles [cut, killed] = [%d, %d].\n", ncut, nkil );
}


static void OnePassTH(
    EZThreadproc	passproc,
    vector<double>	&shXout,
    vector<double>	&shXin,
    int				nXout,
    int				nthr,
    double			shw )
{
    int	nr = vRgn.size(),	// regions total
        nb;					// regions per thread

    if( nthr > nr )
        nthr = nr;

    Xout	= &shXout;
    Xin		= &shXin;
    w		= shw;
    nb		= nr / nthr;

    vthr.resize( nthr );
    Xout->resize( nXout );

    vthr[0].r0		= 0;
    vthr[0].rlim	= nb;

    for( int i = 1; i < nthr; ++i ) {
        CThrdat	&C = vthr[i];
        C.r0	= vthr[i-1].rlim;
        C.rlim	= (i == nthr-1 ? nr : C.r0 + nb);
    }

    if( !EZThreads( passproc, nthr, 1, "passproc" ) )
        exit( 42 );

    Remark();
    vthr.clear();
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::SolveSystem( vector<double> &X, int nTr )
{
    fflush( stdout );
    clock_t	t0 = StartTiming();

#if 0

    HmgphyFromHmgphy( X, nTr );

#else

    vector<double>	Xin;
    AFromTbl( Xin, nTr );
    clock_t			t2 = StartTiming();

    OnePassTH( _OnePass_HFromA, X, Xin, nTr*8, 16, 0.9 );
    printf( "Just HfromA\n" );
    PrintMagnitude( X );
    fflush( stdout );

    for( gpass = 0; gpass < 2000; ++gpass ) {	// 1000 - 2000 good
        Xin = X;
        OnePassTH( _OnePass_HFromH, X, Xin, nTr*8, 16, 0.9 );
        printf( "Done pass %d\n", gpass + 1 ); fflush( stdout );
    }

    PrintMagnitude( X );
    StopTiming( stdout, "Iters", t2 );
    fflush( stdout );
    myc.clear();

#endif

    StopTiming( stdout, "SolveH", t0 );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::WriteTransforms(
    const vector<double>	&X,
    int						bstrings,
    FILE					*FOUT )
{
    printf( "---- Write transforms ----\n" );

    FILE	*ft  = FileOpenOrDie( "THmgphyTable.txt", "w" ),
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

        fprintf( ft,
        "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%.12g\t%.12g\n",
        I.z, I.id, I.rgn,
        X[j  ], X[j+1], X[j+2],
        X[j+3], X[j+4], X[j+5],
        X[j+6], X[j+7] );

        if( !bstrings ) {

            fprintf( FOUT,
            "THMGPHY %d.%d-%d %f %f %f %f %f %f %.12g %.12g\n",
            I.z, I.id, I.rgn,
            X[j  ], X[j+1], X[j+2],
            X[j+3], X[j+4], X[j+5],
            X[j+6], X[j+7] );
        }
        else {
            fprintf( FOUT,
            "THMGPHY '%s::%d' %f %f %f %f %f %f %.12g %.12g\n",
            I.GetName(), I.rgn,
            X[j  ], X[j+1], X[j+2],
            X[j+3], X[j+4], X[j+5],
            X[j+6], X[j+7] );
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

    WriteSideRatios( X );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* TopLeft ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void TopLeft(
    double			&top,
    double			&left,
    const THmgphy	&T,
    int				gW,
    int				gH,
    double			trim )
{
    vector<Point>	cnr( 4 );

    cnr[0] = Point(      trim,      trim );
    cnr[1] = Point( gW-1-trim,      trim );
    cnr[2] = Point( gW-1-trim, gH-1-trim );
    cnr[3] = Point(      trim, gH-1-trim );

    T.Transform( cnr );

    top  = BIGD;
    left = BIGD;

    for( int k = 0; k < 4; ++k ) {
        top  = fmin( top,  cnr[k].y );
        left = fmin( left, cnr[k].x );
    }
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MHmgphy::WriteTrakEM(
    double					xmax,
    double					ymax,
    const vector<double>	&X,
    double					trim,
    int						xml_type,
    int						xml_min,
    int						xml_max )
{
    FILE	*f = FileOpenOrDie( "MultLayHmgphy.xml", "w" );

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
        THmgphy	T( &X[j] );
        double	x_orig;
        double	y_orig;

        TopLeft( y_orig, x_orig, T, gW, gH, trim );

        fprintf( f,
        "\t\t\t<t2_patch\n"
        "\t\t\t\toid=\"%d\"\n"
        "\t\t\t\twidth=\"%d\"\n"
        "\t\t\t\theight=\"%d\"\n"
        "\t\t\t\ttransform=\"matrix(1,0,0,1,%f,%f)\"\n"
        "\t\t\t\ttitle=\"%s\"\n"
        "\t\t\t\ttype=\"%d\"\n"
        "\t\t\t\tfile_path=\"%s\"\n"
        "\t\t\t\to_width=\"%d\"\n"
        "\t\t\t\to_height=\"%d\"\n",
        oid++, gW - offset, gH - offset,
        x_orig, y_orig,
        title, xml_type, path, gW - offset, gH - offset );

        if( xml_min < xml_max ) {

            fprintf( f,
            "\t\t\t\tmin=\"%d\"\n"
            "\t\t\t\tmax=\"%d\"\n"
            "\t\t\t>\n",
            xml_min, xml_max );
        }
        else
            fprintf( f, "\t\t\t>\n" );

        fprintf( f,
        "\t\t\t<ict_transform"
        " class=\"mpicbg.trakem2.transform.HomographyModel2D\""
        " data=\"%f %f %f %f %f %f %.12g %.12g 1\"/>\n"
        "\t\t\t</t2_patch>\n",
        X[j  ], X[j+1], X[j+2], X[j+3],
        X[j+4], X[j+5], X[j+6], X[j+7] );
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

void MHmgphy::WriteJython(
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
        THmgphy	T( &X[j] );
        double	x_orig;
        double	y_orig;

        TopLeft( y_orig, x_orig, T, gW, gH, trim );

        fprintf( f,
        "\"%s\" : [%f, %f, %f, %f, %f, %f, %.12g, %.12g]%s\n",
        path,
        X[j+0], X[j+3], X[j+6],
        X[j+1], X[j+4], X[j+7],
        x_orig, y_orig,
        (itrf == Ntr ? "" : ",") );
    }

    fprintf( f, "}\n" );
    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MHmgphy::G2LPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    THmgphy	I, T( &X[itr * NX] );
    I.InverseOf( T );
    I.Transform( p );
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MHmgphy::L2GPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    THmgphy	T( &X[itr * NX] );
    T.Transform( p );
}


void MHmgphy::L2GPoint(
    vector<Point>			&p,
    const vector<double>	&X,
    int						itr )
{
    THmgphy	T( &X[itr * NX] );
    T.Transform( p );
}


