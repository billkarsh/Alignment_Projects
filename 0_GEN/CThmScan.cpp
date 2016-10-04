

#include	"CThmScan.h"
#include	"Correlation.h"
#include	"EZThreads.h"
#include	"Maths.h"
#include	"Timer.h"

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CThmScan*	ME;






/* --------------------------------------------------------------- */
/* BigEnough ----------------------------------------------------- */
/* --------------------------------------------------------------- */

bool BigEnough( int sx, int sy, void *a )
{
    return	sx >= ME->olap1D &&
            sy >= ME->olap1D &&
            (long)sx * sy > (long)a;
}

/* --------------------------------------------------------------- */
/* EnoughPoints -------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool EnoughPoints( int count1, int count2, void *a )
{
    return count1 > (long)a && count2 > (long)a;
}

/* --------------------------------------------------------------- */
/* NewXFromParabola ---------------------------------------------- */
/* --------------------------------------------------------------- */

static double NewXFromParabola(
    double	x1,
    double	d,
    double	y0,
    double	y1,
    double	y2 )
{
    double D = d * (y0 - y2) / (2 * (y0 + y2 - y1 - y1));

    if( fabs( D ) <= 2 * d )
        x1 += D;

    return x1;
}

/* --------------------------------------------------------------- */
/* _TCDDo1 ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::_TCDDo1( int ic )
{
    CorRec&	C = TCD.vC[ic];
    RFromAngle( C, C.A, *TCD.thm, ic );
}

/* --------------------------------------------------------------- */
/* _TCDGet ------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _TCDGet( void* icstart )
{
    int	nc = ME->TCD.vC.size();

    for( int j = (long)icstart; j < nc; j += ME->TCD.nthr )
        ME->_TCDDo1( j );

    return NULL;
}

/* --------------------------------------------------------------- */
/* TCDGet -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// (1) Caller must previously init public TCD fields.
// (2) Call TCDGet to calculate specified CorRecs.
//
void CThmScan::TCDGet( int nthr )
{
    ME = this;

    int	nc = TCD.vC.size();

    if( nthr > nc )
        nthr = nc;

    TCD.nthr = nthr;

    if( !EZThreads( _TCDGet, nthr, 24, "_TCDGet", flog ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* PTWInterp ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Use NewXFromParabola to interpolate a better NU-magnitude
// where the center guess is (x1,y1) from the sweep and the
// side points are +/1 d from that.
//
// Return xnew.
//
double CThmScan::PTWInterp(
    double	&ynew,
    int		sel,
    double	x1,
    double	y1,
    double	d,
    double	deg,
    ThmRec	&thm )
{
    double	y0, y2, xnew;
    PTWRec&	T = TCD.vT[0];
    CorRec	C;

    T.sel	= sel;
    T.a		= x1 - d;
    RFromAngle( C, deg, thm, 0 );
    y0		= C.R;

    T.a		= x1 + d;
    RFromAngle( C, deg, thm, 0 );
    y2		= C.R;

    xnew = NewXFromParabola( x1, d, y0, y1, y2 );

    T.a		= xnew;
    RFromAngle( C, deg, thm, 0 );
    ynew	= C.R;

    if( ynew < y1 ) {
        xnew = x1;
        ynew = y1;
    }

    fprintf( flog, "PTWIntrp %2d: %5.3f @ %.3f\n", sel, ynew, xnew );

    return xnew;
}

/* --------------------------------------------------------------- */
/* CThmScan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

CThmScan::CThmScan()
{
    ME			= this;
    newAngProc	= NULL;
    flog		= stdout;
    rthresh		= 0.30;
    nbmaxht		= 0.75;
    err			= errOK;
    swpConstXY	= true;
    swpPretweak	= true;
    swpNThreads	= 1;
    useCorrR	= false;
    Ox			= 0;
    Oy			= 0;
    Rx			= -1;
    Ry			= -1;
    olap1D		= 0;
}

/* --------------------------------------------------------------- */
/* CThmScan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CThmScan::Initialize( FILE* flog, CorRec &best )
{
    SetFlog( flog );

    fprintf( flog, "\n---- Thumbnail matching ----\n" );

    best.X	= -999.0;
    best.Y	= -999.0;
    best.A	= -999.0;
    best.R	= -999.0;
}

/* --------------------------------------------------------------- */
/* Pretweaks ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Several near-unity transforms are applied to A to boost initial
// correlation at the given center angle.
//
// Return true if any changes made.
//
bool CThmScan::Pretweaks( double bestR, double deg, ThmRec &thm )
{
    fprintf( flog, "Pretweaks start, best R=%.3f\n", bestR );

    clock_t	t0 = StartTiming();
    bool	anychange = false;

// We will compose a Tptwk by multiplying NU transforms. We examine
// 5 transform types vsel = {Scl, XScl, YScl, XSkw, YSkw}. Each type
// will enter the product at most once, which is tracked with vused.

    vector<int>	vsel( 5 );
    vector<int>	vused( 5, 0 );

    for( int i = tafnuScl; i <= tafnuYSkw; ++i )
        vsel[i] = i;

// To decide which transform type to use, we will do a magnitude
// sweep with each and see which type got the highest peak R. We
// then use interpolation to improve the winner. If using that
// type gives a better R than before, it goes into the product.
// Repeat with other unused types.

    const double	astep	= 0.02;
    const int		nstep	= 5;

    TCD.thm = &thm;

    for( int itype = 0; itype < 5; ++itype ) {

        // Gather sweep data

        TCD.vC.clear();
        TCD.vT.clear();

        for( int sel = 0; sel < 5; ++sel ) {

            if( vused[sel] )
                continue;

            double	abase = 0.0;

            if( sel >= tafnuScl && sel <= tafnuYScl )
                abase = 1.0;

            for( int j = -nstep; j <= nstep; ++j ) {
                TCD.vC.push_back( CorRec( deg ) );
                TCD.vT.push_back( PTWRec( sel, abase + j * astep ) );
            }
        }

        TCDGet( swpNThreads );

        // Print table and pick best

        double	rbest	= 0;
        int		nc		= TCD.vC.size(),
                cursel	= -1,
                ibest	= -1;

        for( int ic = 0; ic < nc; ++ic ) {

            const CorRec&	C = TCD.vC[ic];
            const PTWRec&	T = TCD.vT[ic];

            if( T.sel != cursel ) {

                if( cursel != -1 )
                    fprintf( flog, "\n" );

                fprintf( flog, "PTWSweep %2d:", cursel = T.sel );
            }

            fprintf( flog, " %5.3f", C.R );

            if( C.R > rbest ) {
                rbest = C.R;
                ibest = ic;
            }
        }

        fprintf( flog, "\n" );

        if( ibest == -1 )
            break;

        // Optimize best with interpolator

        PTWRec	Tbest = TCD.vT[ibest];

        Tbest.a = PTWInterp( rbest, Tbest.sel,
                    Tbest.a, rbest, astep, deg, thm );

        // Is it better than before?

        if( rbest > bestR ) {

            TAffine	M;

            fprintf( flog, "PTWUsing %2d\n", Tbest.sel );
            M.NUSelect( Tbest.sel, Tbest.a );
            Tptwk = M * Tptwk;
            vused[Tbest.sel]	= 1;
            bestR				= rbest;
            anychange			= true;
        }
        else
            break;
    }

    TCD.vT.clear();
    TCD.vC.clear();

    StopTiming( flog, "Pretweaks", t0 );

    Tptwk.TPrint( flog, "Approx: Pretweak " );

    return anychange;
}

/* --------------------------------------------------------------- */
/* RFromAngle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::RFromAngle(
    CorRec	&C,
    double	deg,
    ThmRec	&thm,
    int		iaux )
{
    vector<Point>	pts = thm.ap;
    TAffine			R, T( Tptwk );
    int				ox, oy, rx, ry;

    C.A = deg;
    R.NUSetRot( deg * PI/180.0 );

    if( iaux >= 0 && iaux < TCD.vT.size() ) {
        TAffine	M;
        M.NUSelect( TCD.vT[iaux].sel, TCD.vT[iaux].a );
        T = M * Tptwk;
    }

    C.T = R * (Tdfm * T);
    C.T.Apply_R_Part( pts );

    if( newAngProc )
        newAngProc( ox, oy, rx, ry, deg );
    else {
        ox = Ox; oy = Oy;
        rx = Rx; ry = Ry;
    }

    olap1D = thm.olap1D;

    C.R = (useCorrR ? CorrImagesS : CorrImagesF)(
        flog, false, C.X, C.Y,
        pts, thm.av, thm.bp, thm.bv,
        BigEnough, (void*)thm.reqArea,
        EnoughPoints, (void*)thm.reqArea,
        0.0, nbmaxht, ox, oy, rx, ry, thm.ftc );
}

/* --------------------------------------------------------------- */
/* DebugAngs ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::DebugAngs(
    int		alyr,
    int		atil,
    int		blyr,
    int		btil,
    double	center,
    double	hlfwid,
    double	step,
    ThmRec	&thm )
{
    char	file[256];

    sprintf( file, "angs_%d.%d^%d.%d.log", alyr, atil, blyr, btil );
    FILE	*f = fopen( file, "w" );

    fprintf( f, "Deg\tR\tX\tY\n" );

    TCD.thm = &thm;
    TCD.vC.clear();

    for( double a = center-hlfwid; a <= center+hlfwid; a += step )
        TCD.vC.push_back( CorRec( a ) );

    TCDGet( swpNThreads );

    int	nc = TCD.vC.size();

    for( int ic = 0; ic < nc; ++ic ) {

        const CorRec&	C = TCD.vC[ic];
        fprintf( f, "%.3f\t%.4f\t%.3f\t%.3f\n", C.A, C.R, C.X, C.Y );
    }

    TCD.vC.clear();
    fclose( f );
}

/* --------------------------------------------------------------- */
/* RecordAngle --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void RecordAngle(
    FILE			*f,
    const char		*label,
    const CorRec	&C )
{
    fprintf( f,
    "%s: R=%6.3f, A=%8.3f, X=%8.2f, Y=%8.2f\n",
    label, C.R, C.A, C.X, C.Y );
}

/* --------------------------------------------------------------- */
/* AngleScanMaxR ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Original version selecting best angle based on maximal R.
//
// If favorCenter true, sample more densely near center.
//
double CThmScan::AngleScanMaxR(
    CorRec	&best,
    double	center,
    double	hlfwid,
    double	step,
    bool	favorCenter,
    ThmRec	&thm )
{
    fprintf( flog,
    "AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f\n",
    center, hlfwid, step );

    clock_t	t0 = StartTiming();

    TCD.thm = &thm;
    TCD.vC.clear();
    TCD.vC.push_back( CorRec( center ) );

    if( favorCenter ) {

        double	cwide = 2.0 * step,
                cstep = step / 3.0;

        for( double a = center-cwide; a <= center+cwide; a += cstep )
            TCD.vC.push_back( CorRec( a ) );

        // symmetrically opposite if full circle

        if( hlfwid >= 180 ) {

            double	opc = center + 180;

            for( double a = opc-cwide; a <= opc+cwide; a += cstep )
                TCD.vC.push_back( CorRec( a ) );
        }
    }

    for( double a = center-hlfwid; a <= center+hlfwid; a += step )
        TCD.vC.push_back( CorRec( a ) );

    TCDGet( swpNThreads );

    RecordAngle( flog, "Center", best = TCD.vC[0] );

    int	nc = TCD.vC.size();

    for( int ic = 1; ic < nc; ++ic ) {

        const CorRec&	C = TCD.vC[ic];
        RecordAngle( flog, "  Scan", C );

        if( C.R > best.R )
            best = C;
    }

    TCD.vC.clear();

    RecordAngle( flog, "  Best", best );

    StopTiming( flog, "AngleScan", t0 );

    return best.R;
}

/* --------------------------------------------------------------- */
/* AngleScanConstXY ---------------------------------------------- */
/* --------------------------------------------------------------- */

// New version selecting best angle based on maximal R, but only if
// coordinates XY vary smoothly in a neighborhood about that angle.
// Smoothness measured as linear corr. coeff. from fitting X vs A.

class CSort_vC_Dec {
public:
    const vector<CorRec>	&vC;
public:
    CSort_vC_Dec( const vector<CorRec> &vC )
        : vC(vC) {};
    bool operator() ( int a, int b )
        {return vC[a].R > vC[b].R;};
};


double CThmScan::AngleScanConstXY(
    CorRec	&best,
    double	center,
    double	hlfwid,
    double	step,
    ThmRec	&thm )
{
    fprintf( flog,
    "AngleScan: center=%.3f, hlfwid=%.3f, step=%.3f\n",
    center, hlfwid, step );

//RFromAngle( best, center, thm );
//return;

    clock_t	t0 = StartTiming();

    best.X	= 0.0;
    best.Y	= 0.0;
    best.A	= 0.0;
    best.R	= 0.0;

// Sweep and collect

    TCD.thm = &thm;
    TCD.vC.clear();

    for( double a = center-hlfwid; a <= center+hlfwid; a += step )
        TCD.vC.push_back( CorRec( a ) );

    TCDGet( swpNThreads );

    int	nC = TCD.vC.size();

    for( int ic = 0; ic < nC; ++ic )
        RecordAngle( flog, "Scan", TCD.vC[ic] );

// Make indices sorted by decreasing R

    vector<int>	order( nC );

    for( int i = 0; i < nC; ++i )
        order[i] = i;

    CSort_vC_Dec	sorter( TCD.vC );

    sort( order.begin(), order.end(), sorter );

// Evaluate sweep
// Search through decreasing R.
// Take first result for which there are also +- m data points, and
// the X and Y coords vary smoothly over that range (r > rthresh).
// Also note that coords are rounded to nearest pixel, and that
// lincor values that are NAN or INF are generally ok: they arise
// from a zero variance in the coords, meaning that the point is
// stationary, and that's a good thing because it's not noise.

    const double	rthresh = 0.7;
    const			int m = 3;
    const			int M = 2*m + 1;

    for( int i = 0; i < nC; ++i ) {

        int	ic = order[i];

        if( ic < m )
            continue;

        if( nC-1 - ic < m )
            continue;

        vector<double>	A( M ), X( M ), Y( M );
        double			lincor;
        int				n = 0;

        for( int j = ic - m; j <= ic + m; ++j ) {

            A[n] = TCD.vC[j].A;
            X[n] = ROUND( TCD.vC[j].X );
            Y[n] = ROUND( TCD.vC[j].Y );
            ++n;
        }

        LineFit( NULL, NULL, &lincor, &A[0], &X[0], 0, n );
        fprintf( flog, "LCOR: A=%8.3f, RX=%6.3f\n", A[m], lincor );

        if( isfinite( lincor ) && fabs( lincor ) < rthresh )
            continue;

        LineFit( NULL, NULL, &lincor, &A[0], &Y[0], 0, n );
        fprintf( flog, "LCOR: A=%8.3f, RY=%6.3f\n", A[m], lincor );

        if( isfinite( lincor ) && fabs( lincor ) < rthresh )
            continue;

        best = TCD.vC[ic];
        break;
    }

    TCD.vC.clear();

// Report

    RecordAngle( flog, "Best", best );

    StopTiming( flog, "AngleScan", t0 );

    return best.R;
}

/* --------------------------------------------------------------- */
/* AngleScanSel -------------------------------------------------- */
/* --------------------------------------------------------------- */

double CThmScan::AngleScanSel(
    CorRec	&best,
    double	center,
    double	hlfwid,
    double	step,
    bool	favorCenter,
    ThmRec	&thm )
{
    if( swpConstXY )
        return AngleScanConstXY( best, center, hlfwid, step, thm );
    else {
        return AngleScanMaxR(
                best, center, hlfwid, step, favorCenter, thm );
    }
}

/* --------------------------------------------------------------- */
/* AngleScanWithTweaks ------------------------------------------- */
/* --------------------------------------------------------------- */

double CThmScan::AngleScanWithTweaks(
    CorRec	&best,
    double	center,
    double	hlfwid,
    double	step,
    ThmRec	&thm )
{
    if( AngleScanSel( best, center, hlfwid, step, true, thm )
        < rthresh ) {

        if( swpPretweak ) {

            if( Pretweaks( best.R,
                    (best.R > 0.0 ? best.A : center), thm ) ) {

                AngleScanSel( best, center, hlfwid, step, true, thm );
            }
        }
    }

    return best.R;
}

/* --------------------------------------------------------------- */
/* PeakHunt ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Bracket search for peak.
//
double CThmScan::PeakHunt( CorRec &best, double hlfwid, ThmRec &thm )
{
    CorRec	C,
            B0	= best;
    double	L	= best.A - hlfwid,
            R	= best.A + hlfwid,
            M;
    int		k	= 1;

    clock_t	t0 = StartTiming();

    RFromAngle( best, M = (L+R)/2.0, thm );

    while( R - L > 0.0001 ) {

        double	a;

        ++k;

        // move left up
        RFromAngle( C, a = (L+M)/2.0, thm );

        if( C.R >= best.R ) {
            R		= M;
            best	= C;
            M		= a;
        }
        else
            L		= a;

        // move right back
        RFromAngle( C, a = (M+R)/2.0, thm );

        if( C.R >= best.R ) {
            L		= M;
            best	= C;
            M		= a;
        }
        else
            R		= a;
    }

    if( B0.R > best.R )
        best = B0;

    fprintf( flog,
    "PeakHunt: Best: K=%d, R=%.3f, A=%.3f, X=%.3f, Y=%.3f\n",
    k, best.R, best.A, best.X, best.Y );

    StopTiming( flog, "PeakHunt", t0 );

    return best.R;
}

/* --------------------------------------------------------------- */
/* UsePriorAngles ------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CThmScan::UsePriorAngles(
    CorRec	&best,
    int		nprior,
    double	ang0,
    double	hfangpr,
    ThmRec	&thm )
{
    fprintf( flog, "Approx: Using prior angles n=%d, med=%f\n",
    nprior, ang0 );

    if( AngleScanWithTweaks( best, ang0, hfangpr, 0.1, thm )
            < rthresh ||
        PeakHunt( best, 0.1 * 1.5, thm )
            < rthresh ) {

        fprintf( flog,
        "FAIL: Approx: Prior angles R=%g below thresh=%g\n",
        best.R, rthresh );

        err = errLowRPrior;
        return false;
    }

    return true;
}

/* --------------------------------------------------------------- */
/* DenovoBestAngle ----------------------------------------------- */
/* --------------------------------------------------------------- */

bool CThmScan::DenovoBestAngle(
    CorRec	&best,
    double	ang0,
    double	hfangdn,
    double	step,
    ThmRec	&thm,
    bool	failmsg )
{
    if( AngleScanWithTweaks( best, ang0, hfangdn, step, thm )
            < rthresh ||
        AngleScanSel( best, best.A, step*2.0, step*0.05, false, thm )
            < rthresh ||
        PeakHunt( best, step * 0.05 * 1.5, thm )
            < rthresh ) {

        if( failmsg ) {
            fprintf( flog,
            "FAIL: Approx: Denovo R=%g below thresh=%g\n",
            best.R, rthresh );
        }

        err = errLowRDenov;
        return false;
    }

    return true;
}

/* --------------------------------------------------------------- */
/* PostTweaks ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CThmScan::PostTweaks( CorRec &best, ThmRec &thm )
{
    SetUseCorrR( true );
    SetDisc( (int)best.X, (int)best.Y, 40, 40 );
    Pretweaks( best.R, best.A, thm );
    RFromAngle( best, best.A, thm );
}

/* --------------------------------------------------------------- */
/* FinishAtFullRes ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Get final offsets at full resolution.
//
void CThmScan::FinishAtFullRes( CorRec &best, ThmRec &thm )
{
    CorRec	b0 = best;
    int		ok;

    fprintf( flog,
    "Approx: LowRes  R=%.3f, X=%.3f, Y=%.3f\n",
    best.R, best.X, best.Y );

    clock_t	t0 = StartTiming();

    SetUseCorrR( true );
    SetDisc( (int)best.X, (int)best.Y, 80, 80 );
    RFromAngle( best, best.A, thm );

    ok = (fabs( best.X - b0.X ) <= 20)
      && (fabs( best.Y - b0.Y ) <= 20);

    fprintf( flog,
    "Approx: FullRes R=%.3f, X=%.3f, Y=%.3f, use=%c.\n",
    best.R, best.X, best.Y, (ok ? 'Y' : 'N') );

// Always report the thumb-R instead of the fullres value because
// fullres values are scaled differently. We prefer all values in
// the file to be comparable.

    if( ok )
        best.R = b0.R;
    else
        best = b0;

    StopTiming( flog, "FinishAtFullRes", t0 );
}


