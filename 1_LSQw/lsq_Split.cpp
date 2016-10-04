

#include	"lsq_Globals.h"
#include	"lsq_Split.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"

#include	<string.h>

#include	<stack>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// colors per layer
#define	PERZ	10000

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class GDrty {
// Has this worker dirtied his L/R?
public:
    uint8	L, R;
public:
    GDrty() : L(1), R(1) {};
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static Split*			ME;
static const char		*gpath;
static vector<GDrty>	gdrty;
static vector<uint8>	drtyS, drtyD;
static vector<int>		vzd;
static int				nthr,
                        saveclr,
                        maxspan,
                        zs;






/* --------------------------------------------------------------- */
/* Resize -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Allocate space for K vector.
//
void Split::Resize()
{
    int	nz = zohi - zolo + 1;

    K.resize( nz );

    for( int iz = 0; iz < nz; ++iz )
        K[iz].resize( vR[iz].nr, 0 );
}

/* --------------------------------------------------------------- */
/* _ColorMontages ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void* _ColorMontages( void* ithr )
{
    for( int iz = zilo + (long)ithr; iz <= zihi; iz += nthr ) {

        const Rgns&		R = vR[iz];
        vector<int>&	k = ME->K[iz];
        int				clr = (R.z ? PERZ * R.z : 1) - 1;

        for( int ir = 0; ir < R.nr; ++ir ) {

            // Valid tile?
            if( !FLAG_ISUSED( R.flag[ir] ) )
                continue;

            // Already assigned?
            if( k[ir] )
                continue;

            // New seed; new stack; new color
            stack<int>	s;
            s.push( ir );
            ++clr;

            while( !s.empty() ) {

                // Process this guy

                int	jr = s.top();
                s.pop();

                if( k[jr] )
                    continue;

                k[jr] = clr;

                // Push neighbors that are:
                // same-layer, virgin, valid.

                const vector<int>&	P  = R.pts[jr];
                int					np = P.size(), prev = -1;

                for( int ip = 0; ip < np; ++ip ) {

                    const CorrPnt&	C = vC[P[ip]];

                    if( !C.used )
                        continue;

                    if( C.z1 != iz || C.z2 != iz )
                        continue;

                    int other = (C.i1 == jr ? C.i2 : C.i1);

                    if( other == prev )
                        continue;

                    prev = other;

                    if( !k[other] && FLAG_ISUSED( R.flag[other] ) )
                        s.push( other );
                }
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* ColorMontages ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Walk pts to colorize islands per montage.
//
void Split::ColorMontages()
{
    int	nz = zihi - zilo + 1;

    ME		= (Split*)this;
    nthr	= maxthreads;

    if( nthr > nz )
        nthr = nz;

    if( !EZThreads( _ColorMontages, nthr, 1, "_ColorMontages" ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* GUpdt --------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::GUpdt()
{
// Send my GDrty to master; then get global vector

    if( wkid > 0 ) {

        // send mine
        MPISend( &gdrty[wkid], sizeof(GDrty), 0, wkid );

        // get all
        MPIRecv( &gdrty[0], nwks * sizeof(GDrty), 0, wkid );
    }
    else if( nwks > 1 ) {

        // get each worker
        for( int iw = 1; iw < nwks; ++iw )
            MPIRecv( &gdrty[iw], sizeof(GDrty), iw, iw );

        // send each worker
        for( int iw = 1; iw < nwks; ++iw )
            MPISend( &gdrty[0], nwks * sizeof(GDrty), iw, iw );
    }
}

/* --------------------------------------------------------------- */
/* KSend --------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::KSend( int zlo, int zhi, int toLorR )
{
    void	*buf;
    int		bytes,
            wdst;

    if( toLorR == 'L' ) {

        if( !gdrty[wkid].L )
            return true;

        wdst = wkid - 1;

        if( wdst < 0 )
            return true;
    }
    else {

        if( !gdrty[wkid].R )
            return true;

        wdst = wkid + 1;

        if( wdst >= nwks )
            return true;
    }

    for( int iz = zlo; iz <= zhi; ++iz ) {

        buf		= (void*)&K[iz][0];
        bytes	= sizeof(int) * vR[iz].nr;

        MPISend( buf, bytes, wdst, iz - zlo );
    }

    return true;
}

/* --------------------------------------------------------------- */
/* KRecv --------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::KRecv( int zlo, int zhi, int fmLorR )
{
    void	*buf;
    int		bytes,
            wsrc;

    if( fmLorR == 'L' ) {

        wsrc = wkid - 1;

        if( wsrc < 0 )
            return true;

        if( !gdrty[wsrc].R )
            return true;
    }
    else {
        wsrc = wkid + 1;

        if( wsrc >= nwks )
            return true;

        if( !gdrty[wsrc].L )
            return true;
    }

    for( int iz = zlo; iz <= zhi; ++iz ) {

        buf		= (void*)&K[iz][0];
        bytes	= sizeof(int) * vR[iz].nr;

        MPIRecv( buf, bytes, wsrc, iz - zlo );
    }

    return true;
}

/* --------------------------------------------------------------- */
/* KUpdt --------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Split::KUpdt()
{
    if( nwks <= 1 )
        return true;

// To avoid deadlocks, we arrange that at any given time,
// any node is just sending or just receiving. Easily done
// by odd/even role assignment.

// Odd send

    if( wkid & 1 ) {

        KSend( zLlo, zLhi, 'L' );
        KSend( zRlo, zRhi, 'R' );
    }
    else {

        KRecv( zihi + 1, zohi, 'R' );
        KRecv( zolo, zilo - 1, 'L' );
    }

// Even send

    if( !(wkid & 1) ) {

        KSend( zLlo, zLhi, 'L' );
        KSend( zRlo, zRhi, 'R' );
    }
    else {

        KRecv( zihi + 1, zohi, 'R' );
        KRecv( zolo, zilo - 1, 'L' );
    }

    return true;
}

/* --------------------------------------------------------------- */
/* _Propagate1 --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _Propagate1( void* ithr )
{
    int	nd = vzd.size();

    const Rgns&		Rs = vR[zs];
    vector<int>&	ks = ME->K[zs];

    for( int id = (long)ithr; id < nd; id += nthr ) {

        int				zd = vzd[id];
        vector<int>&	kd = ME->K[zd];
        int				nd = vR[zd].nr;

        for( int ir = 0; ir < Rs.nr; ++ir ) {

            int	ksrc = ks[ir];

            // Valid tile?
            if( !ksrc )
                continue;

            const vector<int>&	P  = Rs.pts[ir];
            int					np = P.size();

            for( int ip = 0; ip < np; ++ip ) {

                const CorrPnt&	C = vC[P[ip]];

                if( !C.used )
                    continue;

                int	jr;

                if( C.z1 == zd )
                    jr = C.i1;
                else if( C.z2 == zd )
                    jr = C.i2;
                else
                    continue;

                int	kdst = kd[jr];

                if( !kdst )
                    continue;

                if( kdst == ksrc )
                    continue;

                // Note any difference by setting the
                // drtyD flag now. But at this time,
                // we only update dst colors if the
                // src color is lower. If the dst is
                // lower, it will propagate to us on
                // a later loop iteration.

                drtyD[zd] = 1;

                if( ksrc < kdst ) {
                    for( int ik = 0; ik < nd; ++ik ) {
                        if( kd[ik] == kdst )
                            kd[ik] = ksrc;
                    }
                }
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* Propagate1 ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Propagate colors from source layer zs to its connectees.
//
// To do this with threads that don't collide, all the threads
// will read from same source layer zs, but each will write to
// a unique dst layer zd.
//
void Split::Propagate1()
{
    vzd.clear();

    int	zlo = max( zs - maxspan, zilo ),
        zhi = min( zs + maxspan, zihi ),
        nd  = zhi - zlo + 1;

    for( int i = 0; i < nd; ++i ) {
        if( i + zlo != zs )
            vzd.push_back( i + zlo );
    }

    ME		= (Split*)this;
    nd		= vzd.size();
    nthr	= maxthreads;

    if( nthr > nd )
        nthr = nd;

    if( !EZThreads( _Propagate1, nthr, 1, "_Propagate1" ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* PropagateLocally ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Propagate each dirty layer in a loop until no more dirty
// flags in our local block's interior. However, gdrty flags
// cummulatively record any changes to communicate globally.
//
void Split::PropagateLocally()
{
// clear cummulative dst flags
    gdrty[wkid].L = 0;
    gdrty[wkid].R = 0;

    for(;;) {

        memset( &drtyD[0], 0, zohi - zolo + 1 );

        /* ---------- */
        /* Do a round */
        /* ---------- */

        for( zs = zolo; zs <= zohi; ++zs ) {

            if( drtyS[zs] )
                Propagate1();
        }

        /* ------------------ */
        /* Neighbors touched? */
        /* ------------------ */

        if( wkid > 0 ) {
            for( int iz = zLlo; iz <= zLhi; ++iz )
                gdrty[wkid].L |= drtyD[iz];
        }

        if( wkid + 1 < nwks ) {
            for( int iz = zRlo; iz <= zRhi; ++iz )
                gdrty[wkid].R |= drtyD[iz];
        }

        /* ----- */
        /* Done? */
        /* ----- */

        memcpy( &drtyS[0], &drtyD[0], zohi - zolo + 1 );

        int	drty = 0;
        for( int iz = zilo; iz <= zihi; ++iz )
            drty |= drtyD[iz];

        if( !drty )
            break;
    }
}

/* --------------------------------------------------------------- */
/* ReportCount --------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::ReportCount( const map<int,int>& m )
{
    map<int,int>::const_iterator	it, ie = m.end();

    for( it = m.begin(); it != ie; ++it )
        printf( "Color %9d Count %9d\n", it->first, it->second );
}

/* --------------------------------------------------------------- */
/* CountColors --------------------------------------------------- */
/* --------------------------------------------------------------- */

class mpiclr {
// share color map with workers
public:
    int	first, second;
public:
    mpiclr() {};
    mpiclr( map<int,int>::iterator& it )
    : first(it->first), second(it->second) {};
};

void Split::CountColors( map<int,int>& m )
{
// Cache iterator and end()

    map<int,int>::iterator	it, ie = m.end();
    int						clast = -1;

// Calc my colors

    for( int iz = zilo; iz <= zihi; ++iz ) {

        vector<int>&	k = K[iz];
        int				nr = vR[iz].nr;

        for( int ir = 0; ir < nr; ++ir ) {

            int	clr = k[ir];

            if( !clr )
                continue;

            if( clr != clast )
                it = m.find( clast = clr );

            if( it != ie )
                ++it->second;
            else {	// new entry: reset cache
                m[clr]	= 1;
                ie		= m.end();
                clast	= -1;
            }
        }
    }

// Report my subcounts

    if( nwks > 1 )
        printf( "This worker:\n" );
    else
        printf( "All workers:\n" );

    ReportCount( m );

// Send my colors to master; then get global colors

    if( wkid > 0 ) {

        // send count
        int	nm = m.size();
        MPISend( &nm, sizeof(int), 0, wkid );

        // send elems
        for( it = m.begin(); it != ie; ++it ) {
            mpiclr	c( it );
            MPISend( &c, sizeof(mpiclr), 0, wkid );
        }

        // get count
        m.clear();
        MPIRecv( &nm, sizeof(int), 0, wkid );

        // get elems
        for( int i = 0; i < nm; ++i ) {
            mpiclr	c;
            MPIRecv( &c, sizeof(mpiclr), 0, wkid );
            m[c.first] = c.second;
        }
    }
    else if( nwks > 1 ) {

        int	nm;

        // get each worker
        for( int iw = 1; iw < nwks; ++iw ) {

            MPIRecv( &nm, sizeof(int), iw, iw );

            for( int i = 0; i < nm; ++i ) {

                mpiclr	c;
                MPIRecv( &c, sizeof(mpiclr), iw, iw );

                it = m.find( c.first );

                if( it != m.end() )
                    it->second += c.second;
                else
                    m[c.first] = c.second;
            }
        }

        // send each worker
        nm = m.size();
        ie = m.end();

        for( int iw = 1; iw < nwks; ++iw ) {

            MPISend( &nm, sizeof(int), iw, iw );

            for( it = m.begin(); it != ie; ++it ) {
                mpiclr	c( it );
                MPISend( &c, sizeof(mpiclr), iw, iw );
            }
        }

        printf( "\nAll workers:\n" );
        ReportCount( m );
    }
}

/* --------------------------------------------------------------- */
/* _Save --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SaveXBin( const vector<double> &x, int z )
{
    char	buf[64];
    FILE	*f;

    sprintf( buf, "%s/X_%c_%d.bin",
        gpath, (ME->X.NE == 6 ? 'A' : 'H'), z );

    f = FileOpenOrDie( buf, "wb" );
    fwrite( &x[0], sizeof(double), x.size(), f );
    fclose( f );
}


static void SaveFBin( const vector<uint8> &f, int z )
{
    char	buf[64];
    FILE	*q;

    sprintf( buf, "%s/F_%d.bin", gpath, z );
    q = FileOpenOrDie( buf, "wb" );
    fwrite( &f[0], sizeof(uint8), f.size(), q );
    fclose( q );
}


static void* _Save( void* ithr )
{
    for( int iz = zilo + (long)ithr; iz <= zihi; iz += nthr ) {

        const Rgns&		R = vR[iz];
        vector<int>&	k = ME->K[iz];
        vector<uint8>	f = R.flag;

        // Flag adjustment:
        // Basically, we preserve existing flags,
        // but if it's not the right color, mark
        // it as 'transform missing'.
        //
        // Saveclr is set by caller:
        // >0 => match specific color
        // -1 => consolidate remaining colors
        //  0 => not implemented.

        if( saveclr >= 0 ) {

            for( int j = 0; j < R.nr; ++j ) {

                if( k[j] == saveclr )
                    k[j] = -1;	// mark it removed
                else			// kill all others
                    f[j] = fmRead;
            }
        }
        else {	// save all positive colors together

            for( int j = 0; j < R.nr; ++j ) {

                if( k[j] <= 0 )	// not these
                    f[j] = fmRead;
            }
        }

        SaveXBin( ME->X.X[iz], R.z );
        SaveFBin( f, R.z );
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* Save ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Save tforms and <<modified>> flags...
//
void Split::Save()
{
    char	buf[64];
    int		nz = zihi - zilo + 1;

    ME		= (Split*)this;
    gpath	= buf;
    nthr	= maxthreads;

    if( nthr > nz )
        nthr = nz;

    if( saveclr > 0 ) {
        sprintf( buf, "Splits/X_%c_BIN_CLR%d",
        (X.NE == 6 ? 'A' : 'H'), saveclr );
    }
    else {
        sprintf( buf, "Splits/X_%c_BIN_REM",
        (X.NE == 6 ? 'A' : 'H') );
    }

    DskCreateDir( buf, stdout );

    if( !EZThreads( _Save, nthr, 1, "_SplitSave" ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* BreakOut ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void Split::BreakOut( const map<int,int>& m )
{
    int	nm = m.size();

    if( nm <= 1 ) {
        printf(
        "No splits created;"
        " (single color or all below splitmin).\n" );
        return;
    }

// Splits folder

    DskCreateDir( "Splits", stdout );

// First separately write each color over threshold.
// The Save operation will replace that color by (-1)
// removing it from the K vectors.

    map<int,int>::const_iterator	it, ie = m.end();

    for( it = m.begin(); it != ie; ++it ) {

        if( it->second >= splitmin ) {

            saveclr = it->first;
            Save();
            --nm;
        }
    }

// Now consolidate remaining 'REM' fragments.

    if( nm > 0 ) {
        saveclr = -1;
        Save();
    }
}

/* --------------------------------------------------------------- */
/* Run ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Split::Run()
{
    printf( "\n---- Split ----\n" );

    clock_t	t0 = StartTiming();

// Init color workspace

    Resize();

// First colorize montages alone

    ColorMontages();

// Initially everything dirty

    gdrty.resize( nwks );
    gdrty[0].L			= 0;
    gdrty[nwks - 1].R	= 0;

    drtyS.resize( zohi - zolo + 1, 1 );
    drtyD.resize( zohi - zolo + 1, 0 );

    maxspan = LayerCat_MaxSpan( vL );

// Propagate colors across Z until no changes

    for(;;) {

        /* ---------- */
        /* Do a round */
        /* ---------- */

        // get src data
        KUpdt();

        // set src flags
        if( wkid > 0 && gdrty[wkid - 1].R ) {

            for( int iz = zolo; iz < zilo; ++iz )
                drtyS[iz] = 1;
        }

        if( wkid + 1 < nwks && gdrty[wkid + 1].L ) {

            for( int iz = zihi + 1; iz <= zohi; ++iz )
                drtyS[iz] = 1;
        }

        PropagateLocally();

        /* -------------- */
        /* Done globally? */
        /* -------------- */

        GUpdt();

        int	drty = 0;
        for( int iw = 0; iw < nwks; ++iw )
            drty |= gdrty[iw].L | gdrty[iw].R;

        if( !drty )
            break;
    }

    vzd.clear();
    drtyD.clear();
    drtyS.clear();
    gdrty.clear();

// Split according to final coloring

    map<int,int>	m;

    CountColors( m );
    BreakOut( m );

    StopTiming( stdout, "Split", t0 );
}


