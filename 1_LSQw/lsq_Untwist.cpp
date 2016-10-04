

#include	"lsq_Globals.h"
#include	"lsq_MPI.h"
#include	"lsq_Untwist.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"CRigid.h"
#include	"Timer.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class AIO {
public:
    int		z;
    TAffine	A;
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray		*gX;
static vector<AIO>	vA;
static int			nthr;






/* --------------------------------------------------------------- */
/* _RgdSums ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void* _RgdSums( void* ithr )
{
    int	nb = vA.size();

// For each B-layer...

    for( int ib = (long)ithr; ib < nb; ib += nthr ) {

        // Form the sums

        CRigid					rgd;
        int						ia	= ib + 1;
        const Rgns&				Rb	= vR[ib];
        const Rgns&				Ra	= vR[ia];
        const vector<double>&	xb	= gX->X[ib];
        const vector<double>&	xa	= gX->X[ia];
        AIO&					A	= vA[ib];

        // For each A-layer rgn...

        for( int ir = 0; ir < Ra.nr; ++ir ) {

            if( !FLAG_ISUSED( Ra.flag[ir] ) )
                continue;

            const vector<int>&	P  = Ra.pts[ir];
            const TAffine*		Ta = &X_AS_AFF( xa, ir );
            const TAffine*		Tb;
            int					lastb	= -1,
                                np		= P.size();

            // For each of its points...

            for( int ip = 0; ip < np; ++ip ) {

                const CorrPnt&	C = vC[P[ip]];

                // Want only zA onto zB

                if( C.z1 != ia || C.z2 != ib )
                    continue;

                if( C.i2 != lastb ) {

                    if( !FLAG_ISUSED( Rb.flag[C.i2] ) )
                        continue;

                    Tb = &X_AS_AFF( xb, C.i2 );
                    lastb = C.i2;
                }

                Point	pa = C.p1,
                        pb = C.p2;

                Ta->Transform( pa );
                Tb->Transform( pb );
                rgd.Add( pa, pb );
            }
        }

        // Set transform

        A.z = Ra.z;
        rgd.Solve( A.A );
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* CalcMyPairwiseTForms ------------------------------------------ */
/* --------------------------------------------------------------- */

// Pairwise layer sums and tforms can be done in parallel.
//
static void CalcMyPairwiseTForms( XArray &X )
{
    gX = &X;

    int	nb = zohi - zolo;	// this many b-layers

    vA.resize( nb );

    nthr = maxthreads;

    if( nthr > nb )
        nthr = nb;

    if( !EZThreads( _RgdSums, nthr, 1, "_RgdSums" ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteMyTForms ------------------------------------------------- */
/* --------------------------------------------------------------- */

// The adjustment tform needed for a given layer is the cummulative
// product of all adjustments to layers below. At this point, all
// worker nodes except the last must write their pairwise results
// to files "Untwist/id.bin".
//
static void WriteMyTForms()
{
// If I'm not the last worker...

    if( wkid >= nwks - 1 )
        return;

// Each entry affects layer zA or higher, so write blocks:
// 'zA A'

    char	buf[32];
    DskCreateDir( "Untwist", stdout );
    sprintf( buf, "Untwist/%d.bin", wkid );

    FILE	*f = FileOpenOrDie( buf, "wb" );
    fwrite( &vA[0], sizeof(AIO), vA.size(), f );
    fclose( f );
}

/* --------------------------------------------------------------- */
/* AccumulateBefores --------------------------------------------- */
/* --------------------------------------------------------------- */

// Form product of all affines with z <= my stating z.
//
// We must be careful here because workers have overlapping
// zo ranges, while our product must include each z in order
// and once only. This is fixed just by requiring monotonic
// increase in z.
//
static void AccumulateBefores( TAffine &A0 )
{
    int	z0		= vR[zolo].z,
        zlast	= -1;

    for( int iw = 0; iw < wkid; ++iw ) {

        char	buf[32];
        sprintf( buf, "Untwist/%d.bin", iw );

        int			n = (int)DskBytes( buf ) / sizeof(AIO);
        vector<AIO>	aio( n );

        FILE	*f = FileOpenOrDie( buf, "rb" );
        fread( &aio[0], sizeof(AIO), n, f );
        fclose( f );

        for( int i = 0; i < n; ++i ) {

            const AIO&	A = aio[i];

            if( A.z > z0 )
                return;

            // monotonic z rule

            if( A.z <= zlast )
                continue;

            A0		= A.A * A0;
            zlast	= A.z;
        }
    }
}

/* --------------------------------------------------------------- */
/* Apply --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Apply A0 and pairwise affines.
//
static void Apply( TAffine &A0 )
{
    int	nz = zohi - zolo + 1;

    for( int iz = (wkid ? 0 : 1); iz < nz; ++iz ) {

        if( iz )
            A0 = vA[iz - 1].A * A0;

        const Rgns&		R = vR[iz];
        vector<double>&	x = gX->X[iz];

        for( int ir = 0; ir < R.nr; ++ir ) {

            if( !FLAG_ISUSED( R.flag[ir] ) )
                continue;

            TAffine& T = X_AS_AFF( x, ir );

            T = A0 * T;
        }
    }
}

/* --------------------------------------------------------------- */
/* UntwistAffines ------------------------------------------------ */
/* --------------------------------------------------------------- */

// The standard external scaffold for use with the -prior option
// is 'cross_wkspc/X_A_BIN_scaf'; the result of aligning low-res
// strips from very good montages. However, we can make better
// angle calculations later, after the down correspondence points
// are formed. Adjustments are calculated here as rigid transforms
// T{theta, kx, ky) from layer a to b.
//
void UntwistAffines( XArray &X )
{
    if( X.NE != 6 || zolo >= zohi )
        return;

    clock_t	t0 = StartTiming();

    CalcMyPairwiseTForms( X );
    WriteMyTForms();

    MPIWaitForOthers();

    TAffine	A0;
    AccumulateBefores( A0 );

    Apply( A0 );
    vA.clear();

// Cleanup

    MPIWaitForOthers();

    if( nwks > 1 && !wkid )
        system( "rm -rf Untwist" );

    StopTiming( stdout, "Untwist", t0 );
}


