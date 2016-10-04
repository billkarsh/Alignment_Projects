

#include	"lsq_Dropout.h"
#include	"lsq_Globals.h"
#include	"lsq_MPI.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static vector<Dropout>	vD;
static int				nthr;






/* --------------------------------------------------------------- */
/* _Scan --------------------------------------------------------- */
/* --------------------------------------------------------------- */

void* _Scan( void* ithr )
{
    int	nd = zihi - zilo + 1;

// For each layer...

    for( int id = (long)ithr; id < nd; id += nthr ) {

        int			iz	= id + zilo;
        const Rgns&	R	= vR[iz];
        Dropout&	D	= vD[id];
        FILE		*q	= NULL;

        D.rmax += R.nr;

        // For each rgn...

        for( int ir = 0; ir < R.nr; ++ir ) {

            if( R.flag[ir] ) {

                if( !q ) {
                    DskCreateDir( "Dropouts", stdout );
                    char	buf[64];
                    sprintf( buf, "Dropouts/drop_%d.txt", R.z );
                    q = FileOpenOrDie( buf, "w" );
                }

                int	z, i, r;
                RealZIDR( z, i, r, iz, ir );

                if( FLAG_ISREAD( R.flag[ir] ) ) {
                    fprintf( q, "R %d.%d-%d\n", z, i, r );
                    ++D.read;
                }
                else if( FLAG_ISPNTS( R.flag[ir] ) ) {
                    fprintf( q, "P %d.%d-%d\n", z, i, r );
                    ++D.pnts;
                }
                else if( FLAG_ISKILL( R.flag[ir] ) ) {
                    fprintf( q, "K %d.%d-%d\n", z, i, r );
                    ++D.kill;
                }
                else if( FLAG_ISCUTD( R.flag[ir] ) ) {
                    fprintf( q, "C %d.%d-%d\n", z, i, r );
                    ++D.cutd;
                }
            }
        }

        if( q )
            fclose( q );
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* ScanEachLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScanEachLayer()
{
    int	nd = zihi - zilo + 1;

    vD.resize( nd );

    nthr = maxthreads;

    if( nthr > nd )
        nthr = nd;

    if( !EZThreads( _Scan, nthr, 1, "_ScanEachLayer" ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* GatherCounts -------------------------------------------------- */
/* --------------------------------------------------------------- */

void Dropout::GatherCounts()
{
    int	nd = zihi - zilo + 1;

    for( int id = 0; id < nd; ++id )
        Add( vD[id] );

    if( nwks > 1 ) {
        printf(
        "Worker %03d: MAX %9ld READ %9ld"
        " PNTS %9ld KILL %9ld CUTD %9ld\n",
        wkid, rmax, read, pnts, kill, cutd );
    }
    else {
        printf(
        "All workers: MAX %9ld READ %9ld"
        " PNTS %9ld KILL %9ld CUTD %9ld\n",
        rmax, read, pnts, kill, cutd );
    }

    if( wkid > 0 )
        MPISend( this, sizeof(Dropout), 0, wkid );
    else if( nwks > 1 ) {

        for( int iw = 1; iw < nwks; ++iw ) {

            Dropout	D;
            MPIRecv( &D, sizeof(Dropout), iw, iw );
            Add( D );

            printf(
            "Worker %03d: MAX %9ld READ %9ld"
            " PNTS %9ld KILL %9ld CUTD %9ld\n",
            iw, D.rmax, D.read, D.pnts, D.kill, D.cutd );
        }

        printf(
        "--------------------------------------------"
        "-------------------------------------------\n"
        "     Total: MAX %9ld READ %9ld"
        " PNTS %9ld KILL %9ld CUTD %9ld\n\n",
        rmax, read, pnts, kill, cutd );
    }
}

/* --------------------------------------------------------------- */
/* Scan ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void Dropout::Scan()
{
    printf( "\n---- Dropouts ----\n" );

    clock_t	t0 = StartTiming();

    ScanEachLayer();
    GatherCounts();
    vD.clear();

    StopTiming( stdout, "Drops ", t0 );
}


