

#include	"Disk.h"
#include	"File.h"
#include	"../1_Cross_LowRes/ScapeMeta.h"






/* --------------------------------------------------------------- */
/* ReadLog1 ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool ReadLog1( CLog &L, int z )
{
    char		buf[256];
    sprintf( buf, "scplogs/scp_%d.log", z );

    if( !DskExists( buf ) )
        return false;

    FILE		*f = FileOpenOrDie( buf, "r" );
    CLineScan	LS;
    bool		gotM = false;

    while( LS.Get( f ) > 0 ) {

        if( LS.line[0] != '*' )
            continue;

        char	key = LS.line[1];

        LS.Get( f );

        if( key == 'M' ) {

            CScapeMeta	&E = L.M;

            sscanf( LS.line,
            "%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
            &E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
            &E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );

            gotM = true;
        }
        else if( key == 'A' ) {

            CScapeMeta	&E = L.A;

            sscanf( LS.line,
            "%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
            &E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
            &E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );
        }
        else if( key == 'B' ) {

            CScapeMeta	&E = L.B;

            sscanf( LS.line,
            "%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
            &E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
            &E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );
        }
        else if( key == 'T' ) {

            double	*t = L.T.t;

            sscanf( LS.line + 1,
            "%lf,%lf,%lf,%lf,%lf,%lf",
            &t[0], &t[1], &t[2],
            &t[3], &t[4], &t[5] );

            break;
        }
    }

    fclose( f );
    return gotM;
}

/* --------------------------------------------------------------- */
/* ReadLogs ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int ReadLogs( vector<CLog> &vL, int zmin, int zmax )
{
    for( int z = zmin; z <= zmax; ) {

        CLog	L;

        L.A.z = -1;

        if( ReadLog1( L, z ) ) {

            vL.push_back( L );

            if( L.A.z == -1 )
                break;

            z = L.A.z;
        }
        else
            ++z;
    }

    return vL.size();
}


