

#include	"lsq_Globals.h"
#include	"lsq_MPI.h"
#include	"lsq_XArray.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"Timer.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray*		ME;
static const char	*gpath;
static vector<int>	giz;
static int			nthr;






/* --------------------------------------------------------------- */
/* _AFromIDB ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _AFromIDB( void* ithr )
{
    for( int iz = zolo + (long)ithr; iz <= zohi; iz += nthr ) {

        Rgns&			R = vR[iz];
        vector<double>&	x = ME->X[iz];
        vector<Til2Img>	t2i;

        // Get rgn #1 tforms

        if( !IDBT2IGet_JustIDandT( t2i, idb, R.z ) )
            exit( 42 );

        x.resize( R.nr * 6 );

        int						nt = t2i.size();
        map<int,int>::iterator	en = R.m.end();

        // For each transform in IDB...

        for( int it = 0; it < nt; ++it ) {

            // Get its block start and limit {j0,jlim}

            const Til2Img&			T = t2i[it];
            map<int,int>::iterator	mi = R.m.find( T.id );
            int						j0, jlim;

            if( mi == en )
                continue;

            j0		= mi->second;
            jlim	= (++mi != en ? mi->second : R.nr);

            // Propagate rgn #1 tform to all block members

            for( int j = j0; j < jlim; ++j ) {

                if( R.pts[j].size() >= 3 ) {

                    T.T.CopyOut( &x[j * 6] );
                    FLAG_SETUSED( R.flag[j] );
                }
                else
                    FLAG_SETPNTS( R.flag[j] );
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* _AFromTxt ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// There are two differences in reading 'TXT' vs 'MET' formats:
//
// (1) The folder path name; but that's absorbed in gpath.
//
// (2) We scan each line for {id,r,T} but there is more beyond
//		in the 'MET' case. That's handled in CIDRA::FromFile()
//		by format '%*[^\r\n][\r\n]' which reads and tosses all
//		up to and including the terminator.
//

class CIDRA {
public:
    TAffine	A;
    int		id, r;
public:
    inline bool FromFile( FILE *f )
    {
        return 8 == fscanf( f,
            " %d %d %lf %lf %lf %lf %lf %lf%*[^\r\n][\r\n]",
            &id, &r,
            &A.t[0], &A.t[1], &A.t[2],
            &A.t[3], &A.t[4], &A.t[5] );
    };
};


static bool Read_vA( vector<CIDRA> &vA, int z )
{
    char	buf[2048];
    FILE	*f;
    CIDRA	A;
    bool	nf = true;	// default = no folds

    sprintf( buf, "%s/X_A_%d.txt", gpath, z );
    f = FileOpenOrDie( buf, "r" );

    while( A.FromFile( f ) ) {

        if( A.r > 1 )
            nf = false;

        vA.push_back( A );
    }

    fclose( f );

    return nf;
}


static void* _AFromTxt( void* ithr )
{
    for( int iz = zolo + (long)ithr; iz <= zohi; iz += nthr ) {

        Rgns&			R = vR[iz];
        vector<double>&	x = ME->X[iz];
        vector<CIDRA>	vA;
        int				in = (iz >= zilo && iz <= zihi ),
                        nf = Read_vA( vA, R.z );

        x.resize( R.nr * 6 );

        int						na = vA.size();
        map<int,int>::iterator	en = R.m.end();

        if( !nf ) {	// Propagate rgn #1 to all block members

            // For each transform in vA...

            for( int ia = 0; ia < na; ++ia ) {

                // Get its block start and limit {j0,jlim}

                const CIDRA&			A = vA[ia];
                map<int,int>::iterator	mi = R.m.find( A.id );
                int						j0, jlim;

                if( mi == en )
                    continue;

                j0		= mi->second;
                jlim	= (++mi != en ? mi->second : R.nr);

                // Propagate rgn #1 tform to all block members

                for( int j = j0; j < jlim; ++j ) {

                    if( !in || R.pts[j].size() >= 3 ) {
                        A.A.CopyOut( &x[j * 6] );
                        FLAG_SETUSED( R.flag[j] );
                    }
                    else
                        FLAG_SETPNTS( R.flag[j] );
                }
            }
        }
        else {	// Move each A to its specified position

            // For each transform in vA...

            for( int ia = 0; ia < na; ++ia ) {

                // Get its location j

                const CIDRA&			A = vA[ia];
                map<int,int>::iterator	mi = R.m.find( A.id );
                int						j;

                if( mi == en )
                    continue;

                j = mi->second + A.r - 1;

                if( !in || R.pts[j].size() >= 3 ) {
                    A.A.CopyOut( &x[j * 6] );
                    FLAG_SETUSED( R.flag[j] );
                }
                else
                    FLAG_SETPNTS( R.flag[j] );
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* _HFromTxt ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// There are two differences in reading 'TXT' vs 'MET' formats:
//
// (1) The folder path name; but that's absorbed in gpath.
//
// (2) We scan each line for {id,r,T} but there is more beyond
//		in the 'MET' case. That's handled in CIDRH::FromFile()
//		by format '%*[^\r\n][\r\n]' which reads and tosses all
//		up to and including the terminator.
//

class CIDRH {
public:
    THmgphy	H;
    int		id, r;
public:
    inline bool FromFile( FILE *f )
    {
        return 10 == fscanf( f,
            " %d %d %lf %lf %lf %lf %lf %lf %lf %lf%*[^\r\n][\r\n]",
            &id, &r,
            &H.t[0], &H.t[1], &H.t[2],
            &H.t[3], &H.t[4], &H.t[5],
            &H.t[6], &H.t[7] );
    };
};


static bool Read_vH( vector<CIDRH> &vH, int z )
{
    char	buf[2048];
    FILE	*f;
    CIDRH	H;
    bool	nf = true;	// default = no folds

    sprintf( buf, "%s/X_H_%d.txt", gpath, z );
    f = FileOpenOrDie( buf, "r" );

    while( H.FromFile( f ) ) {

        if( H.r > 1 )
            nf = false;

        vH.push_back( H );
    }

    fclose( f );

    return nf;
}


static void* _HFromTxt( void* ithr )
{
    for( int iz = zolo + (long)ithr; iz <= zohi; iz += nthr ) {

        Rgns&			R = vR[iz];
        vector<double>&	x = ME->X[iz];
        vector<CIDRH>	vH;
        int				in = (iz >= zilo && iz <= zihi ),
                        nf = Read_vH( vH, R.z );

        x.resize( R.nr * 8 );

        int						nh = vH.size();
        map<int,int>::iterator	en = R.m.end();

        if( !nf ) {	// Propagate rgn #1 to all block members

            // For each transform in vH...

            for( int ih = 0; ih < nh; ++ih ) {

                // Get its block start and limit {j0,jlim}

                const CIDRH&			H = vH[ih];
                map<int,int>::iterator	mi = R.m.find( H.id );
                int						j0, jlim;

                if( mi == en )
                    continue;

                j0		= mi->second;
                jlim	= (++mi != en ? mi->second : R.nr);

                // Propagate rgn #1 tform to all block members

                for( int j = j0; j < jlim; ++j ) {

                    if( !in || R.pts[j].size() >= 4 ) {
                        H.H.CopyOut( &x[j * 8] );
                        FLAG_SETUSED( R.flag[j] );
                    }
                    else
                        FLAG_SETPNTS( R.flag[j] );
                }
            }
        }
        else {	// Move each H to its specified position

            // For each transform in vH...

            for( int ih = 0; ih < nh; ++ih ) {

                // Get its location j

                const CIDRH&			H = vH[ih];
                map<int,int>::iterator	mi = R.m.find( H.id );
                int						j;

                if( mi == en )
                    continue;

                j = mi->second + H.r - 1;

                if( !in || R.pts[j].size() >= 4 ) {
                    H.H.CopyOut( &x[j * 8] );
                    FLAG_SETUSED( R.flag[j] );
                }
                else
                    FLAG_SETPNTS( R.flag[j] );
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* _XFromBin ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadXBin( vector<double> &x, int z )
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/X_%c_%d.bin",
        gpath, (ME->NE == 6 ? 'A' : 'H'), z );

    f = FileOpenOrDie( buf, "rb" );
    fread( &x[0], sizeof(double), x.size(), f );
    fclose( f );
}


static void ReadFBin( vector<uint8> &f, int z )
{
    char	buf[2048];
    FILE	*q;

    sprintf( buf, "%s/F_%d.bin", gpath, z );
    q = FileOpenOrDie( buf, "rb" );
    fread( &f[0], sizeof(uint8), f.size(), q );
    fclose( q );
}


static void* _XFromBin( void* ithr )
{
    for( int iz = zolo + (long)ithr; iz <= zohi; iz += nthr ) {

        Rgns&			R = vR[iz];
        vector<double>&	x = ME->X[iz];
        int				minpts = ME->NE / 2;

        x.resize( R.nr * ME->NE );
        ReadXBin( x, R.z );
        ReadFBin( R.flag, R.z );

        // If read tforms are in wings we adopt their
        // flags as is and our only concern is whether
        // they are used or not via FLAG_ISUSED().
        //
        // If in our inner range, we test point count.
        //
        if( iz >= zilo && iz <= zihi ) {

            for( int j = 0; j < R.nr; ++j ) {

                if( R.pts[j].size() < minpts )
                    FLAG_ADDPNTS( R.flag[j] );
            }
        }
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* _Save --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SaveXBin( const vector<double> &x, int z )
{
    char	buf[64];
    FILE	*f;

    sprintf( buf, "%s/X_%c_%d.bin",
        gpath, (ME->NE == 6 ? 'A' : 'H'), z );

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

        const Rgns&				R = vR[iz];
        const vector<double>&	x = ME->X[iz];

        SaveXBin( x, R.z );
        SaveFBin( R.flag, R.z );
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* _UpdtFS ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _UpdtFS( void* ithr )
{
    int	nz = giz.size();

    for( int iz = (long)ithr; iz < nz; iz += nthr ) {

        Rgns&			R = vR[giz[iz]];
        vector<double>&	x = ME->X[giz[iz]];

        ReadXBin( x, R.z );
        ReadFBin( R.flag, R.z );
    }

    return NULL;
}

/* --------------------------------------------------------------- */
/* Resize -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Resize( int ne )
{
    NE = ne;

    int	nz = zohi - zolo + 1;

    X.resize( nz );

    for( int iz = 0; iz < nz; ++iz )
        X[iz].resize( vR[iz].nr * ne );
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Load( const char *path )
{
    clock_t	t0 = StartTiming();

    int	nz = zohi - zolo + 1;

    X.resize( nz );

    ME		= this;
    gpath	= path;
    nthr	= maxthreads;

    if( nthr > nz )
        nthr = nz;

    EZThreadproc	proc;
    const char		*sproc;

    if( !path || !path[0] ) {
        NE		= 6;
        proc	= _AFromIDB;
        sproc	= "AFromIDB";
    }
    else {

        const char *name = FileNamePtr( path );

        if( strstr( name, "X_A" ) ) {

            NE = 6;

            if( strstr( name, "X_A_TXT" ) ||
                strstr( name, "X_A_MET" ) ) {

                proc	= _AFromTxt;
                sproc	= "AFromTxt";
            }
            else if( strstr( name, "X_A_BIN" ) ) {
                proc	= _XFromBin;
                sproc	= "XFromBin";
            }
            else
                goto error;
        }
        else if( strstr( name, "X_H" ) ) {

            NE = 8;

            if( strstr( name, "X_H_TXT" ) ||
                strstr( name, "X_H_MET" ) ) {

                proc	= _HFromTxt;
                sproc	= "HFromTxt";
            }
            else if( strstr( name, "X_H_BIN" ) ) {
                proc	= _XFromBin;
                sproc	= "XFromBin";
            }
            else
                goto error;
        }
        else {
error:
            printf( "Unknown prior type [%s].\n", name );
            exit( 42 );
        }
    }

    if( !EZThreads( proc, nthr, 1, sproc ) )
        exit( 42 );

    StopTiming( stdout, sproc, t0 );
}

/* --------------------------------------------------------------- */
/* Save ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Save tforms and <<modified>> flags...
//
void XArray::Save() const
{
    clock_t	t0 = StartTiming();

    int	nz = nz = zihi - zilo + 1;

    char	buf[32];
    sprintf( buf, "X_%c_BIN", (NE == 6 ? 'A' : 'H') );
    DskCreateDir( buf, stdout );

    ME		= (XArray*)this;
    gpath	= buf;
    nthr	= maxthreads;

    if( nthr > nz )
        nthr = nz;

    if( !EZThreads( _Save, nthr, 1, "_XArraySave" ) )
        exit( 42 );

    StopTiming( stdout, "Save", t0 );
}

/* --------------------------------------------------------------- */
/* UpdtFS -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Older flavor of updating across nodes in which everyone
// writes their zi range and then reads their wings using
// this function. Now replaced by MPI send/receive.
//
void XArray::UpdtFS()
{
    clock_t	t0 = StartTiming();

    giz.clear();

    for( int i = zolo; i < zilo; ++i )
        giz.push_back( i );

    for( int i = zihi + 1; i <= zohi; ++i )
        giz.push_back( i );

    int	nz = giz.size();

    if( !nz )
        return;

    char	buf[32];
    sprintf( buf, "X_%c_BIN", (NE == 6 ? 'A' : 'H') );

    ME		= this;
    gpath	= buf;
    nthr	= maxthreads;

    if( nthr > nz )
        nthr = nz;

    if( !EZThreads( _UpdtFS, nthr, 1, "_UpdtFS" ) )
        exit( 42 );

    giz.clear();

    StopTiming( stdout, "Updt", t0 );
}

/* --------------------------------------------------------------- */
/* Send ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool XArray::Send( int zlo, int zhi, int XorF, int toLorR )
{
    void	*buf;
    int		bytes,
            wdst;

    if( toLorR == 'L' ) {

        wdst = wkid - 1;

        if( wdst < 0 )
            return true;
    }
    else {
        wdst = wkid + 1;

        if( wdst >= nwks )
            return true;
    }

    for( int iz = zlo; iz <= zhi; ++iz ) {

        if( XorF == 'X' ) {
            buf		= (void*)&X[iz][0];
            bytes	= sizeof(double) * X[iz].size();
        }
        else {
            buf		= (void*)&vR[iz].flag[0];
            bytes	= vR[iz].nr;
        }

        MPISend( buf, bytes, wdst, iz - zlo );
    }

    return true;
}

/* --------------------------------------------------------------- */
/* Recv ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool XArray::Recv( int zlo, int zhi, int XorF, int fmLorR )
{
    void	*buf;
    int		bytes,
            wsrc;

    if( fmLorR == 'L' ) {

        wsrc = wkid - 1;

        if( wsrc < 0 )
            return true;
    }
    else {
        wsrc = wkid + 1;

        if( wsrc >= nwks )
            return true;
    }

    for( int iz = zlo; iz <= zhi; ++iz ) {

        if( XorF == 'X' ) {
            buf		= (void*)&X[iz][0];
            bytes	= sizeof(double) * X[iz].size();
        }
        else {
            buf		= (void*)&vR[iz].flag[0];
            bytes	= vR[iz].nr;
        }

        MPIRecv( buf, bytes, wsrc, iz - zlo );
    }

    return true;
}

/* --------------------------------------------------------------- */
/* Updt ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool XArray::Updt()
{
    if( nwks <= 1 )
        return true;

// To avoid deadlocks, we arrange that at any given time,
// any node is just sending or just receiving. Easily done
// by odd/even role assignment.

// Odd send

    if( wkid & 1 ) {

        Send( zLlo, zLhi, 'X', 'L' );
        Send( zLlo, zLhi, 'F', 'L' );

        Send( zRlo, zRhi, 'X', 'R' );
        Send( zRlo, zRhi, 'F', 'R' );
    }
    else {

        Recv( zihi + 1, zohi, 'X', 'R' );
        Recv( zihi + 1, zohi, 'F', 'R' );

        Recv( zolo, zilo - 1, 'X', 'L' );
        Recv( zolo, zilo - 1, 'F', 'L' );
    }

// Even send

    if( !(wkid & 1) ) {

        Send( zLlo, zLhi, 'X', 'L' );
        Send( zLlo, zLhi, 'F', 'L' );

        Send( zRlo, zRhi, 'X', 'R' );
        Send( zRlo, zRhi, 'F', 'R' );
    }
    else {

        Recv( zihi + 1, zohi, 'X', 'R' );
        Recv( zihi + 1, zohi, 'F', 'R' );

        Recv( zolo, zilo - 1, 'X', 'L' );
        Recv( zolo, zilo - 1, 'F', 'L' );
    }

    return true;
}


