

#include	"File.h"
#include	"PipeFiles.h"

#include	<string.h>

using namespace ns_pipergns;


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* AFromIDB ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool Rgns::AFromIDB()
{
    vector<Til2Img>	t2i;

    // Get rgn #1 tforms

    if( !IDBT2IGet_JustIDandT( t2i, *idb, z ) )
        return false;

    x.resize( nr * 6 );

    int						nt = t2i.size();
    map<int,int>::iterator	en = m.end();

    // For each transform in IDB...

    for( int it = 0; it < nt; ++it ) {

        // Get its block start and limit {j0,jlim}

        const Til2Img&			T = t2i[it];
        map<int,int>::iterator	mi = m.find( T.id );
        int						j0, jlim;

        if( mi == en )
            continue;

        j0		= mi->second;
        jlim	= (++mi != en ? mi->second : nr);

        // Propagate rgn #1 tform to all block members

        for( int j = j0; j < jlim; ++j ) {

            T.T.CopyOut( &x[j * 6] );
            FLAG_SETUSED( flag[j] );
        }
    }

    return true;
}

/* --------------------------------------------------------------- */
/* AFromTxt ------------------------------------------------------ */
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


static bool Read_vA(
    vector<CIDRA>	&vA,
    const char		*path,
    int				z,
    FILE			*flog )
{
    char	buf[2048];
    FILE	*f;
    CIDRA	A;
    bool	nf = true;	// default = no folds

    vA.clear();

    sprintf( buf, "%s/X_A_%d.txt", path, z );
    f = fopen( buf, "r" );

    if( f ) {

        while( A.FromFile( f ) ) {

            if( A.r > 1 )
                nf = false;

            vA.push_back( A );
        }

        fclose( f );
    }
    else
        fprintf( flog, "Rgns: Can't open [%s].\n", buf );

    return nf;
}


bool Rgns::AFromTxt( const char *path )
{
    vector<CIDRA>	vA;
    int				nf = Read_vA( vA, path, z, flog );

    x.resize( nr * 6 );

    int						na = vA.size();
    map<int,int>::iterator	en = m.end();

    if( !nf ) {	// Propagate rgn #1 to all block members

        // For each transform in vA...

        for( int ia = 0; ia < na; ++ia ) {

            // Get its block start and limit {j0,jlim}

            const CIDRA&			A = vA[ia];
            map<int,int>::iterator	mi = m.find( A.id );
            int						j0, jlim;

            if( mi == en )
                continue;

            j0		= mi->second;
            jlim	= (++mi != en ? mi->second : nr);

            // Propagate rgn #1 tform to all block members

            for( int j = j0; j < jlim; ++j ) {

                A.A.CopyOut( &x[j * 6] );
                FLAG_SETUSED( flag[j] );
            }
        }
    }
    else {	// Move each A to its specified position

        // For each transform in vA...

        for( int ia = 0; ia < na; ++ia ) {

            // Get its location j

            const CIDRA&			A = vA[ia];
            map<int,int>::iterator	mi = m.find( A.id );
            int						j;

            if( mi == en )
                continue;

            j = mi->second + A.r - 1;

            A.A.CopyOut( &x[j * 6] );
            FLAG_SETUSED( flag[j] );
        }
    }

    return true;
}

/* --------------------------------------------------------------- */
/* HFromTxt ------------------------------------------------------ */
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


static bool Read_vH(
    vector<CIDRH>	&vH,
    const char		*path,
    int				z,
    FILE*			flog )
{
    char	buf[2048];
    FILE	*f;
    CIDRH	H;
    bool	nf = true;	// default = no folds

    vH.clear();

    sprintf( buf, "%s/X_H_%d.txt", path, z );
    f = fopen( buf, "r" );

    if( f ) {

        while( H.FromFile( f ) ) {

            if( H.r > 1 )
                nf = false;

            vH.push_back( H );
        }

        fclose( f );
    }
    else
        fprintf( flog, "Rgns: Can't open [%s].\n", buf );

    return nf;
}


bool Rgns::HFromTxt( const char *path )
{
    vector<CIDRH>	vH;
    int				nf = Read_vH( vH, path, z, flog );

    x.resize( nr * 8 );

    int						nh = vH.size();
    map<int,int>::iterator	en = m.end();

    if( !nf ) {	// Propagate rgn #1 to all block members

        // For each transform in vH...

        for( int ih = 0; ih < nh; ++ih ) {

            // Get its block start and limit {j0,jlim}

            const CIDRH&			H = vH[ih];
            map<int,int>::iterator	mi = m.find( H.id );
            int						j0, jlim;

            if( mi == en )
                continue;

            j0		= mi->second;
            jlim	= (++mi != en ? mi->second : nr);

            // Propagate rgn #1 tform to all block members

            for( int j = j0; j < jlim; ++j ) {

                H.H.CopyOut( &x[j * 8] );
                FLAG_SETUSED( flag[j] );
            }
        }
    }
    else {	// Move each H to its specified position

        // For each transform in vH...

        for( int ih = 0; ih < nh; ++ih ) {

            // Get its location j

            const CIDRH&			H = vH[ih];
            map<int,int>::iterator	mi = m.find( H.id );
            int						j;

            if( mi == en )
                continue;

            j = mi->second + H.r - 1;

            H.H.CopyOut( &x[j * 8] );
            FLAG_SETUSED( flag[j] );
        }
    }

    return true;
}

/* --------------------------------------------------------------- */
/* Rgns::ReadXBin ------------------------------------------------ */
/* --------------------------------------------------------------- */

bool Rgns::ReadXBin( const char *path )
{
    int	nx = nr * NE;

    if( !nx )
        return false;

    char	buf[2048];
    FILE	*f;
    sprintf( buf, "%s/X_%c_%d.bin", path, (NE == 6 ? 'A' : 'H'), z );

    if( f = fopen( buf, "rb" ) ) {

        x.resize( nx );
        fread( &x[0], sizeof(double), nx, f );
        fclose( f );
        return true;
    }
    else {
        fprintf( flog, "Rgns: Can't open [%s].\n", buf );
        return false;
    }
}

/* --------------------------------------------------------------- */
/* Rgns::ReadFBin ------------------------------------------------ */
/* --------------------------------------------------------------- */

bool Rgns::ReadFBin( const char *path )
{
    if( !nr )
        return false;

    char	buf[2048];
    FILE	*f;
    sprintf( buf, "%s/F_%d.bin", path, z );

    if( f = fopen( buf, "rb" ) ) {

        fread( &flag[0], sizeof(uint8), nr, f );
        fclose( f );
        return true;
    }
    else {
        fprintf( flog, "Rgns: Can't open [%s].\n", buf );
        return false;
    }
}

/* --------------------------------------------------------------- */
/* Rgns::Init ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Rgns::Init( const string &idb, int iz, FILE *flog )
{
    this->flog	= flog;
    this->idb	= &idb;
    z			= iz;
    nr			= IDBGetIDRgnMap( m, idb, z, flog );

    flag.assign( nr, fmRead );
    return (nr != 0);
}

/* --------------------------------------------------------------- */
/* Rgns::Load ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Rgns::Load( const char *path )
{
    if( !path || !path[0] ) {

        NE = 6;
        return AFromIDB();
    }
    else {

        const char	*name = FileNamePtr( path );

        if( strstr( name, "X_A" ) ) {

            NE = 6;

            if( strstr( name, "X_A_TXT" ) ||
                strstr( name, "X_A_MET" ) ) {

                return AFromTxt( path );
            }
            else if( strstr( name, "X_A_BIN" ) )
                return ReadXBin( path ) && ReadFBin( path );
            else
                goto error;
        }
        else if( strstr( name, "X_H" ) ) {

            NE = 8;

            if( strstr( name, "X_H_TXT" ) ||
                strstr( name, "X_H_MET" ) ) {

                return HFromTxt( path );
            }
            else if( strstr( name, "X_H_BIN" ) )
                return ReadXBin( path ) && ReadFBin( path );
            else
                goto error;
        }
        else {
error:
            fprintf( flog,
            "Rgns: Unknown folder name pattern [%s].\n", name );
            exit( 42 );
        }
    }
}

/* --------------------------------------------------------------- */
/* Rgns::SaveBIN ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Save currently loaded tforms (and optionally flags)
// in binary format.
//
// Caller must create containing folder.
//
// Return true if success.
//
bool Rgns::SaveBIN( const char *path, bool writeflags )
{
    const char	*name	= FileNamePtr( path );
    const char	*isaff	= strstr( name, "X_A_BIN" );
    bool		ok		= false;

    if( isaff || strstr( name, "X_H_BIN" ) ) {

        if( (isaff && NE != 6) || (!isaff && NE != 8) ) {
            fprintf( flog,
            "Rgns: SaveBIN type mismatch NE=%d -> [%s].\n",
            NE, name );
        }
        else {

            char	buf[2048];
            FILE	*f;

            sprintf( buf, "%s/X_%c_%d.bin",
                path, (NE == 6 ? 'A' : 'H'), z );

            f = FileOpenOrDie( buf, "wb" );
            fwrite( &x[0], sizeof(double), x.size(), f );
            fclose( f );

            if( writeflags ) {

                sprintf( buf, "%s/F_%d.bin", path, z );
                f = FileOpenOrDie( buf, "wb" );
                fwrite( &flag[0], sizeof(uint8), nr, f );
                fclose( f );
            }

            ok = true;
        }
    }
    else
        fprintf( flog, "Rgns: Unknown SaveBIN type [%s].\n", name );

    return ok;
}

/* --------------------------------------------------------------- */
/* Rgns::SaveTXT ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Save currently loaded tforms as text. Text-style folders
// never contain flag files.
//
// Caller must create containing folder.
//
// Return true if success.
//
bool Rgns::SaveTXT( const char *path )
{
    const char	*name = FileNamePtr( path );
    FILE		*f;
    char		buf[2048];
    bool		ok = false;

    if( strstr( name, "X_A_TXT" ) ) {

        if( NE != 6 ) {
            fprintf( flog,
            "Rgns: SaveTXT type mismatch NE=%d -> [%s].\n",
            NE, name );
        }
        else {

            sprintf( buf, "%s/X_A_%d.txt", path, z );
            f = FileOpenOrDie( buf, "w", flog );

            map<int,int>::iterator	mi, en = m.end();

            for( mi = m.begin(); mi != en; ) {

                int	id		= mi->first,
                    j0		= mi->second,
                    jlim	= (++mi == en ? nr : mi->second);

                for( int j = j0; j < jlim; ++j ) {

                    if( !FLAG_ISUSED( flag[j] ) )
                        continue;

                    TAffine&	T = X_AS_AFF( x, j );

                    fprintf( f, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
                    id, j - j0 + 1,
                    T.t[0], T.t[1], T.t[2],
                    T.t[3], T.t[4], T.t[5] );
                }
            }

            fclose( f );
            ok = true;
        }
    }
    else if( strstr( name, "X_H_TXT" ) ) {

        if( NE != 8 ) {
            fprintf( flog,
            "Rgns: SaveTXT type mismatch NE=%d -> [%s].\n",
            NE, name );
        }
        else {

            sprintf( buf, "%s/X_H_%d.txt", path, z );
            f = FileOpenOrDie( buf, "w", flog );

            map<int,int>::iterator	mi, en = m.end();

            for( mi = m.begin(); mi != en; ) {

                int	id		= mi->first,
                    j0		= mi->second,
                    jlim	= (++mi == en ? nr : mi->second);

                for( int j = j0; j < jlim; ++j ) {

                    if( !FLAG_ISUSED( flag[j] ) )
                        continue;

                    THmgphy&	T = X_AS_HMY( x, j );

                    fprintf( f,
                    "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%.12g\t%.12g\n",
                    id, j - j0 + 1,
                    T.t[0], T.t[1], T.t[2],
                    T.t[3], T.t[4], T.t[5],
                    T.t[6], T.t[7] );
                }
            }

            fclose( f );
            ok = true;
        }
    }
    else
        fprintf( flog, "Rgns: Unknown SaveTXT type [%s].\n", name );

    return ok;
}


