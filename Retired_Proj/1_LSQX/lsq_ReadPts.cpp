

#include	"lsq_ReadPts.h"

#include	"File.h"
#include	"Maths.h"

#include	<string.h>

//#include	<unordered_map>
//using namespace std;


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

map<MZID,int>	nConRgn;	// # connected-rgn this image
int				gW = 4056,	// image coords
                gH = 4056;






/* --------------------------------------------------------------- */
/* CMapZIDR ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// URL <http://www.drdobbs.com/windows/user-defined-hash-functions-for-unordere/231600210?pgno=3>

struct Hasher {
    size_t operator()( const MZIDR &m ) const
    {
        return SuperFastHash( (char*)&m.z, 3*sizeof(int) );
    };
};

class CMapZIDR {
private:
    map<MZIDR,int>	mi;
//	unordered_map<MZIDR,int,Hasher>	mi;
    int								nr;
public:
    CMapZIDR() : nr(0) {};
    int Find( const RGN &R );
};

// Return zero-based index for given RGN.
//
// If already stored, that index is returned. Else, new
// entries are created with index (nr); nr is incremented.
//
int CMapZIDR::Find( const RGN &R )
{
    MZIDR						key( R.z, R.id, R.rgn );
    map<MZIDR,int>::iterator	it = mi.find( key );

    if( it != mi.end() )
        return it->second;

    mi[key] = nr;
    vRgn.push_back( R );
    return nr++;
}

/* --------------------------------------------------------------- */
/* ReadPts_StrTags ----------------------------------------------- */
/* --------------------------------------------------------------- */

void ReadPts_StrTags(
    FILE		*FOUT,
    CNX			*cnx,
    SML			*sml,
    int			(*IDFromName)( const char *name ),
    const char	*dirfile,
    const char	*ptsfile )
{
    printf( "---- Read pts ----\n" );

    FILE		*f = FileOpenOrDie( ptsfile, "r" );
    CLineScan	LS;
    DIR			dir;	// map name strings to z layers

    dir.ReadDIRFile( dirfile, FOUT );

    CMapZIDR	CM;
    int			nlines = 0;

    for(;;) {

        char	name1[2048], name2[2048];

        if( LS.Get( f ) <= 0 )
            break;

        ++nlines;

        if( !strncmp( LS.line, "CPOINT", 6 ) ) {

            char	key1[32], key2[32];
            Point	p1, p2;

            if( 8 != sscanf( LS.line + 7,
                        "'%[^']' %s %lf %lf '%[^']' %s %lf %lf",
                        name1, key1, &p1.x, &p1.y,
                        name2, key2, &p2.x, &p2.y ) ) {

                printf(
                "WARNING: 'CPOINT' format error; line %d\n",
                nlines );

                continue;
            }

            RGN	R1( name1, key1 );
            RGN	R2( name2, key2 );
            int r1 = CM.Find( R1 );
            int r2 = CM.Find( R2 );

            cnx->AddCorrespondence( r1, r2 );
            sml->AddPOINTPair( r1, p1, r2, p2 );

            vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
        }
        else if( !strncmp( LS.line, "POINT", 5 ) ) {

            Point	p1, p2;

            if( 6 != sscanf( LS.line + 6,
                        "%s %lf %lf %s %lf %lf",
                        name1, &p1.x, &p1.y,
                        name2, &p2.x, &p2.y ) ) {

                printf(
                "WARNING: 'POINT' format error; line %d\n",
                nlines );

                continue;
            }

            RGN	R1( name1, dir, IDFromName( name1 ) );
            RGN	R2( name2, dir, IDFromName( name2 ) );
            int r1 = CM.Find( R1 );
            int r2 = CM.Find( R2 );

            cnx->AddCorrespondence( r1, r2 );
            sml->AddPOINTPair( r1, p1, r2, p2 );

            vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
        }
        else if( !strncmp( LS.line, "FOLDMAP", 7 ) ) {

            int	z, id, nrgn = -1;

            sscanf( LS.line + 8, "'%*[^']' %s %d", name1, &nrgn );
            ZIDFromFMPath( z, id, name1 );
            nConRgn[MZID( z, id )] = nrgn;

            fprintf( FOUT, LS.line );
        }
        else if( !strncmp( LS.line, "IMAGESIZE", 9 ) ) {

            if( 2 != sscanf( LS.line + 10, "%d %d", &gW, &gH ) ) {
                printf( "Bad IMAGESIZE line '%s'.\n", LS.line );
                exit( 42 );
            }

            fprintf( FOUT, LS.line );
            printf( LS.line );
        }
        else {

            char	*s = strchr( LS.line, ' ' );

            if( s )
                *s = 0;

            printf(
            "WARNING: Unknown entry type; '%s' line %d\n",
            LS.line, nlines );
        }
    }

    fclose( f );

    printf( "\n" );
}

/* --------------------------------------------------------------- */
/* IsDaviCorner -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true to reject corner-corner matches.
//
static bool IsDaviCorner( const RGN &R1, const RGN &R2 )
{
    if( R1.z != R2.z )
        return false;

    const Til2Img	*m1, *m2;
    RGN::GetMeta( &m1, &m2, R1, R2 );

    return ((m2->row - m1->row) * (m2->col - m1->col) != 0);
}

/* --------------------------------------------------------------- */
/* RejectPair ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Chance to optionally apply rejection criteria against
// point pairs.
//
// Return true to reject.
//
static bool RejectPair( const RGN &R1, const RGN &R2 )
{
// ------------------------------------
// reject layer 30 and 38
#if 0
    if( R1.z == 30 || R1.z == 38 || R2.z == 30 || R2.z == 38 )
        return true;
#endif
// ------------------------------------

// ------------------------------------
// accept only col/row subsection
#if 0
    const Til2Img	*m1, *m2;
    RGN::GetMeta( &m1, &m2, R1, R2 );
    int		col = m1->col, row = m1->row;

    if( col != -999 ) {

//		if( col < 33 || col > 39 || row < 85 || row > 91 )
        if( col < 25-4 || col > 25+4 || row < 19-4 || row > 19+4 )
//		if( row > 2 || col > 2 )
            return true;

//		if( col == 44 && row == 19 )
//		if( col == 48 && row == 16 )
//			return true;

// repeat tests for R2
        col = p2->col, row = p2->row;

//		if( col < 33 || col > 39 || row < 85 || row > 91 )
        if( col < 25-4 || col > 25+4 || row < 19-4 || row > 19+4 )
//		if( row > 2 || col > 2 )
            return true;

//		if( col == 44 && row == 19 )
//		if( col == 48 && row == 16 )
//			return true;
    }
#endif
// ------------------------------------

    return false;
}

/* --------------------------------------------------------------- */
/* ReadPts_NumTags ----------------------------------------------- */
/* --------------------------------------------------------------- */

void ReadPts_NumTags(
    FILE		*FOUT,
    CNX			*cnx,
    SML			*sml,
    const char	*ptsfile,
    int			davinocorn )
{
    printf( "---- Read pts ----\n" );

    FILE		*f = FileOpenOrDie( ptsfile, "r" );
    CLineScan	LS;

    CMapZIDR	CM;
    int			nlines = 0;

    for(;;) {

        char	name1[2048], name2[2048];

        if( LS.Get( f ) <= 0 )
            break;

        ++nlines;

        if( !strncmp( LS.line, "CPOINT2", 7 ) ) {

            char	key1[32], key2[32];
            Point	p1, p2;

            if( 6 != sscanf( LS.line + 8,
                        "%s %lf %lf %s %lf %lf",
                        key1, &p1.x, &p1.y,
                        key2, &p2.x, &p2.y ) ) {

                printf(
                "WARNING: 'CPOINT2' format error; line %d\n",
                nlines );

                continue;
            }

            RGN	R1( key1 );
            RGN	R2( key2 );

            if( davinocorn && IsDaviCorner( R1, R2 ) )
                continue;

            if( RejectPair( R1, R2 ) )
                continue;

            int r1 = CM.Find( R1 );
            int r2 = CM.Find( R2 );

            cnx->AddCorrespondence( r1, r2 );
//			sml->AddPOINTPair( r1, p1, r2, p2 );

            vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
        }
        else if( !strncmp( LS.line, "FOLDMAP2", 8 ) ) {

            int	z, id, nrgn;

            sscanf( LS.line + 9, "%d.%d %d", &z, &id, &nrgn );
            nConRgn[MZID( z, id )] = nrgn;
        }
        else if( !strncmp( LS.line, "IDBPATH", 7 ) ) {

            char	buf[2048];

            if( !sscanf( LS.line + 8, "%s", buf ) ) {

                printf( "Bad IDBPATH line '%s'.\n", LS.line );
                exit( 42 );
            }

            idb = buf;
            fprintf( FOUT, LS.line );
            printf( LS.line );
        }
        else if( !strncmp( LS.line, "IMAGESIZE", 9 ) ) {

            if( 2 != sscanf( LS.line + 10, "%d %d", &gW, &gH ) ) {
                printf( "Bad IMAGESIZE line '%s'.\n", LS.line );
                exit( 42 );
            }

            fprintf( FOUT, LS.line );
            printf( LS.line );
        }
        else {

            char	*s = strchr( LS.line, ' ' );

            if( s )
                *s = 0;

            printf(
            "WARNING: Unknown entry type; '%s' line %d\n",
            LS.line, nlines );
        }
    }

    fclose( f );

    printf( "\n" );

    IDBT2ICacheClear();
}


