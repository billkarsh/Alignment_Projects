

#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"

#include	<string.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static T2ICache	_C1,
                _C2;






/* --------------------------------------------------------------- */
/* OpenPairLog --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rename stdout using image labels.
//
void OpenPairLog( int alr, int atl, int blr, int btl )
{
    char	slog[256];

    sprintf( slog, "pair_%d.%d^%d.%d.log", alr, atl, blr, btl );
    freopen( slog, "a", stdout );
}

/* --------------------------------------------------------------- */
/* NamePtsFile --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Standard file capturing stdout from ptest.
//
char *NamePtsFile( char *buf, int alr, int blr )
{
    const char	*sud;

    if( alr < blr )
        sud = "up";
    else if( alr == blr )
        sud = "same";
    else
        sud = "down";

    sprintf( buf, "pts.%s", sud );

    return buf;
}

/* --------------------------------------------------------------- */
/* NameLogFile --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Standard file capturing stderr from ptest.
//
char *NameLogFile( char *buf, int alr, int atl, int blr, int btl )
{
    sprintf( buf, "pair_%d.%d^%d.%d.log", alr, atl, blr, btl );
    return buf;
}

/* --------------------------------------------------------------- */
/* GetPrm -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool GetPrm(
    void		*v,
    const char	*fmt,
    const char	*caller,
    FILE		*fin,
    FILE		*flog )
{
    CLineScan	LS;
    bool		ok = false;

    while( LS.Get( fin ) > 0 ) {

        if( isspace( LS.line[0] ) || LS.line[0] == '#' )
            continue;
        else if( 1 == sscanf( LS.line, fmt, v ) ) {

            fprintf( flog, LS.line );

            // rectify empty path spec
            if( strchr( fmt, 's' ) && ((char*)v)[0] == '<' )
                ((char*)v)[0] = 0;

            ok = true;
            goto exit;
        }
        else {
            fprintf( flog, "%s: Expected [%s], found [%.23s...].\n",
            caller, fmt, LS.line );
            goto exit;
        }
    }

    fprintf( flog, "%s: Can't find [%s].\n", caller, fmt );

exit:
    return ok;
}

/* --------------------------------------------------------------- */
/* ReadScriptParams ---------------------------------------------- */
/* --------------------------------------------------------------- */

#define	GETPRM_SCR( addr, pat )									\
    if( !GetPrm( addr, pat, "ReadScriptParams", f, flog ) )		\
        goto exit

// Read parameter file governing alignment pipeline.
//
bool ReadScriptParams(
    ScriptParams	&S,
    const char		*scriptpath,
    FILE			*flog )
{
    FILE	*f = fopen( scriptpath, "r" );
    int		ok = false;

    if( f ) {

        char	buf[2048];

        fprintf( flog, "\n---- Script parameters ----\n" );

        GETPRM_SCR( &S.slotspernode, "slotspernode=%d" );

        GETPRM_SCR( &S.usingfoldmasks, "usingfoldmasks=%c" );
        GETPRM_SCR( buf, "croprectfile=%s" ); S.croprectfile = buf;
        GETPRM_SCR( &S.makefmjparam, "makefmjparam=%d" );
        GETPRM_SCR( &S.makefmslots, "makefmslots=%d" );

        GETPRM_SCR( &S.createauxdirs, "createauxdirs=%c" );
        GETPRM_SCR( &S.montageblocksize, "montageblocksize=%d" );
        GETPRM_SCR( &S.mintileolapfrac, "mintileolapfrac=%lf" );
        GETPRM_SCR( &S.ignorecorners, "ignorecorners=%c" );
        GETPRM_SCR( &S.makesamejparam, "makesamejparam=%d" );
        GETPRM_SCR( &S.makesameslots, "makesameslots=%d" );

        GETPRM_SCR( &S.crossscale, "crossscale=%d" );
        GETPRM_SCR( &S.legendremaxorder, "legendremaxorder=%d" );
        GETPRM_SCR( &S.rendersdevcnts, "rendersdevcnts=%d" );
        GETPRM_SCR( &S.maskoutresin, "maskoutresin=%c" );

        GETPRM_SCR( &S.stripwidth, "stripwidth=%d" );
        GETPRM_SCR( &S.stripsweepspan, "stripsweepspan=%lf" );
        GETPRM_SCR( &S.stripsweepstep, "stripsweepstep=%lf" );
        GETPRM_SCR( &S.stripmincorr, "stripmincorr=%lf" );
        GETPRM_SCR( &S.stripslots, "stripslots=%d" );

        GETPRM_SCR( &S.crossblocksize, "crossblocksize=%d" );
        GETPRM_SCR( &S.blocksweepspan, "blocksweepspan=%lf" );
        GETPRM_SCR( &S.blocksweepstep, "blocksweepstep=%lf" );
        GETPRM_SCR( &S.blockxyconf, "blockxyconf=%lf" );
        GETPRM_SCR( &S.blockmincorr, "blockmincorr=%lf" );
        GETPRM_SCR( &S.blocknomcorr, "blocknomcorr=%lf" );
        GETPRM_SCR( &S.blocknomcoverage, "blocknomcoverage=%lf" );
        GETPRM_SCR( &S.blockreqdz, "blockreqdz=%d" );
        GETPRM_SCR( &S.blockmaxdz, "blockmaxdz=%d" );
        GETPRM_SCR( &S.blockslots, "blockslots=%d" );

        GETPRM_SCR( &S.makedownjparam, "makedownjparam=%d" );
        GETPRM_SCR( &S.makedownslots, "makedownslots=%d" );

        GETPRM_SCR( &S.xmlpixeltype, "xmlpixeltype=%d" );
        GETPRM_SCR( &S.xmlsclmin, "xmlsclmin=%d" );
        GETPRM_SCR( &S.xmlsclmax, "xmlsclmax=%d" );

        // finish Y/N booleans
        S.usingfoldmasks	= (toupper( S.usingfoldmasks ) == 'Y');
        S.createauxdirs		= (toupper( S.createauxdirs ) == 'Y');
        S.ignorecorners		= (toupper( S.ignorecorners ) == 'Y');
        S.maskoutresin		= (toupper( S.maskoutresin ) == 'Y');

        fprintf( flog, "\n" );

        ok = true;
    }
    else {
        fprintf( flog,
        "ReadScriptParams: Can't open [%s].\n", scriptpath );
    }

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* ReadMatchParams ----------------------------------------------- */
/* --------------------------------------------------------------- */

#define	GETPRM_MCH( addr, pat )									\
    if( !GetPrm( addr, pat, "ReadMatchParams", f, flog ) )		\
        goto exit

// Read parameter file governing thumbnail and mesh matching.
//
// Path precedence:
//
// 1. Use explicit matchparamspath if given.
// 2. Use (alr,blr) params to try loading override file from
//		temp dir with name pattern 'matchparams_alr_blr.txt'.
// 3. Load the standard temp/matchparams.txt.
//
bool ReadMatchParams(
    MatchParams		&M,
    int				alr,
    int				blr,
    const char		*matchparamspath,
    FILE			*flog )
{
    char	name[1024];
    FILE	*f = NULL;
    int		ok = false;

    if( matchparamspath )
        strcpy( name, matchparamspath );
    else {

        sprintf( name, "../../matchparams_%d_%d.txt", alr, blr );
        f = fopen( name, "r" );

        if( !f )
            sprintf( name, "../../matchparams.txt" );
    }

    if( f || (f = fopen( name, "r" )) ) {

        fprintf( flog, "\n---- Match parameters ----\n" );

        GETPRM_MCH( &M.PXBRO, "PXBRO=%d" );
        GETPRM_MCH( &M.PXLENS, "PXLENS=%c" );
        GETPRM_MCH( &M.PXRESMSK, "PXRESMSK=%c" );
        GETPRM_MCH( &M.PXDOG, "PXDOG=%c" );
        GETPRM_MCH( &M.PXDOG_R1, "PXDOG_R1=%d" );
        GETPRM_MCH( &M.PXDOG_R2, "PXDOG_R2=%d" );
        GETPRM_MCH( &M.FLD, "FLD=%c" );
        GETPRM_MCH( &M.PRETWEAK, "PRETWEAK=%c" );
        GETPRM_MCH( &M.SCALE, "SCALE=%lf" );
        GETPRM_MCH( &M.XSCALE, "XSCALE=%lf" );
        GETPRM_MCH( &M.YSCALE, "YSCALE=%lf" );
        GETPRM_MCH( &M.SKEW, "SKEW=%lf" );

        GETPRM_MCH( &M.MODE_SL, "MODE_SL=%c" );
        GETPRM_MCH( &M.MODE_XL, "MODE_XL=%c" );
        GETPRM_MCH( &M.TAB2DFM_SL, "TAB2DFM_SL=%c" );
        GETPRM_MCH( &M.TAB2DFM_XL, "TAB2DFM_XL=%c" );
        GETPRM_MCH( &M.XYCONF_SL, "XYCONF_SL=%lf" );
        GETPRM_MCH( &M.XYCONF_XL, "XYCONF_XL=%lf" );
        GETPRM_MCH( &M.THMDEC_SL, "THMDEC_SL=%d" );
        GETPRM_MCH( &M.THMDEC_XL, "THMDEC_XL=%d" );
        GETPRM_MCH( &M.OLAP1D_SL, "OLAP1D_SL=%d" );
        GETPRM_MCH( &M.OLAP1D_XL, "OLAP1D_XL=%d" );
        GETPRM_MCH( &M.OLAP2D_SL, "OLAP2D_SL=%d" );
        GETPRM_MCH( &M.OLAP2D_XL, "OLAP2D_XL=%d" );
        GETPRM_MCH( &M.NBMXHT_SL, "NBMXHT_SL=%lf" );
        GETPRM_MCH( &M.NBMXHT_XL, "NBMXHT_XL=%lf" );
        GETPRM_MCH( &M.HFANGDN_SL, "HFANGDN_SL=%lf" );
        GETPRM_MCH( &M.HFANGDN_XL, "HFANGDN_XL=%lf" );
        GETPRM_MCH( &M.HFANGPR_SL, "HFANGPR_SL=%lf" );
        GETPRM_MCH( &M.HFANGPR_XL, "HFANGPR_XL=%lf" );
        GETPRM_MCH( &M.RTRSH_SL, "RTRSH_SL=%lf" );
        GETPRM_MCH( &M.RTRSH_XL, "RTRSH_XL=%lf" );
        GETPRM_MCH( &M.TWEAKS, "TWEAKS=%c" );
        GETPRM_MCH( &M.LIMXY_SL, "LIMXY_SL=%d" );
        GETPRM_MCH( &M.LIMXY_XL, "LIMXY_XL=%d" );
        GETPRM_MCH( &M.WTHMPR, "WTHMPR=%c" );

        GETPRM_MCH( &M.OPT_SL, "OPT_SL=%c" );
        GETPRM_MCH( &M.RIT_SL, "RIT_SL=%lf" );
        GETPRM_MCH( &M.RIT_XL, "RIT_XL=%lf" );
        GETPRM_MCH( &M.RFA_SL, "RFA_SL=%lf" );
        GETPRM_MCH( &M.RFA_XL, "RFA_XL=%lf" );
        GETPRM_MCH( &M.RFT_SL, "RFT_SL=%lf" );
        GETPRM_MCH( &M.RFT_XL, "RFT_XL=%lf" );
        GETPRM_MCH( &M.TMC, "TMC=%lf" );
        GETPRM_MCH( &M.TSC, "TSC=%lf" );
        GETPRM_MCH( &M.MNL, "MNL=%d" );
        GETPRM_MCH( &M.MTA, "MTA=%d" );
        GETPRM_MCH( &M.MMA, "MMA=%d" );
        GETPRM_MCH( &M.ONE, "ONE=%c" );
        GETPRM_MCH( &M.IFM, "IFM=%lf" );
        GETPRM_MCH( &M.FFM, "FFM=%lf" );
        GETPRM_MCH( &M.FYL, "FYL=%lf" );
        GETPRM_MCH( &M.CPD, "CPD=%lf" );
        GETPRM_MCH( &M.EMM, "EMM=%c" );
        GETPRM_MCH( &M.EMT, "EMT=%lf" );
        GETPRM_MCH( &M.WDI, "WDI=%c" );
        GETPRM_MCH( &M.LDA, "LDA=%lf" );
        GETPRM_MCH( &M.LDR, "LDR=%lf" );
        GETPRM_MCH( &M.LDC, "LDC=%lf" );
        GETPRM_MCH( &M.DXY, "DXY=%lf" );
        GETPRM_MCH( &M.WMT, "WMT=%c" );
        GETPRM_MCH( &M.WTT, "WTT=%c" );

        // ensure upper case
        M.FLD		= toupper( M.FLD );
        M.MODE_SL	= toupper( M.MODE_SL );
        M.MODE_XL	= toupper( M.MODE_XL );

        // finish Y/N booleans
        M.PXLENS		= (toupper( M.PXLENS ) == 'Y');
        M.PXRESMSK		= (toupper( M.PXRESMSK ) == 'Y');
        M.PXDOG			= (toupper( M.PXDOG ) == 'Y');
        M.PRETWEAK		= (toupper( M.PRETWEAK ) == 'Y');
        M.TAB2DFM_SL	= (toupper( M.TAB2DFM_SL ) == 'Y');
        M.TAB2DFM_XL	= (toupper( M.TAB2DFM_XL ) == 'Y');
        M.TWEAKS		= (toupper( M.TWEAKS ) == 'Y');
        M.WTHMPR		= (toupper( M.WTHMPR ) == 'Y');
        M.OPT_SL		= (toupper( M.OPT_SL ) == 'Y');
        M.ONE			= (toupper( M.ONE ) == 'Y');
        M.EMM			= (toupper( M.EMM ) == 'Y');
        M.WDI			= (toupper( M.WDI ) == 'Y');
        M.WMT			= (toupper( M.WMT ) == 'Y');
        M.WTT			= (toupper( M.WTT ) == 'Y');

        fprintf( flog, "\n" );

        ok = true;
    }
    else
        fprintf( flog, "ReadMatchParams: Can't open [%s].\n", name );

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBFromTemp --------------------------------------------------- */
/* --------------------------------------------------------------- */

void IDBFromTemp(
    string		&idbpath,
    const char	*tempdir,
    FILE		*flog )
{
    idbpath.clear();

    char	buf[2048];
    sprintf( buf, "%s/imageparams.txt", tempdir );
    FILE	*f = fopen( buf, "r" );

    if( f ) {

        CLineScan	LS;

        while( LS.Get( f ) > 0 ) {

            if( 1 == sscanf( LS.line, "IDBPATH %[^\n]", buf ) ) {

                idbpath = buf;
                goto close;
            }
        }

        fprintf( flog,
        "IDBFromTemp: imageparams.txt missing IDBPATH tag.\n" );

close:
        fclose( f );
    }
    else {
        fprintf( flog,
        "IDBFromTemp: Can't open imageparams.txt.\n" );
    }
}

/* --------------------------------------------------------------- */
/* IDBGetImageDims ----------------------------------------------- */
/* --------------------------------------------------------------- */

bool IDBGetImageDims(
    int				&w,
    int				&h,
    const string	&idb,
    FILE			*flog )
{
    char	name[2048];
    int		ok = false;
    sprintf( name, "%s/imageparams.txt", idb.c_str() );
    FILE	*f = fopen( name, "r" );

    if( f ) {

        CLineScan	LS;

        while( LS.Get( f ) > 0 ) {

            if( 2 == sscanf( LS.line, "IMAGESIZE %d %d", &w, &h ) ) {
                ok = true;
                goto close;
            }
        }

        fprintf( flog,
        "IDBGetImageDims: imageparams.txt missing IMAGESIZE tag.\n" );

close:
        fclose( f );
    }
    else {
        fprintf( flog,
        "IDBGetImageDims: Can't open [%s].\n", name );
    }

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBGetIDRgnMap ------------------------------------------------ */
/* --------------------------------------------------------------- */

class Cfmline {
public:
    int	id, nr;
public:
    Cfmline()	{};
    inline bool FromFile( FILE *f )
    {
        int	z;
        return 3 == fscanf( f, "FOLDMAP2"
            " %d.%d %d\n", &z, &id, &nr );
    };
    bool operator < ( const Cfmline &rhs ) const
    {
        return id < rhs.id;
    };
};

// For given z, create mapping from (tileID,rgn) to unique
// zero based array index. The caller uses the map<> to get
// an index for (id,rgn) as follows:
//
// --- create the map
// map<int,int>	idmap;
// IDBGetIDRgnMap( idmap, idb, z );
//
// --- use map to get index for (id,rgn)
// map<int,int>::iterator	it = idmap.find( id );
// int						index = it->second + rgn - 1;
//
// Return array size.
//
// Notes:
// The entries are sorted before mapping, so one can recover
// nr from a map iterator 'it': ((it+1)->second - it->second).
//
int IDBGetIDRgnMap(
    map<int,int>	&m,
    const string	&idb,
    int				z,
    FILE			*flog )
{
    char			name[2048];
    FILE			*f;
    vector<Cfmline>	vline;
    int				nelem = 0;

    m.clear();

// Scan entries into vector

    sprintf( name, "%s/%d/fm.same", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        Cfmline	line;

        while( line.FromFile( f ) )
            vline.push_back( line );
    }
    else {

        if( flog ) {
            fprintf( flog,
            "IDBGetIDRgnMap: Can't open [%s].\n", name );
        }

        return 0;
    }

    if( f )
        fclose( f );

// Sort entries by id, then map to 0-based indices

    int	n = vline.size();

    if( n ) {
        sort( vline.begin(), vline.end() );
        for( int i = 0; i < n; ++i ) {
            m[vline[i].id] = nelem;
            nelem         += vline[i].nr;
        }
    }
    else if( flog )
        fprintf( flog, "IDBGetIDRgnMap: Empty file [%s].\n", name );

    return nelem;
}

/* --------------------------------------------------------------- */
/* IDBT2IGet1 ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Fetch t2i clone for this z/id.
//
bool IDBT2IGet1(
    Til2Img			&t2i,
    const string	&idb,
    int				z,
    int				id,
    const char		*forcepath,
    FILE			*flog )
{
// override provided

    if( forcepath ) {
        t2i.id		= id;
        t2i.T		= TAffine( 1,0,0,0,1,0 );
        t2i.col		= -999;
        t2i.row		= -999;
        t2i.cam		= 0;
        t2i.path	= forcepath;
        fprintf( flog, "IDBT2IGet1: forcepath = [%s].\n", forcepath );
        return true;
    }

// from idb

    char	name[2048];
    FILE	*f;
    int		ok = false;

    if( idb.empty() )
        sprintf( name, "../%d/TileToImage.txt", z );
    else
        sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 ) {
            fprintf( flog, "IDBT2IGet1: Empty file [%s].\n", name );
            goto exit;
        }

        while( LS.Get( f ) > 0 ) {

            char	buf[2048];

            t2i.id = -1;

            sscanf( LS.line, "%d", &t2i.id );

            if( t2i.id != id )
                continue;

            sscanf( LS.line,
            "%d"
            "\t%lf\t%lf\t%lf"
            "\t%lf\t%lf\t%lf"
            "\t%d\t%d\t%d"
            "\t%[^\t\n]",
            &t2i.id,
            &t2i.T.t[0], &t2i.T.t[1], &t2i.T.t[2],
            &t2i.T.t[3], &t2i.T.t[4], &t2i.T.t[5],
            &t2i.col, &t2i.row, &t2i.cam,
            buf );

            t2i.path	= buf;
            ok			= true;
            goto exit;
        }

        fprintf( flog, "IDBT2IGet1: No entry for [%d %d].\n", z, id );
    }
    else
        fprintf( flog, "IDBT2IGet1: Can't open [%s].\n", name );

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBT2IGetAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Get vector of clones.
//
bool IDBT2IGetAll(
    vector<Til2Img>	&t2i,
    const string	&idb,
    int				z,
    FILE			*flog )
{
    char	name[2048];
    FILE	*f;
    int		ok = false;

    t2i.clear();

    if( idb.empty() )
        sprintf( name, "../%d/TileToImage.txt", z );
    else
        sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 ) {
            fprintf( flog,
            "IDBT2IGetAll: Empty file [%s].\n", name );
            goto exit;
        }

        while( LS.Get( f ) > 0 ) {

            Til2Img	E;
            char	buf[2048];

            sscanf( LS.line,
            "%d"
            "\t%lf\t%lf\t%lf"
            "\t%lf\t%lf\t%lf"
            "\t%d\t%d\t%d"
            "\t%[^\t\n]",
            &E.id,
            &E.T.t[0], &E.T.t[1], &E.T.t[2],
            &E.T.t[3], &E.T.t[4], &E.T.t[5],
            &E.col, &E.row, &E.cam,
            buf );

            E.path = buf;

            t2i.push_back( E );
        }

        ok = true;
    }
    else
        fprintf( flog, "IDBT2IGetAll: Can't open [%s].\n", name );

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBT2IGet_JustIDandT ------------------------------------------ */
/* --------------------------------------------------------------- */

// Get uncached vector of Til2Img with valid {id,T}.
//
// Entries happen to be ordered by id because idb files are
// always maintained that way.
//
bool IDBT2IGet_JustIDandT(
    vector<Til2Img>	&t2i,
    const string	&idb,
    int				z,
    FILE			*flog )
{
    char	name[2048];
    FILE	*f;
    int		ok = false;

    t2i.clear();

    sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 ) {
            fprintf( flog,
            "IDBT2IGet_JustIDandT: Empty file [%s].\n", name );
            goto exit;
        }

        while( LS.Get( f ) > 0 ) {

            Til2Img	E;

            sscanf( LS.line,
            "%d"
            "\t%lf\t%lf\t%lf"
            "\t%lf\t%lf\t%lf",
            &E.id,
            &E.T.t[0], &E.T.t[1], &E.T.t[2],
            &E.T.t[3], &E.T.t[4], &E.T.t[5] );

            t2i.push_back( E );
        }

        ok = true;
    }
    else {
        fprintf( flog,
        "IDBT2IGet_JustIDandT: Can't open [%s].\n", name );
    }

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheLoad ----------------------------------------------- */
/* --------------------------------------------------------------- */

bool IDBT2ICacheLoad(
    T2ICache		&C,
    const string	&idb,
    int				z,
    FILE			*flog )
{
    char	name[2048];
    FILE	*f;

    IDBT2ICacheClear( C );

    if( idb.empty() )
        sprintf( name, "../%d/TileToImage.txt", z );
    else
        sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 ) {
            fprintf( flog,
            "IDBT2ICacheLoad: Empty file [%s].\n", name );
            goto exit;
        }

        while( LS.Get( f ) > 0 ) {

            Til2Img	E;
            char	buf[2048];

            sscanf( LS.line,
            "%d"
            "\t%lf\t%lf\t%lf"
            "\t%lf\t%lf\t%lf"
            "\t%d\t%d\t%d"
            "\t%[^\t\n]",
            &E.id,
            &E.T.t[0], &E.T.t[1], &E.T.t[2],
            &E.T.t[3], &E.T.t[4], &E.T.t[5],
            &E.col, &E.row, &E.cam,
            buf );

            E.path = buf;

            C.m[E.id] = E;
        }

        C.z = z;
    }
    else
        fprintf( flog, "IDBT2ICacheLoad: Can't open [%s].\n", name );

exit:
    if( f )
        fclose( f );

    return (C.z >= 0);
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheClear ---------------------------------------------- */
/* --------------------------------------------------------------- */

void IDBT2ICacheClear( T2ICache &C )
{
    C.m.clear();
    C.z = -1;
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheClear ---------------------------------------------- */
/* --------------------------------------------------------------- */

void IDBT2ICacheClear()
{
    IDBT2ICacheClear( _C1 );
    IDBT2ICacheClear( _C2 );
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheFetch ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Get pointer into cache.
//
static bool IDBT2ICacheFetch(
    const Til2Img*&	t2i,
    T2ICache		&C,
    int				z,
    int				id,
    FILE			*flog )
{
    map<int,Til2Img>::iterator	it = C.m.find( id );

    if( it == C.m.end() ) {

        fprintf( flog,
        "IDBT2ICacheFetch: No entry for [%d %d].\n", z, id );

        t2i = NULL;
        return false;
    }
    else
        t2i = &it->second;

    return true;
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheNGet1 ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Cache data for this z and get pointer to this tile.
//
// Call IDBT2ICacheClear() to explicitly recover
// cache memory, but see notes there.
//
bool IDBT2ICacheNGet1(
    const Til2Img*&	t2i,
    const string	&idb,
    int				z,
    int				id,
    FILE			*flog )
{
// load cache

    if( z != _C1.z && !IDBT2ICacheLoad( _C1, idb, z, flog ) )
        return false;

// read from cache

    return IDBT2ICacheFetch( t2i, _C1, z, id, flog );
}

/* --------------------------------------------------------------- */
/* IDBT2ICacheNGet2 ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Cache data for up to two z's and get pointers to tile data.
//
// Call IDBT2ICacheClear() to explicitly recover
// cache memory, but see notes there.
//
bool IDBT2ICacheNGet2(
    const Til2Img*&	t2i1,
    const Til2Img*&	t2i2,
    const string	&idb,
    int				z1,
    int				id1,
    int				z2,
    int				id2,
    FILE			*flog )
{
// assign/load caches

    T2ICache	*p1 = &_C1, *p2;

    if( z1 == z2 ) {

        if( z1 == _C2.z )
            p1 = &_C2;

        p2 = p1;

        if( z1 != p1->z && !IDBT2ICacheLoad( *p1, idb, z1, flog ) )
            return false;
    }
    else {

        p2 = &_C2;

        if( z1 == _C2.z || z2 == _C1.z ) {
            p1 = &_C2;
            p2 = &_C1;
        }

        if( z1 != p1->z && !IDBT2ICacheLoad( *p1, idb, z1, flog ) )
            return false;

        if( z2 != p2->z && !IDBT2ICacheLoad( *p2, idb, z2, flog ) )
            return false;
    }

    return	IDBT2ICacheFetch( t2i1, *p1, z1, id1, flog ) &&
            IDBT2ICacheFetch( t2i2, *p2, z2, id2, flog );
}

/* --------------------------------------------------------------- */
/* IDBTil2FM ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given IDBPATH/TileToFM file for this tile's data.
//
bool IDBTil2FM(
    Til2FM			&t2f,
    const string	&idb,
    int				z,
    int				id,
    FILE			*flog )
{
    char	name[2048];

// old-style z/id dir hierarchy

    if( idb.empty() ) {

        int	len;

        // try name as tif
        len = sprintf( name, "../%d/%d/fm.tif", z, id );

        if( !DskExists( name ) ) {
            // assume png
            name[len-3] = 'p';
            name[len-2] = 'n';
            name[len-1] = 'g';
        }

        t2f.path	= name;
        t2f.id		= id;
        return true;
    }

// new way using idb

    FILE	*f;
    int		ok = false;

    sprintf( name, "%s/%d/TileToFM.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 ) {
            fprintf( flog, "IDBTil2FM: Empty file [%s].\n", name );
            goto exit;
        }

        while( LS.Get( f ) > 0 ) {

            char	buf[2048];

            t2f.id = -1;

            sscanf( LS.line, "%d", &t2f.id );

            if( t2f.id != id )
                continue;

            sscanf( LS.line, "%d\t%[^\t\n]", &t2f.id, buf );
            t2f.path	= buf;
            ok			= true;
            goto exit;
        }

        fprintf( flog, "IDBTil2FM: No entry for [%d %d].\n", z, id );
    }
    else
        fprintf( flog, "IDBTil2FM: Can't open [%s].\n", name );

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* IDBTil2FMD ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given IDBPATH/TileToFMd file for this tile's data.
//
// Unlike IDBTil2FM, this function does not print error messages
// upon failure; it just returns false. In that case, the caller
// should attempt to get the standard fm.
//
bool IDBTil2FMD(
    Til2FM			&t2f,
    const string	&idb,
    int				z,
    int				id )
{
    char	name[2048];

// old-style z/id dir hierarchy

    if( idb.empty() ) {

        int	len;

        // try name as tif
        len = sprintf( name, "../%d/%d/fmd.tif", z, id );

        if( !DskExists( name ) ) {
            // assume png
            name[len-3] = 'p';
            name[len-2] = 'n';
            name[len-1] = 'g';
        }

        t2f.path	= name;
        t2f.id		= id;
        return true;
    }

// new way using idb

    FILE	*f;
    int		ok = false;

    sprintf( name, "%s/%d/TileToFMD.txt", idb.c_str(), z );

    if( f = fopen( name, "r" ) ) {

        CLineScan	LS;

        if( LS.Get( f ) <= 0 )
            goto exit;

        while( LS.Get( f ) > 0 ) {

            char	buf[2048];

            t2f.id = -1;

            sscanf( LS.line, "%d", &t2f.id );

            if( t2f.id != id )
                continue;

            sscanf( LS.line, "%d\t%[^\t\n]", &t2f.id, buf );
            t2f.path	= buf;
            ok			= true;
            goto exit;
        }
    }

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* PrintTil2Img -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print Til2Img entry; cAB = {'A', 'B'}.
//
void PrintTil2Img( FILE *flog, int cAB, const Til2Img &t2i )
{
    char	buf[128];

    if( t2i.col != -999 ) {
        sprintf( buf, " c.r.cam=[%d.%d.%d]",
            t2i.col, t2i.row, t2i.cam );
    }
    else
        buf[0] = 0;

    fprintf( flog, "Til2Img entry: %c"
    " T=[%7.4f %7.4f %8.2f %7.4f %7.4f %8.2f]"
    "%s"
    " path=[%s].\n",
    cAB,
    t2i.T.t[0], t2i.T.t[1], t2i.T.t[2],
    t2i.T.t[3], t2i.T.t[4], t2i.T.t[5],
    buf,
    t2i.path.c_str() );
}

/* --------------------------------------------------------------- */
/* PrintTil2FM --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print Til2FM entry; cAB = {'A', 'B'}.
//
void PrintTil2FM( FILE *flog, int cAB, const Til2FM &t2f )
{
    fprintf( flog, "Til2FM entry: %c path=[%s].\n",
    cAB, t2f.path.c_str() );
}

/* --------------------------------------------------------------- */
/* ReadThmPair --------------------------------------------------- */
/* --------------------------------------------------------------- */

bool ReadThmPair(
    ThmPair	&tpr,
    int		alr,
    int		atl,
    int		acr,
    int		blr,
    int		btl,
    int		bcr,
    FILE	*flog )
{
    CMutex	M;
    char	name[256];
    FILE	*f;
    int		ok = false;

    sprintf( name, "tpr_%d_%d", alr, blr );

    if( M.Get( name ) ) {

        sprintf( name, "ThmPair_%d^%d.txt", alr, blr );
        f = fopen( name, "r" );

        if( f ) {

            CLineScan	LS;

            if( LS.Get( f ) <= 0 ) {

                fprintf( flog,
                "ReadThmPair: Empty file [%s].\n", name );
                goto exit;
            }

            while( LS.Get( f ) > 0 ) {

                sscanf( LS.line,
                "%d\t%d\t%d\t%d\t%d"
                "\t%lf\t%lf"
                "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                &tpr.atl, &tpr.acr, &tpr.btl, &tpr.bcr, &tpr.err,
                &tpr.A, &tpr.R,
                &tpr.T.t[0], &tpr.T.t[1], &tpr.T.t[2],
                &tpr.T.t[3], &tpr.T.t[4], &tpr.T.t[5] );

                if( tpr.atl != atl || tpr.btl != btl ||
                    tpr.acr != acr || tpr.bcr != bcr ) {

                    continue;
                }

                fprintf( flog, "ReadThmPair: Got entry: "
                "A=%f, R=%f, T=[%f %f %f %f %f %f].\n",
                tpr.A, tpr.R,
                tpr.T.t[0], tpr.T.t[1], tpr.T.t[2],
                tpr.T.t[3], tpr.T.t[4], tpr.T.t[5] );

                ok = true;
                goto exit;
            }

            fprintf( flog,
            "ReadThmPair: No entry for %d.%d-%d^%d.%d-%d\n",
            alr, atl, acr, blr, btl, bcr );
        }
        else
            fprintf( flog, "ReadThmPair: Can't open [%s].\n", name );

exit:
        if( f )
            fclose( f );
    }

    M.Release();

    return ok;
}

/* --------------------------------------------------------------- */
/* ReadAllThmPair ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return true if successfully read file. Even so, entry count
// may be zero, so always check tpr.size().
//
bool ReadAllThmPair(
    vector<ThmPair>	&tpr,
    int				alr,
    int				blr,
    FILE			*flog )
{
    CMutex	M;
    char	name[256];
    FILE	*f;
    int		ok = false;

    tpr.clear();

    sprintf( name, "tpr_%d_%d", alr, blr );

    if( M.Get( name ) ) {

        sprintf( name, "ThmPair_%d^%d.txt", alr, blr );
        f = fopen( name, "r" );

        if( f ) {

            CLineScan	LS;

            if( LS.Get( f ) <= 0 ) {

                fprintf( flog,
                "ReadThmPair: Empty file [%s].\n", name );
                goto exit;
            }

            while( LS.Get( f ) > 0 ) {

                ThmPair		P;

                sscanf( LS.line,
                "%d\t%d\t%d\t%d\t%d"
                "\t%lf\t%lf"
                "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                &P.atl, &P.acr, &P.btl, &P.bcr, &P.err,
                &P.A, &P.R,
                &P.T.t[0], &P.T.t[1], &P.T.t[2],
                &P.T.t[3], &P.T.t[4], &P.T.t[5] );

                tpr.push_back( P );
            }

            ok = true;
        }
        else
            fprintf( flog, "ReadThmPair: Can't open [%s].\n", name );

exit:
        if( f )
            fclose( f );
    }

    M.Release();

    return ok;
}

/* --------------------------------------------------------------- */
/* WriteThmPairHdr ----------------------------------------------- */
/* --------------------------------------------------------------- */

void WriteThmPairHdr( FILE *f )
{
    fprintf( f,
    "Atl\tAcr\tBtl\tBcr\tErr\tDeg\tR"
    "\tT0\tT1\tX\tT3\tT4\tY\n" );
}

/* --------------------------------------------------------------- */
/* CreateJobsDir ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CreateJobsDir(
    const char	*lyrdir,
    int			ix,
    int			iy,
    int			za,
    int			zb,
    FILE		*flog )
{
    char	name[2048];
    int		len;

// Create dir
    len = sprintf( name, "%s/%c%d_%d",
            lyrdir, (za == zb ? 'S' : 'D'), ix, iy );
    DskCreateDir( name, flog );

// Create ThmPair file
    if( zb >= 0 ) {
        sprintf( name + len, "/ThmPair_%d^%d.txt", za, zb );
        FILE	*f = FileOpenOrDie( name, "w", flog );
        WriteThmPairHdr( f );
        fclose( f );
    }
}

/* --------------------------------------------------------------- */
/* WriteThmPair -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Note: Table writing bypassed if invalid tile ID.
//
void WriteThmPair(
    const ThmPair	&tpr,
    int				alr,
    int				atl,
    int				acr,
    int				blr,
    int				btl,
    int				bcr )
{
    if( atl < 0 || btl < 0 )
        return;

    CMutex	M;
    char	name[256];

    sprintf( name, "tpr_%d_%d", alr, blr );

    if( M.Get( name ) ) {

        sprintf( name, "ThmPair_%d^%d.txt", alr, blr );
        FILE *f = fopen( name, "a" );

        if( f ) {
            fprintf( f,
                "%d\t%d\t%d\t%d\t%d"
                "\t%f\t%f"
                "\t%f\t%f\t%f\t%f\t%f\t%f\n",
                atl, acr, btl, bcr, tpr.err,
                tpr.A, tpr.R,
                tpr.T.t[0], tpr.T.t[1], tpr.T.t[2],
                tpr.T.t[3], tpr.T.t[4], tpr.T.t[5] );
            fflush( f );
            fclose( f );
        }
    }

    M.Release();
}

/* --------------------------------------------------------------- */
/* ZIDFromFMPath ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given path to a foldmask (or other) file located in the
// Lou-style working directory hierarchy:
//
//	e.g. '/groups/.../temp/z/id/fm.tif'
//	e.g. '../z/id/fm.png'
//
// attempt to extract the z and id values.
//
// Return true if success.
//
bool ZIDFromFMPath( int &z, int &id, const char *path )
{
    const char	*s;
    int			nprev;

// back over filename

    s = path + (nprev = strlen( path ));

    while( nprev-- > 0 && *--s != '/' )
        ;

    if( *s != '/' )
        return false;

// back over id and translate

    while( nprev-- > 0 && *--s != '/' )
        ;

    if( nprev <= 0 )
        return false;

    if( 1 != sscanf( s + 1, "%d", &id ) )
        return false;

// back over z and translate

    while( nprev-- > 0 && *--s != '/' )
        ;

    if( *s != '/' )
        return false;

    if( 1 != sscanf( s + 1, "%d", &z ) )
        return false;

    return true;
}

/* --------------------------------------------------------------- */
/* LoadTAffineTbl_AllZ ------------------------------------------- */
/* --------------------------------------------------------------- */

// From an LSQ-style TAffineTable file, fill out a map (z,id,rgn)
// of TAffine and a set of all unique z-values.
//
void LoadTAffineTbl_AllZ(
    map<MZIDR,TAffine>	&Tmap,
    set<int>			&Zset,
    const char			*path,
    FILE				*flog )
{
    FILE		*f	= FileOpenOrDie( path, "r", flog );
    CLineScan	LS;

    for(;;) {

        if( LS.Get( f ) <= 0 )
            break;

        MZIDR	zir;
        TAffine	T;

        sscanf( LS.line, "%d\t%d\t%d"
        "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &zir.z, &zir.id, &zir.rgn,
        &T.t[0], &T.t[1], &T.t[2],
        &T.t[3], &T.t[4], &T.t[5] );

        Zset.insert( zir.z );
        Tmap[zir] = T;
    }

    fclose( f );
}

/* --------------------------------------------------------------- */
/* LoadTAffineTbl_RngZ ------------------------------------------- */
/* --------------------------------------------------------------- */

// Load entries in range [zi,zf] inclusive.
//
void LoadTAffineTbl_RngZ(
    map<MZIDR,TAffine>	&Tmap,
    int					zi,
    int					zf,
    const char			*path,
    FILE				*flog )
{
    FILE		*f	= FileOpenOrDie( path, "r", flog );
    CLineScan	LS;

    for(;;) {

        if( LS.Get( f ) <= 0 )
            break;

        MZIDR	zir;
        TAffine	T;

        zir.z = atoi( LS.line );

        if( zir.z < zi )
            continue;

        if( zir.z > zf )
            break;

        sscanf( LS.line, "%d\t%d\t%d"
        "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &zir.z, &zir.id, &zir.rgn,
        &T.t[0], &T.t[1], &T.t[2],
        &T.t[3], &T.t[4], &T.t[5] );

        Tmap[zir] = T;
    }

    fclose( f );
}

/* --------------------------------------------------------------- */
/* LoadTHmgphyTbl_AllZ ------------------------------------------- */
/* --------------------------------------------------------------- */

// From an LSQ-style THmgphyTable file, fill out a map (z,id,rgn)
// of TAffine and a set of all unique z-values.
//
void LoadTHmgphyTbl_AllZ(
    map<MZIDR,THmgphy>	&Tmap,
    set<int>			&Zset,
    const char			*path,
    FILE				*flog )
{
    FILE		*f	= FileOpenOrDie( path, "r", flog );
    CLineScan	LS;

    for(;;) {

        if( LS.Get( f ) <= 0 )
            break;

        MZIDR	zir;
        THmgphy	T;

        sscanf( LS.line, "%d\t%d\t%d"
        "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &zir.z, &zir.id, &zir.rgn,
        &T.t[0], &T.t[1], &T.t[2], &T.t[3],
        &T.t[4], &T.t[5], &T.t[6], &T.t[7] );

        Zset.insert( zir.z );
        Tmap[zir] = T;
    }

    fclose( f );
}

/* --------------------------------------------------------------- */
/* LoadTHmgphyTbl_RngZ ------------------------------------------- */
/* --------------------------------------------------------------- */

// Load entries in range [zi,zf] inclusive.
//
void LoadTHmgphyTbl_RngZ(
    map<MZIDR,THmgphy>	&Tmap,
    int					zi,
    int					zf,
    const char			*path,
    FILE				*flog )
{
    FILE		*f	= FileOpenOrDie( path, "r", flog );
    CLineScan	LS;

    for(;;) {

        if( LS.Get( f ) <= 0 )
            break;

        MZIDR	zir;
        THmgphy	T;

        zir.z = atoi( LS.line );

        if( zir.z < zi )
            continue;

        if( zir.z > zf )
            break;

        sscanf( LS.line, "%d\t%d\t%d"
        "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &zir.z, &zir.id, &zir.rgn,
        &T.t[0], &T.t[1], &T.t[2], &T.t[3],
        &T.t[4], &T.t[5], &T.t[6], &T.t[7] );

        Tmap[zir] = T;
    }

    fclose( f );
}


