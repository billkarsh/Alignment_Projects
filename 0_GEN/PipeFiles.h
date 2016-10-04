

#pragma once


#include	"GenDefs.h"
#include	"TAffine.h"
#include	"THmgphy.h"

#include	<map>
#include	<set>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
    string	croprectfile;
    double	mintileolapfrac,
            stripsweepspan,
            stripsweepstep,
            stripmincorr,
            blocksweepspan,
            blocksweepstep,
            blockxyconf,
            blockmincorr,
            blocknomcorr,
            blocknomcoverage;
    int		slotspernode,
            usingfoldmasks,
            makefmjparam,
            makefmslots,
            createauxdirs,
            montageblocksize,
            ignorecorners,
            makesamejparam,
            makesameslots,
            crossscale,
            legendremaxorder,
            rendersdevcnts,
            maskoutresin,
            stripwidth,
            stripslots,
            crossblocksize,
            blockreqdz,
            blockmaxdz,
            blockslots,
            makedownjparam,
            makedownslots,
            xmlpixeltype,
            xmlsclmin,
            xmlsclmax;
} ScriptParams;

typedef struct {
    double	SCALE,
            XSCALE,
            YSCALE,
            SKEW,
            XYCONF_SL,
            XYCONF_XL,
            NBMXHT_SL,
            NBMXHT_XL,
            HFANGDN_SL,
            HFANGDN_XL,
            HFANGPR_SL,
            HFANGPR_XL,
            RTRSH_SL,
            RTRSH_XL,
            RIT_SL,
            RIT_XL,
            RFA_SL,
            RFA_XL,
            RFT_SL,
            RFT_XL,
            TMC,
            TSC,
            IFM,
            FFM,
            FYL,
            CPD,
            EMT,
            LDA,
            LDR,
            LDC,
            DXY;
    int		PXBRO,
            PXLENS,
            PXRESMSK,
            PXDOG,
            PXDOG_R1,
            PXDOG_R2,
            FLD,
            PRETWEAK,
            MODE_SL,
            MODE_XL,
            TAB2DFM_SL,
            TAB2DFM_XL,
            THMDEC_SL,
            THMDEC_XL,
            OLAP1D_SL,
            OLAP1D_XL,
            OLAP2D_SL,
            OLAP2D_XL,
            TWEAKS,
            LIMXY_SL,
            LIMXY_XL,
            WTHMPR,
            OPT_SL,
            MNL,
            MTA,
            MMA,
            ONE,
            EMM,
            WDI,
            WMT,
            WTT;
} MatchParams;

typedef struct {
// entry: TileToImage.txt
    string	path;
    TAffine	T;
    int		id,
            col,
            row,
            cam;
} Til2Img;

typedef struct {
// entry: TileToFM(D).txt
    string	path;
    int		id;
} Til2FM;

typedef struct {
    Til2Img		t2i;
    int			z,
                id;
} PicSpec;

typedef struct {
// entry: ThmPair.txt
    TAffine	T;
    double	A, R;
    int		atl, acr,
            btl, bcr,
            err;
} ThmPair;

class T2ICache {
// all t2i for this z
public:
    map<int,Til2Img>	m;
    int					z;
public:
    T2ICache() : z(-1) {};
};

/* --------------------------------------------------------------- */
/* Read tforms organized as Z-files ------------------------------ */
/* --------------------------------------------------------------- */

namespace ns_pipergns {

enum RgnFlags {
    fbDead		= 0x01,
    fbRead		= 0x02,
    fbPnts		= 0x04,
    fbKill		= 0x08,
    fbCutd		= 0x10,

    fmRead		= fbRead + fbDead
};

#define	FLAG_ISUSED( f )	(((f) & fbDead) == 0)
#define	FLAG_SETUSED( f )	(f = 0)

class Rgns {
// The rgns for given layer
// indexed by 0-based 'idx0'
private:
    FILE			*flog;
    const string	*idb;
public:
    vector<double>	x;		// tform elems
    vector<uint8>	flag;	// rgn flags
    map<int,int>	m;		// map id -> idx0
    int				nr,		// num rgns
                    z,		// common z
                    NE;		// num elems/tform
private:
    bool AFromIDB();
    bool AFromTxt( const char *path );
    bool HFromTxt( const char *path );
    bool ReadXBin( const char *path );
    bool ReadFBin( const char *path );
public:
    bool Init( const string &idb, int iz, FILE *flog );
    bool Load( const char *path );
    // Caller must create dst path BEFORE calling SaveXXX()
    bool SaveBIN( const char *path, bool writeflags );
    bool SaveTXT( const char *path );
};

};	// ns_pipergns

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void OpenPairLog( int alr, int atl, int blr, int btl );

char *NamePtsFile( char *buf, int alr, int blr );
char *NameLogFile( char *buf, int alr, int atl, int blr, int btl );

bool ReadScriptParams(
    ScriptParams	&S,
    const char		*scriptpath,
    FILE			*flog = stdout );

bool ReadMatchParams(
    MatchParams		&M,
    int				alr,
    int				blr,
    const char		*matchparamspath = NULL,
    FILE			*flog = stdout );

void IDBFromTemp(
    string		&idbpath,
    const char	*tempdir,
    FILE		*flog = stdout );

bool IDBGetImageDims(
    int				&w,
    int				&h,
    const string	&idb,
    FILE			*flog = stdout );

int IDBGetIDRgnMap(
    map<int,int>	&m,
    const string	&idb,
    int				z,
    FILE			*flog = stdout );

bool IDBT2IGet1(
    Til2Img			&t2i,
    const string	&idb,
    int				z,
    int				id,
    const char		*forcepath = NULL,
    FILE			*flog = stdout );

bool IDBT2IGetAll(
    vector<Til2Img>	&t2i,
    const string	&idb,
    int				z,
    FILE			*flog = stdout );

bool IDBT2IGet_JustIDandT(
    vector<Til2Img>	&t2i,
    const string	&idb,
    int				z,
    FILE			*flog = stdout );

bool IDBT2ICacheLoad(
    T2ICache		&C,
    const string	&idb,
    int				z,
    FILE			*flog = stdout );

void IDBT2ICacheClear( T2ICache &C );
void IDBT2ICacheClear();

bool IDBT2ICacheNGet1(
    const Til2Img*&	t2i,
    const string	&idb,
    int				z,
    int				id,
    FILE			*flog = stdout );

bool IDBT2ICacheNGet2(
    const Til2Img*&	t2i1,
    const Til2Img*&	t2i2,
    const string	&idb,
    int				z1,
    int				id1,
    int				z2,
    int				id2,
    FILE			*flog = stdout );

bool IDBTil2FM(
    Til2FM			&t2f,
    const string	&idb,
    int				z,
    int				id,
    FILE			*flog = stdout );

bool IDBTil2FMD(
    Til2FM			&t2f,
    const string	&idb,
    int				z,
    int				id );

void PrintTil2Img( FILE *flog, int cAB, const Til2Img &t2i );
void PrintTil2FM( FILE *flog, int cAB, const Til2FM &t2f );

bool ReadThmPair(
    ThmPair	&tpr,
    int		alr,
    int		atl,
    int		acr,
    int		blr,
    int		btl,
    int		bcr,
    FILE	*flog = stdout );

bool ReadAllThmPair(
    vector<ThmPair>	&tpr,
    int				alr,
    int				blr,
    FILE			*flog = stdout );

void WriteThmPairHdr( FILE *f );

void CreateJobsDir(
    const char	*lyrdir,
    int			ix,
    int			iy,
    int			za,
    int			zb,
    FILE		*flog );

void WriteThmPair(
    const ThmPair	&tpr,
    int				alr,
    int				atl,
    int				acr,
    int				blr,
    int				btl,
    int				bcr );

bool ZIDFromFMPath( int &z, int &id, const char *path );

void LoadTAffineTbl_AllZ(
    map<MZIDR,TAffine>	&Tmap,
    set<int>			&Zset,
    const char			*path,
    FILE				*flog = NULL );

void LoadTAffineTbl_RngZ(
    map<MZIDR,TAffine>	&Tmap,
    int					zi,
    int					zf,
    const char			*path,
    FILE				*flog = NULL );

void LoadTHmgphyTbl_AllZ(
    map<MZIDR,THmgphy>	&Tmap,
    set<int>			&Zset,
    const char			*path,
    FILE				*flog = NULL );

void LoadTHmgphyTbl_RngZ(
    map<MZIDR,THmgphy>	&Tmap,
    int					zi,
    int					zf,
    const char			*path,
    FILE				*flog = NULL );


