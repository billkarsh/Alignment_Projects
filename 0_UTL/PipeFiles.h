

#pragma once


#include	"GenDefs.h"
#include	"CTForm.h"

#include	<map>
#include	<set>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

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
			DIT_SL,
			DIT_XL,
			DFA_SL,
			DFA_XL,
			DFT_SL,
			DFT_XL,
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
			PXDOG,
			PXDOG_R1,
			PXDOG_R2,
			FLD,
			PRETWEAK,
			MODE_SL,
			MODE_XL,
			THMDEC,
			OLAP1D_SL,
			OLAP1D_XL,
			OLAP2D_SL,
			OLAP2D_XL,
			TWEAKS,
			LIMXY_SL,
			LIMXY_XL,
			DSL,
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
	TForm	T;
	int		tile;
	string	path;
} Til2Img;

typedef struct {
// entry: TileToFM(D).txt
	int		tile;
	string	path;
} Til2FM;

typedef struct {
// entry: ThmPair.txt
	TForm	T;
	double	A, R;
	int		atl, btl,
			acr, bcr,
			err;
} ThmPair;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void OpenPairLog( int alr, int atl, int blr, int btl );

bool ReadMatchParams(
	MatchParams		&M,
	int				alr,
	int				blr,
	FILE			*flog = stdout );

void IDBReadImgParams( string &idbpath, FILE *flog = stdout );

bool IDBAllTil2Img(
	vector<Til2Img>	&t2i,
	const string	&idb,
	int				layer,
	FILE			*flog = stdout );

bool IDBTil2Img(
	Til2Img			&t2i,
	const string	&idb,
	int				layer,
	int				tile,
	const char		*forcepath = NULL,
	FILE			*flog = stdout );

bool IDBTil2FM(
	Til2FM			&t2f,
	const string	&idb,
	int				layer,
	int				tile,
	FILE			*flog = stdout );

bool IDBTil2FMD(
	Til2FM			&t2f,
	const string	&idb,
	int				layer,
	int				tile );

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

void LoadTFormTbl_AllZ(
	map<MZID,TForm>	&Tmap,
	set<int>		&Zset,
	const char		*path,
	FILE			*flog = NULL );

void LoadTFormTbl_ThisZ(
	map<MZIDR,TForm>	&Tmap,
	int					z,
	const char			*path,
	FILE				*flog = NULL );


