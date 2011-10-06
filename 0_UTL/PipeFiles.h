

#pragma once


#include	"CTForm.h"

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	int			PXBRO,
				PXDOG,
				PXDOG_R1,
				PXDOG_R2;
} PxPairPrm;

typedef struct {
	PxPairPrm	pxp;
	double		SCALE,
				XSCALE,
				YSCALE,
				SKEW,
				HFANGDN_SL,
				HFANGDN_XL,
				HFANGPR_SL,
				HFANGPR_XL,
				QTRSH_SL,
				QTRSH_XL,
				RTRSH_SL,
				RTRSH_XL;
	int			FLD,
				SLOPPY_SL,
				OLAP1D,
				OLAP2D_SL,
				OLAP2D_XL,
				TWEAKS,
				INPALN,
				DINPUT;
} ThmParams;

typedef struct {
	PxPairPrm	pxp;
	double		DIT,
				DAF,
				DFT,
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
	int			FLD,
				DSL,
				MNL,
				MTA,
				MMA,
				ONE,
				EMM,
				WDI,
				WMT,
				WTT;
} MeshParams;

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

bool ReadThmParams( ThmParams &T, int layer, FILE *flog = stdout );
bool ReadMeshParams( MeshParams &M, int layer, FILE *flog = stdout );

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

void WriteThmPair(
	const ThmPair	&tpr,
	int				alr,
	int				atl,
	int				acr,
	int				blr,
	int				btl,
	int				bcr );

bool ZIDFromFMPath( int &z, int &id, const char *path );


