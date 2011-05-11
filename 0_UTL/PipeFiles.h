

#pragma once


#include	"CTForm.h"


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
				HFANG_SL,
				HFANG_XL,
				QTRSH_SL,
				QTRSH_XL,
				RTRSH_SL,
				RTRSH_XL;
	int			FLD,
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
	char	path[2048];
} Til2Img;

typedef struct {
// entry: ThmPair.txt
	TForm	T;
	double	A, Q, R;
	int		atl, btl,
			acr, bcr,
			err;
} ThmPair;

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void OpenPairLog( int alr, int atl, int blr, int btl );

bool ReadThmParams( ThmParams &T, FILE *flog );
bool ReadMeshParams( MeshParams &M, FILE *flog );

bool ReadTil2Img(
	Til2Img	&t2i,
	int		layer,
	int		tile,
	FILE	*flog );

void PrintTil2Img( FILE *flog, int cAB, const Til2Img &t2i );

bool ReadThmPair(
	ThmPair	&tpr,
	int		alr,
	int		atl,
	int		acr,
	int		blr,
	int		btl,
	int		bcr,
	FILE	*flog );

bool ReadAllThmPair(
	vector<ThmPair>	&tpr,
	int				alr,
	int				blr,
	FILE			*flog );

void WriteThmPair(
	const ThmPair	&tpr,
	int				alr,
	int				atl,
	int				acr,
	int				blr,
	int				btl,
	int				bcr );


