

#pragma once


#include	"GenDefs.h"
#include	"CTForm.h"
#include	"CRegexID.h"

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CUTile {

public:
	string	name;
	int		z;		// z layer
	int		id;		// inlayer id
	int		ix;		// idx to aux data
	TForm	T;		// local to global

public:
	CUTile()	{ix = 0;};
};


typedef struct TSAux {
	double	r;
	TForm	inv;
} TSAux;

/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CTileSet {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID		re_id;
	FILE*			flog;
	int				gW, gH;	// each tile dims

public:
	vector<CUTile>	vtil;
	vector<TSAux>	vaux;

public:
	void SetLogFile( FILE* fout )
		{flog = fout;};

	void SetDecoderPat( const char *pat );
	int	 DecodeID( const char *name );

	void FillFromRickFile( const char *path, int zmin, int zmax );
	void FillFromTrakEM2( const char *path, int zmin, int zmax );
	void FillFromIDB( const string &idb, int zmin, int zmax );

	void SetTileDims( int w, int h )	{gW = w; gH = h;};
	void SetTileDimsFromImageFile();
	void SetTileDimsFromIDB( const string &idb );

	void GetTileDims( int &w, int &h )
		{w = gW; h = gH;};

	void InitAuxData();

	void SortAll_z();
	void SortAll_z_id();
	void SortAll_z_r();

	void GetOrder_id( vector<int> &order, int is0, int isN );

	void GetLayerLimits( int &i0, int &iN );

	void LayerBounds( DBox &B, int is0, int isN );
	void LayerAssignR( int is0, int isN, const DBox &B );

	double ABOlap( int a, int b );
};


