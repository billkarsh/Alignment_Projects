

#pragma once


#include	"GenDefs.h"
#include	"TAffine.h"

#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum TSConst {
// ApplyClix()::tfType
	tsClixSimilar	= 0,
	tsClixAffine	= 1
};

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CUTile {

public:
	string	name;
	TAffine	T;		// local to global
	int		z,		// z layer
			id,		// inlayer id
			col,	// grid col index
			row,	// grid row index
			cam,	// camera index
			ix;		// idx to aux data

public:
	CUTile() : col(-999), row(-999), cam(0), ix(0) {};
};


typedef struct TSAux {
	double	r;
} TSAux;


typedef struct TSClix {
	int				Az, Bz;
	vector<Point>	A,  B;

	bool operator < ( const TSClix &rhs ) const
		{return Az < rhs.Az;};
} TSClix;

/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CTileSet {

	friend void* _Scape_Paint( void *ithr );

private:
	FILE*			flog;
	int				gW, gH;	// each tile dims

public:
	vector<CUTile>	vtil;
	vector<TSAux>	vaux;

public:
	CTileSet() : flog(stdout), gW(0), gH(0) {};

	void SetLogFile( FILE* fout )
		{flog = fout;};

	void FillFromRickFile( const char *path, int zmin, int zmax );
	void FillFromTrakEM2( const char *path, int zmin, int zmax );
	void FillFromIDB( const string &idb, int zmin, int zmax );
	void FillFromRgns(
		const char		*path,
		const string	&idb,
		int				zmin,
		int				zmax );

	void SetTileDims( int w, int h )	{gW = w; gH = h;};
	void SetTileDimsFromImageFile();
	void SetTileDimsFromIDB( const string &idb );

	void GetTileDims( int &w, int &h )
		{w = gW; h = gH;};

	void InitAuxData();

	void SortAll_z();
	void SortAll_z_id();
	void SortAll_z_r();

	int  GetOrder_id( vector<int> &order, int is0, int isN );

	void GetLayerLimits( int &i0, int &iN );
	void GetLayerLimitsR( int &i0, int &iN );

	void ReadClixFile( vector<TSClix> &clk, const char *path );
	TAffine SimilarityFromClix( const TSClix &clk );
	TAffine AffineFromClix( const TSClix &clk );
	void ApplyClix( int tfType, const char *path );

	void BoundsPlus1( DBox &B, const vector<Point> &cnr, int i );
	void LayerBounds( DBox &B, int is0, int isN );
	void AllBounds( DBox &B );
	void Reframe( DBox &B );

	void LayerAssignR( int is0, int isN, const DBox &B );

	double ABOlap( int a, int b, const TAffine *Tab = NULL );

	void WriteTrakEM2Layer(
		FILE*	f,
		int		&oid,
		int		xmltype,
		int		xmlmin,
		int		xmlmax,
		int		is0,
		int		isN );

	void WriteTrakEM2(
		const char	*path,
		DBox		&B,
		int			xmltype,
		int			xmlmin,
		int			xmlmax );

	void WriteTrakEM2_EZ(
		const char	*path,
		int			xmltype,
		int			xmlmin,
		int			xmlmax );

	void WriteTileToImage(
		const char	*topdir,
		bool		create_zdir,
		bool		create_nmrcdir,
		int			is0,
		int			isN );

private:
	void Scape_AdjustBounds(
		uint32				&ws,
		uint32				&hs,
		double				&x0,
		double				&y0,
		vector<TAffine>		&vTadj,
		const vector<int>	&vid,
		double				scale,
		int					szmult ) const;

	void Scape_PaintTH( int nthr ) const;

public:
	uint8* Scape(
		uint32				&ws,
		uint32				&hs,
		double				&x0,
		double				&y0,
		const vector<int>	&vid,
		double				scale,
		int					szmult,
		int					bkval,
		int					lgord,
		int					sdnorm,
		bool				resmask,
		int					nthr ) const;
};


