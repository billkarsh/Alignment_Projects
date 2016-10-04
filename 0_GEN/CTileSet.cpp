

#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"PipeFiles.h"
#include	"CTileSet.h"
#include	"LinEqu.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Geometry.h"

#include	<algorithm>
using namespace std;
using namespace ns_pipergns;


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* RejectTile ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Chance to optionally apply rejection criteria against tile.
//
// Return true to reject.
//
static bool RejectTile( const CUTile &til )
{
// ------------------------------------
// accept only col/row subsection
#if 0
    if( til.col != -999 ) {
        int	col = til.col, row = til.row;
//		if( col < 56 || col > 68 || row < 55 || row > 67 )
        if( col < 33 || col > 55 || row < 73 || row > 103 )
            return true;
    }
#endif
// ------------------------------------

    return false;
}

/* --------------------------------------------------------------- */
/* FillFromRickFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::FillFromRickFile( const char *path, int zmin, int zmax )
{
    FILE	*fp = FileOpenOrDie( path, "r", flog );

/* ---------- */
/* Scan lines */
/* ---------- */

    for( ;; ) {

        CUTile	til;
        char	name[2048];
        int		z;

        /* ---------- */
        /* Get a line */
        /* ---------- */

        if( fscanf( fp,
            "%d%d"
            "%lf%lf%lf%lf%lf%lf"
            "%d%d%d%s\n",
            &z, &til.id,
            &til.T.t[0], &til.T.t[1], &til.T.t[2],
            &til.T.t[3], &til.T.t[4], &til.T.t[5],
            &til.col, &til.row, &til.cam, name ) != 12 ) {

            break;
        }

        if( z > zmax )
            continue;

        if( z < zmin )
            continue;

        /* --------- */
        /* Set entry */
        /* --------- */

        til.name	= name;
        til.z		= z;

        if( !RejectTile( til ) )
            vtil.push_back( til );
    }

/* ----- */
/* Close */
/* ----- */

    fclose( fp );
}

/* --------------------------------------------------------------- */
/* XMLGetTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static TiXmlElement* XMLGetTiles(
    CTileSet		*TS,
    TiXmlElement*	layer,
    int				z )
{
    TiXmlElement*	pfirst	= layer->FirstChildElement( "t2_patch" );
    TiXmlElement*	ptch	= pfirst;

    for( ; ptch; ptch = ptch->NextSiblingElement() ) {

        CUTile		til;
        const char	*name = ptch->Attribute( "file_path" );

        til.name	= name;
        til.z		= z;
        til.id		= IDFromPatch( ptch );
        til.T.ScanTrackEM2( ptch->Attribute( "transform" ) );

        // extract col, row, cam if present
        {
            int	rgn, col, row, cam;

            if( 6 == sscanf( ptch->Attribute( "title" ),
                "%d.%d-%d_%d.%d.%d",
                &z, &til.id, &rgn,
                &col, &row, &cam ) ) {

                til.col = col;
                til.row = row;
                til.cam = cam;
            }
        }

        if( !RejectTile( til ) )
            TS->vtil.push_back( til );
    }

    return pfirst;
}

/* --------------------------------------------------------------- */
/* FillFromTrakEM2 ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Automatically sets (gW,gH) from xml file.
//
void CTileSet::FillFromTrakEM2( const char *path, int zmin, int zmax )
{
/* ---- */
/* Open */
/* ---- */

    XML_TKEM		xml( path, flog );
    TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* For each layer */
/* -------------- */

    for( ; layer; layer = layer->NextSiblingElement() ) {

        int	z = atoi( layer->Attribute( "z" ) );

        if( z > zmax )
            break;

        if( z < zmin )
            continue;

        TiXmlElement* p0 = XMLGetTiles( this, layer, z );

        if( p0 && !gW ) {
            gW = atoi( p0->Attribute( "width" ) );
            gH = atoi( p0->Attribute( "height" ) );
        }
    }
}

/* --------------------------------------------------------------- */
/* FillFromIDB --------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::FillFromIDB( const string &idb, int zmin, int zmax )
{
    for( int z = zmin; z <= zmax; ++z ) {

        vector<Til2Img>	t2i;

        if( IDBT2IGetAll( t2i, idb, z, flog ) ) {

            int	nt = t2i.size();

            for( int i = 0; i < nt; ++i ) {

                const Til2Img&	E = t2i[i];
                CUTile			til;

                til.name	= E.path;
                til.T		= E.T;
                til.z		= z;
                til.id		= E.id;
                til.col		= E.col;
                til.row		= E.row;
                til.cam		= E.cam;

                if( !RejectTile( til ) )
                    vtil.push_back( til );
            }
        }
    }
}

/* --------------------------------------------------------------- */
/* FillFromRgns -------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::FillFromRgns(
    const char		*path,
    const string	&idb,
    int				zmin,
    int				zmax )
{
    if( !path ) {
        FillFromIDB( idb, zmin, zmax );
        return;
    }

    for( int z = zmin; z <= zmax; ++z ) {

        Rgns	R;

        if( !R.Init( idb, z, flog ) )
            continue;

        if( !R.Load( path ) )
            continue;

        map<int,int>::iterator	mi, en = R.m.end();

        for( mi = R.m.begin(); mi != en; ) {

            const Til2Img	*t2i;

            int	id		= mi->first,
                j0		= mi->second,
                jlim	= (++mi == en ? R.nr : mi->second);

            if( !IDBT2ICacheNGet1( t2i, idb, z, id, flog ) )
                continue;

            for( int j = j0; j < jlim; ++j ) {

                if( !FLAG_ISUSED( R.flag[j] ) )
                    continue;

                CUTile	til;

                til.name	= t2i->path;
                til.T		= X_AS_AFF( R.x, j );
                til.z		= z;
                til.id		= id;
                til.col		= t2i->col;
                til.row		= t2i->row;
                til.cam		= t2i->cam;

                if( !RejectTile( til ) )
                    vtil.push_back( til );
            }
        }
    }
}

/* --------------------------------------------------------------- */
/* SetTileDimsFromImageFile -------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::SetTileDimsFromImageFile()
{
    uint32	w, h;
    uint8	*ras = Raster8FromAny( vtil[0].name.c_str(), w, h, flog );

    if( !ras || !w ) {
        fprintf( flog, "Error loading [%s].\n", vtil[0].name.c_str() );
        exit( 42 );
    }

    RasterFree( ras );

    SetTileDims( w, h );
}

/* --------------------------------------------------------------- */
/* SetTileDimsFromIDB -------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::SetTileDimsFromIDB( const string &idb )
{
    if( !IDBGetImageDims( gW, gH, idb, flog ) )
        exit( 42 );
}

/* --------------------------------------------------------------- */
/* InitAuxData --------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::InitAuxData()
{
    int	nt = vtil.size();

    if( nt ) {

        vaux.resize( nt );

        for( int i = 0; i < nt; ++i ) {

            vtil[i].ix	= i;
            vaux[i].r	= 0.0;
        }
    }
}

/* --------------------------------------------------------------- */
/* Sort_z_inc ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_z_inc( const CUTile &A, const CUTile &B )
{
    return A.z < B.z;
}

/* --------------------------------------------------------------- */
/* CSort_id_inc -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CSort_id_Inc {
public:
    const vector<CUTile>	&vtil;
public:
    CSort_id_Inc( const vector<CUTile> &vtil )
        : vtil(vtil) {};
    bool operator() ( int a, int b )
        {return vtil[a].id < vtil[b].id;};
};

/* --------------------------------------------------------------- */
/* Sort_z_id_inc ------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_z_id_inc( const CUTile &A, const CUTile &B )
{
    return A.z < B.z || (A.z == B.z && A.id < B.id);
}

/* --------------------------------------------------------------- */
/* CSort_r_id_inc ------------------------------------------------ */
/* --------------------------------------------------------------- */

class CSort_r_id_inc {
public:
    const vector<TSAux>	&vaux;
public:
    CSort_r_id_inc( const vector<TSAux> &vaux )
        : vaux(vaux) {};
    bool operator() ( const CUTile &A, const CUTile &B )
        {return
            vaux[A.ix].r < vaux[B.ix].r ||
            (vaux[A.ix].r == vaux[B.ix].r && A.id < B.id);};
};

/* --------------------------------------------------------------- */
/* SortAll_z ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::SortAll_z()
{
    sort( vtil.begin(), vtil.end(), Sort_z_inc );
}

/* --------------------------------------------------------------- */
/* SortAll_z_id -------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::SortAll_z_id()
{
    sort( vtil.begin(), vtil.end(), Sort_z_id_inc );
}

/* --------------------------------------------------------------- */
/* SortAll_z_r --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Within each layer/montage, sort the tiles into ascending
// order of distance from center.
//
void CTileSet::SortAll_z_r()
{
    if( !vaux.size() )
        InitAuxData();

// First, just sort into layers

    SortAll_z();

// For each layer...

    int	is0, isN;

    GetLayerLimits( is0 = 0, isN );

    while( isN != -1 ) {

        // Get montage bounds
        // Assign radii from montage center
        // Sort this layer by r

        DBox			B;
        CSort_r_id_inc	sorter( vaux );

        LayerBounds( B, is0, isN );
        LayerAssignR( is0, isN, B );
        sort( vtil.begin() + is0, vtil.begin() + isN, sorter );

        GetLayerLimits( is0 = isN, isN );
    }
}

/* --------------------------------------------------------------- */
/* GetOrder_id --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Size and fill caller's order vector with indices of tiles
// in layer [is0,isN) and in order of tile id.
//
// Return updated loop limit (isN - is0).
//
int CTileSet::GetOrder_id( vector<int> &order, int is0, int isN )
{
    isN -= is0;

    order.resize( isN );

    for( int i = 0; i < isN; ++i )
        order[i] = is0 + i;

    CSort_id_Inc	sorter( vtil );

    sort( order.begin(), order.end(), sorter );

    return isN;
}

/* --------------------------------------------------------------- */
/* GetLayerLimits ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
// Given starting index i0, which selects a layer z, set iN to be
// one beyond the highest index of a picture having same z. This
// makes loop limits [i0,iN) exclusive.
//
// If i0 or iN are out of bounds, both are set to -1.
//
// Usage illustration:
// -------------------
// int	is0, isN;
//
// TS.GetLayerLimits( is0 = 0, isN );
//
// while( isN != -1 ) {
//	/* do something */
//	TS.GetLayerLimits( is0 = isN, isN );
// }
//
void CTileSet::GetLayerLimits( int &i0, int &iN )
{
    int	nt = vtil.size();

    if( i0 < 0 || i0 >= nt ) {

        i0 = -1;
        iN = -1;
    }
    else {

        int	Z = vtil[i0].z;

        for( iN = i0 + 1; iN < nt && vtil[iN].z == Z; ++iN )
            ;
    }
}

/* --------------------------------------------------------------- */
/* GetLayerLimitsR ----------------------------------------------- */
/* --------------------------------------------------------------- */

// This variant is for walking the layers in reverse order.
//
// Tiles must already be sorted by z (subsort doesn't matter).
//
// Given starting index iN, which is one beyond a desired tile
// range, that is, [i0,iN) exclusive, walk backwards and set i0
// accordingly.
//
// If i0 or iN are out of bounds, both are set to -1.
//
// Usage illustration:
// -------------------
// int	is0, isN;
//
// TS.GetLayerLimitsR( is0, isN = TS.vtil.size() );
//
// while( is0 != -1 ) {
//	/* do something */
//	TS.GetLayerLimitsR( is0, isN = is0 );
// }
//
void CTileSet::GetLayerLimitsR( int &i0, int &iN )
{
    if( iN <= 0 ) {

        i0 = -1;
        iN = -1;
    }
    else {

        int	Z = vtil[i0 = iN - 1].z;

        for( ; i0 > 0 && vtil[i0 - 1].z == Z; --i0 )
            ;
    }
}

/* --------------------------------------------------------------- */
/* ReadClixFile -------------------------------------------------- */
/* --------------------------------------------------------------- */

// File format:
// Az Bz Ax1 Ay1 Bx1 By1 Ax2 Ay2 Bx2 By2 [...Axn Ayn Bxn Byn]
//
void CTileSet::ReadClixFile( vector<TSClix> &clk, const char *path )
{
    CLineScan	LS;
    FILE*		f = FileOpenOrDie( path, "r" );

    while( LS.Get( f ) > 0 ) {

        TSClix	C;
        int		N = 0, k;

        if( 2 == sscanf( LS.line, "%d%d%n", &C.Az, &C.Bz, &k ) ) {

            for(;;) {

                Point	A, B;

                N += k + 1;

                if( 4 == sscanf( LS.line + N, "%lf%lf%lf%lf%n",
                            &A.x, &A.y, &B.x, &B.y, &k ) ) {

                    C.A.push_back( A );
                    C.B.push_back( B );
                }
                else
                    break;
            }

            if( C.A.size() >= 2 )
                clk.push_back( C );
        }
    }

    fclose( f );

    sort( clk.begin(), clk.end() );
}

/* --------------------------------------------------------------- */
/* SimilarityFromClix -------------------------------------------- */
/* --------------------------------------------------------------- */

// Return 4 parameter similarity transform (rot + trans): T(A) = B:
//
//		a Ax  -  b Ay  +  c  =  Bx
//		b Ax  +  a Ay  +  d  =  By
//
TAffine CTileSet::SimilarityFromClix( const TSClix &clk )
{
// Create system of normal equations

    double	RHS[4];
    double	LHS[4*4];
    int		np    = clk.A.size(),
            i1[3] = { 0, 1, 2 },
            i2[3] = { 1, 0, 3 };

    Zero_Quick( LHS, RHS, 4 );

    for( int i = 0; i < np; ++i ) {

        const Point&	A = clk.A[i];
        const Point&	B = clk.B[i];

        double	v1[3] = { A.x, -A.y, 1.0 };
        double	v2[3] = { A.x,  A.y, 1.0 };

        AddConstraint_Quick( LHS, RHS, 4, 3, i1, v1, B.x );
        AddConstraint_Quick( LHS, RHS, 4, 3, i2, v2, B.y );
    }

// Solve

    Solve_Quick( LHS, RHS, 4 );

// Return

    return TAffine(
        RHS[0], -RHS[1], RHS[2],
        RHS[1],  RHS[0], RHS[3] );
}

/* --------------------------------------------------------------- */
/* AffineFromClix ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return 6 parameter affine transformation: T(A) = B:
//
//		a Ax  +  b Ay  +  c  =  Bx
//		d Ax  +  e Ay  +  f  =  By
//
TAffine CTileSet::AffineFromClix( const TSClix &clk )
{
    int	np = clk.A.size();

    if( np < 3 )
        return SimilarityFromClix( clk );

// Create system of normal equations

    double	RHS[6];
    double	LHS[6*6];
    int		i1[3] = { 0, 1, 2 },
            i2[3] = { 3, 4, 5 };

    Zero_Quick( LHS, RHS, 6 );

    for( int i = 0; i < np; ++i ) {

        const Point&	A = clk.A[i];
        const Point&	B = clk.B[i];

        double	v[3] = { A.x, A.y, 1.0 };

        AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
        AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
    }

// Solve

    Solve_Quick( LHS, RHS, 6 );

// Return

    return TAffine( &RHS[0] );
}

/* --------------------------------------------------------------- */
/* ApplyClix ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
// Read and apply click pairs file given by path.
//
// The file has following format:
// Az Bz Ax1 Ay1 Bx1 By1 Ax2 Ay2 Bx2 By2 [...Axn Ayn Bxn Byn]
//
// That is, each line specifies a set of points {Ai} on layer
// Az, and a matched set of points {Bi} on BELOW layer Bz;
// again, Bz < Az.
//
void CTileSet::ApplyClix( int tfType, const char *path )
{
    vector<TSClix>	clk;

    ReadClixFile( clk, path );

    int	cmin,
        cmax,
        nc		= clk.size(),
        nt		= vtil.size(),
        zmin	= vtil[0].z,
        zmax	= vtil[nt - 1].z;

// Find lowest applicable click

    for( cmin = 0; cmin < nc - 1; ++cmin ) {
        if( clk[cmin].Bz >= zmin )
            break;
    }

    zmin = clk[cmin].Az;

// Find highest applicable click

    for( cmax = nc - 1; cmax >= cmin; --cmax ) {
        if( clk[cmax].Az <= zmax )
            break;
    }

// Each click affects all layers from clk[i].Az through zmax.
// Assign TAffine vT[k] for each possible layer in range [zmin, zmax].
// The vT[k] accumulate effects of sequential click transforms.

    int				nz = zmax - zmin + 1;
    vector<TAffine>	vT( nz );

    for( int i = cmin; i <= cmax; ++i ) {

        TAffine	R;

        if( tfType == tsClixAffine )
            R = AffineFromClix( clk[i] );
        else
            R = SimilarityFromClix( clk[i] );

        for( int k = clk[i].Az - zmin; k < nz; ++k )
            vT[k] = R * vT[k];
    }

// Walk actual layers and apply corresponding vT[k]

    int	is0, isN;

    GetLayerLimits( is0 = 0, isN );

    while( isN != -1 ) {

        int	z = vtil[is0].z;

        if( z >= zmin ) {

            for( int i = is0; i < isN; ++i ) {

                TAffine&	T = vtil[i].T;

                T = vT[z-zmin] * T;
            }
        }

        GetLayerLimits( is0 = isN, isN );
    }
}

/* --------------------------------------------------------------- */
/* BoundsPlus1 --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Update existing bounds B by including tile i.
//
void CTileSet::BoundsPlus1( DBox &B, const vector<Point> &cnr, int i )
{
    vector<Point>	c( 4 );
    memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
    vtil[i].T.Transform( c );

    for( int k = 0; k < 4; ++k ) {
        B.L = fmin( B.L, c[k].x );
        B.R = fmax( B.R, c[k].x );
        B.B = fmin( B.B, c[k].y );
        B.T = fmax( B.T, c[k].y );
    }
}

/* --------------------------------------------------------------- */
/* LayerBounds --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
void CTileSet::LayerBounds( DBox &B, int is0, int isN )
{
    vector<Point>	cnr;
    Set4Corners( cnr, gW, gH );

    B.L = BIGD, B.R = -BIGD,
    B.B = BIGD, B.T = -BIGD;

    for( int i = is0; i < isN; ++i )
        BoundsPlus1( B, cnr, i );
}

/* --------------------------------------------------------------- */
/* AllBounds ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::AllBounds( DBox &B )
{
    vector<Point>	cnr;
    Set4Corners( cnr, gW, gH );

    B.L = BIGD, B.R = -BIGD,
    B.B = BIGD, B.T = -BIGD;

    int	nt = vtil.size();

    for( int i = 0; i < nt; ++i )
        BoundsPlus1( B, cnr, i );
}

/* --------------------------------------------------------------- */
/* Reframe ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given result of AllBounds, calculate a neat BBox for whole
// stack and translate all transforms to reference new origin.
//
void CTileSet::Reframe( DBox &B )
{
    int	nt = vtil.size();

    for( int i = 0; i < nt; ++i )
        vtil[i].T.AddXY( -B.L, -B.B );

    B.R = ceil( B.R - B.L + 1 );
    B.T = ceil( B.T - B.B + 1 );
    B.L = 0;
    B.B = 0;
}

/* --------------------------------------------------------------- */
/* LayerAssignR -------------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::LayerAssignR( int is0, int isN, const DBox &B )
{
    Point	C( (B.L + B.R) / 2, (B.B + B.T) / 2 );

    for( int i = is0; i < isN; ++i ) {

        CUTile&	U = vtil[i];
        Point	cnr;	// (0,0)

        U.T.Transform( cnr );
        vaux[U.ix].r = cnr.DistSqr( C );
    }
}

/* --------------------------------------------------------------- */
/* ABOlap -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// This is an ABOlap helper class for ordering vertices by angle

class OrdVert {

public:
    double	a;
    vertex	v;

public:
    OrdVert( const vertex& v ) : v(v) {};

    bool operator < ( const OrdVert &rhs ) const
        {return a < rhs.a;};
};


// Return area(intersection) / (gW*gH).
//
// We construct the polygon enclosing the intersection in
// b-space. Its vertices comprise the set of rectangle edge
// crossings + any rectangle vertices that lie interior to
// the other.
//
// Once we've collected the vertices we order them by angle
// about their common centroid so they can be assembled into
// directed and ordered line segments for the area calculator.
//
double CTileSet::ABOlap( int a, int b, const TAffine *Tab )
{
// Quick proximity check

    TAffine	T;	// map a->b
    Point	pb( gW/2, gH/2 ),
            pa = pb;

    if( Tab )
        T = *Tab;
    else
        T.FromAToB( vtil[a].T, vtil[b].T );

    T.Transform( pa );

    if( pa.DistSqr( pb ) >= gW*gW + gH*gH )
        return 0.0;

// Detailed calculation

    vector<vertex>	va( 4 ), vb( 4 );

// Set b vertices (b-system) CCW ordered for LeftSide calcs

    vb[0] = vertex( 0   , 0 );
    vb[1] = vertex( gW-1, 0 );
    vb[2] = vertex( gW-1, gH-1 );
    vb[3] = vertex( 0   , gH-1 );

// Set a vertices (b-system) CCW ordered for LeftSide calcs

    {
        vector<Point>	c( 4 );
        Set4Corners( c, gW, gH );
        T.Transform( c );

        for( int i = 0; i < 4; ++i ) {

            va[i].x = int(c[i].x);
            va[i].y = int(c[i].y);

            // make close coords same

            for( int j = 0; j < i; ++j ) {

                if( iabs( va[i].x - va[j].x ) <= 2 )
                    va[i].x = va[j].x;

                if( iabs( va[i].y - va[j].y ) <= 2 )
                    va[i].y = va[j].y;
            }
        }
    }

// Gather all side intersections and all
// corners of A within B and vice versa.

    vector<OrdVert>	verts;

    for( int i = 0; i < 4; ++i ) {

        // count LeftSide successes for corner [i]

        int	ainb = 0,
            bina = 0;

        for( int j = 0; j < 4; ++j ) {

            // side crosses

            vertex	pi, pj;
            int		n;

            n = ClosedSegIsects( pi, pj,
                    va[i], va[(i+1)%4],
                    vb[j], vb[(j+1)%4] );

            switch( n ) {
                case 2:
                    verts.push_back( OrdVert( pj ) );
                case 1:
                    verts.push_back( OrdVert( pi ) );
                default:
                    break;
            }

            // corner counts

            ainb += LeftSide( vb[j], vb[(j+1)%4], va[i] );
            bina += LeftSide( va[j], va[(j+1)%4], vb[i] );
        }

        if( ainb == 4 )
            verts.push_back( OrdVert( va[i] ) );

        if( bina == 4 )
            verts.push_back( OrdVert( vb[i] ) );
    }

// Intersection?

    int	nv = verts.size();

    if( nv < 3 )
        return 0.0;

// Centroid

    Point	C;

    for( int i = 0; i < nv; ++i ) {
        C.x += verts[i].v.x;
        C.y += verts[i].v.y;
    }

    C.x /= nv;
    C.y /= nv;

// Set OrdVert angles w.r.t centroid,
// and sort by angle.

    for( int i = 0; i < nv; ++i ) {

        OrdVert&	V = verts[i];

        V.a = atan2( V.v.y - C.y, V.v.x - C.x );
    }

    sort( verts.begin(), verts.end() );

// Remove any duplicate verts

    for( int i = nv - 1; i > 0; --i ) {

        for( int j = 0; j < i; ++j ) {

            if( verts[i].v == verts[j].v ||
                verts[i].a == verts[j].a ) {

                verts.erase( verts.begin() + i );
                --nv;
                break;
            }
        }
    }

// Assemble vertices into polygon

    vector<lineseg>	pgon;

    for( int i = 0; i < nv; ++i )
        pgon.push_back( lineseg( verts[i].v, verts[(i+1)%nv].v ) );

// Return area

    return AreaOfPolygon( pgon ) / (gW * gH);
}

/* --------------------------------------------------------------- */
/* WriteTrakEM2Layer --------------------------------------------- */
/* --------------------------------------------------------------- */

// - xmltype is an ImagePlus pixel type code:
//		AUTO		= -1
//		GRAY8		= 0
//		GRAY16		= 1
//		GRAY32		= 2
//		COLOR_256	= 3
//		COLOR_RGB	= 4
//
// - xmlmin, xmlmax are used iff xmlmin < xmlmax. To let TrakEM2
// autoscale, set both to zero.
//
void CTileSet::WriteTrakEM2Layer(
    FILE*	f,
    int		&oid,
    int		xmltype,
    int		xmlmin,
    int		xmlmax,
    int		is0,
    int		isN )
{
// Layer prologue

    fprintf( f,
    "\t\t<t2_layer\n"
    "\t\t\toid=\"%d\"\n"
    "\t\t\tthickness=\"0\"\n"
    "\t\t\tz=\"%d\"\n"
    "\t\t>\n",
    oid++, vtil[is0].z );

// Tiles

    vector<int> order;

    isN = GetOrder_id( order, is0, isN );

    for( int i = 0; i < isN; ++i ) {

        const CUTile&	U = vtil[order[i]];
        char			title[128];

        if( U.col != -999 ) {

            sprintf( title, "%d.%d-1_%d.%d.%d",
                U.z, U.id, U.col, U.row, U.cam );
        }
        else
            sprintf( title, "%d.%d-1", U.z, U.id );

        fprintf( f,
        "\t\t\t<t2_patch\n"
        "\t\t\t\toid=\"%d\"\n"
        "\t\t\t\twidth=\"%d\"\n"
        "\t\t\t\theight=\"%d\"\n"
        "\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
        "\t\t\t\ttitle=\"%s\"\n"
        "\t\t\t\ttype=\"%d\"\n"
        "\t\t\t\tfile_path=\"%s\"\n"
        "\t\t\t\to_width=\"%d\"\n"
        "\t\t\t\to_height=\"%d\"\n",
        oid++, gW, gH,
        U.T.t[0], U.T.t[3], U.T.t[1], U.T.t[4], U.T.t[2], U.T.t[5],
        title, xmltype, U.name.c_str(), gW, gH );

        if( xmlmin < xmlmax ) {

            fprintf( f,
            "\t\t\t\tmin=\"%d\"\n"
            "\t\t\t\tmax=\"%d\"\n"
            "\t\t\t/>\n",
            xmlmin, xmlmax );
        }
        else
            fprintf( f, "\t\t\t/>\n" );
    }

// Layer epilogue

    fprintf( f, "\t\t</t2_layer>\n" );
}

/* --------------------------------------------------------------- */
/* WriteTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
// - BBox B previously obtained from AllBounds, Reframe.
//
// - xmltype is an ImagePlus pixel type code:
//		AUTO		= -1
//		GRAY8		= 0
//		GRAY16		= 1
//		GRAY32		= 2
//		COLOR_256	= 3
//		COLOR_RGB	= 4
//
// - xmlmin, xmlmax are used iff xmlmin < xmlmax. To let TrakEM2
// autoscale, set both to zero.
//
void CTileSet::WriteTrakEM2(
    const char	*path,
    DBox		&B,
    int			xmltype,
    int			xmlmin,
    int			xmlmax )
{
// Open file

    FILE	*f = FileOpenOrDie( path, "w" );

// Prologue + bounds

    int	oid = 3;

    fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

    TrakEM2WriteDTD( f );

    fprintf( f, "<trakem2>\n" );

    fprintf( f,
    "\t<project\n"
    "\t\tid=\"0\"\n"
    "\t\ttitle=\"Project\"\n"
    "\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
    "\t\tn_mipmap_threads=\"8\"\n"
    "\t/>\n" );

    fprintf( f,
    "\t<t2_layer_set\n"
    "\t\toid=\"%d\"\n"
    "\t\ttransform=\"matrix(1,0,0,1,0,0)\"\n"
    "\t\ttitle=\"Top level\"\n"
    "\t\tlayer_width=\"%.2f\"\n"
    "\t\tlayer_height=\"%.2f\"\n"
    "\t>\n",
    oid++, B.R, B.T );

// Layers

    int	is0, isN;

    GetLayerLimits( is0 = 0, isN );

    while( isN != -1 ) {

        WriteTrakEM2Layer( f, oid, xmltype, xmlmin, xmlmax, is0, isN );
        GetLayerLimits( is0 = isN, isN );
    }

// Epilogue

    fprintf( f, "\t</t2_layer_set>\n" );
    fprintf( f, "</trakem2>\n" );
    fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteTrakEM2_EZ ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
// - xmltype is an ImagePlus pixel type code:
//		AUTO		= -1
//		GRAY8		= 0
//		GRAY16		= 1
//		GRAY32		= 2
//		COLOR_256	= 3
//		COLOR_RGB	= 4
//
// - xmlmin, xmlmax are used iff xmlmin < xmlmax. To let TrakEM2
// autoscale, set both to zero.
//
void CTileSet::WriteTrakEM2_EZ(
    const char	*path,
    int			xmltype,
    int			xmlmin,
    int			xmlmax )
{
    DBox	B;

    AllBounds( B );
    Reframe( B );
    WriteTrakEM2( path, B, xmltype, xmlmin, xmlmax );
}

/* --------------------------------------------------------------- */
/* WriteTileToImage ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Tiles must already be sorted by z (subsort doesn't matter).
//
// Optionally create subfolder: topdir/z.
// Optionally create subfolder: topdir/z/nmrc.
//
// Write one TileToImage.txt file for given layer.
//
void CTileSet::WriteTileToImage(
    const char	*topdir,
    bool		create_zdir,
    bool		create_nmrcdir,
    int			is0,
    int			isN )
{
// Manage dirs

    char	path[2048];
    FILE	*f;
    int		len;

    len = sprintf( path, "%s/%d", topdir, vtil[is0].z );

    if( create_zdir )
        DskCreateDir( path, flog );

    if( create_nmrcdir ) {

        sprintf( path + len, "/nmrc" );
        DskCreateDir( path, flog );
    }

// Open file

    sprintf( path + len, "/TileToImage.txt" );
    f = FileOpenOrDie( path, "w", flog );

// Header

    fprintf( f, "ID\tT0\tT1\tX\tT3\tT4\tY\tCol\tRow\tCam\tPath\n" );

// Write sorted entries

    vector<int> order;

    isN = GetOrder_id( order, is0, isN );

    for( int i = 0; i < isN; ++i ) {

        const CUTile&	U = vtil[order[i]];
        const double*	T = U.T.t;

        fprintf( f,
            "%d"
            "\t%f\t%f\t%f\t%f\t%f\t%f"
            "\t%d\t%d\t%d\t%s\n",
            U.id,
            T[0], T[1], T[2], T[3], T[4], T[5],
            U.col, U.row, U.cam, U.name.c_str() );
    }

    fclose( f );
}


