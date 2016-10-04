

#pragma once


#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stack>
using namespace std;


/* --------------------------------------------------------------- */
/* Lines --------------------------------------------------------- */
/* --------------------------------------------------------------- */

double SegPointDist(
    int		x0,
    int		y0,
    int		x1,
    int		y1,
    int		x2,
    int		y2 );

/* --------------------------------------------------------------- */
/* Boxes --------------------------------------------------------- */
/* --------------------------------------------------------------- */

void BBoxFromPoints( IBox &B, const vector<Point> &pts );

void BBoxFromPoints( DBox &B, const vector<Point> &pts );

void BoxesFromShifts(
    IBox	&B1,
    IBox	&B2,
    int		w,
    int		h,
    int		w2,
    int		h2,
    int		x,
    int		y );

int TightestBBox( DBox &B, const vector<Point> &pts );

/* --------------------------------------------------------------- */
/* Regions ------------------------------------------------------- */
/* --------------------------------------------------------------- */

int Propagate(
    vector<double>	&v,
    int				w,
    int				h,
    int				first,
    double			tmin,
    double			tmax,
    double			set_to );

int Propagate(
    vector<Point>	&plist,
    vector<double>	&v,
    int				w,
    int				h,
    int				first,
    double			thresh,
    double			set_to );

int MapBlobRng(
    vector<uint8>			&map,
    const vector<double>	&I,
    int						w,
    int						h,
    int						first,
    double					tmin,
    double					tmax );

int MapBlobRng(
    vector<uint8>	&map,
    const uint8		*I,
    int				w,
    int				h,
    int				first,
    int				tmin,
    int				tmax );

int MapBlobVar(
    vector<uint8>	&map,
    const uint8		*I,
    int				w,
    int				h,
    int				first,
    int				size,
    int				tol );

void DilateMap1Pix( vector<uint8> &map, int w, int h );
void ErodeMap1Pix( vector<uint8> &map, int w, int h );

/* --------------------------------------------------------------- */
/* Polygons ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void PixelListFromPolygon(
    vector<Point>			&Plist,
    const vector<Point>		&Pgon );

/* --------------------------------------------------------------- */
/* Images -------------------------------------------------------- */
/* --------------------------------------------------------------- */

void ImageFromValuesAndPoints(
    vector<double>			&I,
    int						w,
    int						h,
    const vector<double>	&v,
    const vector<Point>		&plist,
    int						xmin,
    int						ymin );

void ValuesFromImageAndPoints(
    vector<double>			&v,
    const uint8*			raster,
    int						w,
    const vector<Point>		&plist );

void ValuesFromImageAndPoints(
    vector<double>			&v,
    double					&dvx,
    double					&dvy,
    const vector<double>	&raster,
    int						w,
    const vector<Point>		&plist,
    const vector<double>	&spv );

/* --------------------------------------------------------------- */
/* Meshes -------------------------------------------------------- */
/* --------------------------------------------------------------- */

class vertex {

public:
    int	x, y;	// coordinates
    int	orig;	// index in original array @@@ not used
    int	dir;	// compass-pt direction to next vertex

public:
    vertex()
    : x(0), y(0), orig(0), dir(0) {};

    vertex( int x, int y )
    : x(x), y(y), orig(0), dir(0) {};

    vertex( int x, int y, int dir )
    : x(x), y(y), orig(0), dir(dir) {};

    vertex( const vertex &v )
    : x(v.x), y(v.y), orig(v.orig), dir(v.dir) {};

    bool operator == ( const vertex &rhs ) const
        {return x == rhs.x && y == rhs.y;};

    bool operator < ( const vertex &rhs ) const
        {return x < rhs.x || (x == rhs.x && y < rhs.y);};

    int   IDistSqr( const vertex& rhs ) const;
    double DistSqr( const vertex& rhs ) const;
};


class lineseg {

public:
    vertex	v[2];

public:
    lineseg( const vertex &v1, const vertex &v2 )
        {v[0] = v1; v[1] = v2;};

    lineseg( int x0, int y0, int x1, int y1 )
        {v[0].x = x0; v[0].y = y0; v[1].x = x1; v[1].y = y1;};

    bool operator == ( const lineseg &rhs ) const
        {return v[0] == rhs.v[0] && v[1] == rhs.v[1];};

    double LenSqr() const
        {return v[1].DistSqr( v[0] );};
};


class triangle {

public:
    double	a[3][3];	// matrix turning point (x,y,1) into
                        //   a barycentric representation
                        //   (linear combination of vertices)
    int		v[3];		// indices into vector of control points

public:
    double Area( const vector<Point> &ctl ) const;
};


double SegPointDist(
    const vertex	&v0,
    const vertex	&v1,
    const vertex	&v2 );

bool LeftSide( const Point &a, const Point &b, const Point &c );
bool LeftSide( const vertex &a, const vertex &b, const vertex &c );

bool LinesCross(
    const vertex	&p1,
    const vertex	&p2,
    const vertex	&p3,
    const vertex	&p4 );

int ClosedSegIsects(
    vertex			&pi,
    vertex			&pj,
    const vertex	&p1,
    const vertex	&p2,
    vertex			p3,
    vertex			p4 );

bool OpenSegsCross(
    const vertex	&p1,
    const vertex	&p2,
    const vertex	&p3,
    const vertex	&p4 );

bool AnyCrossing(
    const vector<lineseg>	&s,
    const vertex			&a,
    const vertex			&b );

int CountCrossings(
    const vector<lineseg>	&s,
    const vertex			&a,
    const vertex			&b );

bool IsSubseg(
    stack<lineseg>	&stk,
    const lineseg	&e,
    const lineseg	&W );

bool InTriangle(
    const vertex	&va,
    const vertex	&vb,
    const vertex	&vc,
    const vertex	&p );

bool AnyInside(
    const vertex			&va,
    const vertex			&vb,
    const vertex			&vc,
    const vector<lineseg>	&s,
    const vector<vertex>	&ips );

double AreaOfTriangle(
    const Point		&v0,
    const Point		&v1,
    const Point		&v2 );

double AreaOfTriangle(
    const vertex	&v0,
    const vertex	&v1,
    const vertex	&v2 );

int BestTriangle(
    const vector<triangle>	&T,
    const vector<vertex>	&C,
    const Point				&p );

double AreaOfPolygon( const vector<lineseg> &edges );

void RemoveFromMap(
    vector<uint8>	&map,
    int				w,
    int				h,
    const vertex	&v,
    int				dist );


