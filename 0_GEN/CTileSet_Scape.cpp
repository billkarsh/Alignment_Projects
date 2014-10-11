

#include	"CTileSet.h"
#include	"EZThreads.h"
#include	"ImageIO.h"
#include	"Maths.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Scape_AdjustBounds -------------------------------------------- */
/* --------------------------------------------------------------- */

void CTileSet::Scape_AdjustBounds(
	uint32				&ws,
	uint32				&hs,
	double				&x0,
	double				&y0,
	vector<TAffine>		&vTadj,
	const vector<int>	&vid,
	double				scale,
	int					szmult ) const
{
// maximum extents over all image corners

	vector<Point>	cnr;
	double			xmin, xmax, ymin, ymax;
	int				nt = vid.size();

	xmin =  BIGD;
	xmax = -BIGD;
	ymin =  BIGD;
	ymax = -BIGD;

	Set4Corners( cnr, gW, gH );

	for( int i = 0; i < nt; ++i ) {

		vector<Point>	c( 4 );
		memcpy( &c[0], &cnr[0], 4*sizeof(Point) );
		vtil[vid[i]].T.Transform( c );

		for( int k = 0; k < 4; ++k ) {
			xmin = fmin( xmin, c[k].x );
			xmax = fmax( xmax, c[k].x );
			ymin = fmin( ymin, c[k].y );
			ymax = fmax( ymax, c[k].y );
		}
	}

// scale, and expand out to integer bounds

	x0 = xmin * scale;
	y0 = ymin * scale;
	ws = (int)ceil( (xmax - xmin + 1) * scale );
	hs = (int)ceil( (ymax - ymin + 1) * scale );

// ensure dims divisible by szmult

	int	rem;

	if( rem = ws % szmult )
		ws += szmult - rem;

	if( rem = hs % szmult )
		hs += szmult - rem;

// propagate new dims

	TAffine	A;

	A.NUSetScl( scale );

	for( int i = 0; i < nt; ++i ) {

		TAffine	T;

		T = A * vtil[vid[i]].T;
		T.AddXY( -x0, -y0 );
		vTadj.push_back( T );
	}
}

/* --------------------------------------------------------------- */
/* ScanLims ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Fill in the range of coords [x0,xL); [y0,yL) in the scape
// that the given tile will occupy. This tells the painter
// which region to fill in.
//
// We expand the tile by one pixel to be conservative, and
// transform the tile bounds to scape coords. In no case are
// the scape coords to exceed [0,ws); [0,hs).
//
static void ScanLims(
	int				&x0,
	int				&xL,
	int				&y0,
	int				&yL,
	int				ws,
	int				hs,
	const TAffine	&T,
	int				wi,
	int				hi )
{
	double	xmin, xmax, ymin, ymax;

	xmin =  BIGD;
	xmax = -BIGD;
	ymin =  BIGD;
	ymax = -BIGD;

	vector<Point>	cnr( 4 );

// generous box (outset 1 pixel) for scanning

	cnr[0] = Point( -1.0, -1.0 );
	cnr[1] = Point(   wi, -1.0 );
	cnr[2] = Point(   wi,   hi );
	cnr[3] = Point( -1.0,   hi );

	T.Transform( cnr );

	for( int k = 0; k < 4; ++k ) {
		xmin = fmin( xmin, cnr[k].x );
		xmax = fmax( xmax, cnr[k].x );
		ymin = fmin( ymin, cnr[k].y );
		ymax = fmax( ymax, cnr[k].y );
	}

	x0 = max( 0, (int)floor( xmin ) );
	y0 = max( 0, (int)floor( ymin ) );
	xL = min( ws, (int)ceil( xmax ) );
	yL = min( hs, (int)ceil( ymax ) );
}

/* --------------------------------------------------------------- */
/* Downsample ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Downsample( uint8 *ras, int &w, int &h, int iscl )
{
	int	n  = iscl * iscl,
		ws = (int)ceil( (double)w / iscl ),
		hs = (int)ceil( (double)h / iscl ),
		w0 = w,
		xo, yo, xr, yr;

	yr = iscl - h % iscl;
	xr = iscl - w % iscl;

	for( int iy = 0; iy < hs; ++iy ) {

		yo = 0;

		if( iy == hs - 1 )
			yo = yr;

		for( int ix = 0; ix < ws; ++ix ) {

			double	sum = 0.0;

			xo = 0;

			if( ix == ws - 1 )
				xo = xr;

			for( int dy = 0; dy < iscl; ++dy ) {

				for( int dx = 0; dx < iscl; ++dx )
					sum += ras[ix*iscl-xo+dx + w0*(iy*iscl-yo+dy)];
			}

			ras[ix+ws*iy] = int(sum / n);
		}
	}

	w = ws;
	h = hs;
}

/* --------------------------------------------------------------- */
/* NormRas ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Force image mean to 127 and sd to sdnorm.
//
static void NormRas( uint8 *r, int w, int h, int lgord, int sdnorm )
{
// flatfield & convert to doubles

	int				n = w * h;
	vector<double>	v;

	LegPolyFlatten( v, r, w, h, lgord );

// rescale to mean=127, sd=sdnorm

	for( int i = 0; i < n; ++i ) {

		int	pix = 127 + int(v[i] * sdnorm);

		if( pix < 0 )
			pix = 0;
		else if( pix > 255 )
			pix = 255;

		r[i] = pix;
	}
}

/* --------------------------------------------------------------- */
/* Scape_Paint --------------------------------------------------- */
/* --------------------------------------------------------------- */

class CPaintPrms {
// Parameters for _Scape_Paint()
public:
	uint8					*scp;
	uint32					ws;
	uint32					hs;
	const vector<TAffine>	&vTadj;
	const vector<int>		&vid;
	int						iscl;
	int						bkval;
	int						lgord;
	int						sdnorm;
	bool					resmask;
public:
	CPaintPrms(
		uint8					*scp,
		uint32					ws,
		uint32					hs,
		const vector<TAffine>	&vTadj,
		const vector<int>		&vid,
		int						iscl,
		int						bkval,
		int						lgord,
		int						sdnorm,
		bool					resmask )
	: scp(scp), ws(ws), hs(hs),
	vTadj(vTadj), vid(vid), iscl(iscl), bkval(bkval),
	lgord(lgord), sdnorm(sdnorm), resmask(resmask)
	{};
};

class CThrdat {
public:
	int	i0, ilim;
};

static const CTileSet	*ME;
static const CPaintPrms	*GP;
static vector<CThrdat>	vthr;

void* _Scape_Paint( void *ithr )
{
	CThrdat	&me = vthr[(long)ithr];

	for( int i = me.i0; i < me.ilim; ++i ) {

		vector<uint8>	msk;
		uint8*			src;
		TAffine			inv;
		uint32			w,  h;
		int				x0, xL, y0, yL,
						wL, hL,
						wi, hi;

		src = Raster8FromAny(
				ME->vtil[GP->vid[i]].name.c_str(),
				w, h, ME->flog );

		if( GP->resmask )
			ResinMask8( msk, src, w, h, false );

		if( GP->sdnorm > 0 )
			NormRas( src, w, h, GP->lgord, GP->sdnorm );

		if( GP->resmask ) {

			int	n = w * h;

			for( int j = 0; j < n; ++j ) {
				if( !msk[j] )
					src[j] = GP->bkval;
			}
		}

		ScanLims( x0, xL, y0, yL, GP->ws, GP->hs, GP->vTadj[i], w, h );
		wi = w;
		hi = h;

		inv.InverseOf( GP->vTadj[i] );

		if( GP->iscl > 1 ) {	// Scaling down

			// actually downsample src image
			Downsample( src, wi, hi, GP->iscl );

			// and point at the new pixels
			TAffine	A;
			A.NUSetScl( 1.0/GP->iscl );
			inv = A * inv;
		}

		wL = wi - 1;
		hL = hi - 1;

		for( int iy = y0; iy < yL; ++iy ) {

			for( int ix = x0; ix < xL; ++ix ) {

				Point	p( ix, iy );

				inv.Transform( p );

				if( p.x >= 0 && p.x < wL &&
					p.y >= 0 && p.y < hL ) {

					int	pix =
					(int)SafeInterp( p.x, p.y, src, wi, hi );

					if( pix != GP->bkval )
						GP->scp[ix+GP->ws*iy] = pix;
				}
			}
		}

		RasterFree( src );
	}

	return NULL;
}


void CTileSet::Scape_PaintTH( int nthr ) const
{
	int	nt = GP->vid.size(),	// tiles total
		nb;						// tiles per thread

	if( nthr > nt )
		nthr = nt;

	nb = nt / nthr;

	vthr.resize( nthr );

	vthr[0].i0		= 0;
	vthr[0].ilim	= nb;

	for( int i = 1; i < nthr; ++i ) {
		CThrdat	&C = vthr[i];
		C.i0	= vthr[i-1].ilim;
		C.ilim	= (i == nthr-1 ? nt : C.i0 + nb);
	}

	if( !EZThreads( _Scape_Paint, nthr, 2, "_Scape_Paint", flog ) )
		exit( 42 );

	vthr.clear();
}

/* --------------------------------------------------------------- */
/* Scape --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Allocate and return pointer to new montage built by painting
// the listed tiles, and return its dims (ws, hs), and the top-left
// of the scaled bounding box (x0, y0).
//
// Return NULL if unsuccessful.
//
// scale	- for example, 0.25 reduces by 4X.
// szmult	- scape dims made divisible by szmult.
// bkval	- default scape value where no data.
// lgord	- Legendre poly max order.
// sdnorm	- if > 0, image normalized to mean=127, sd=sdnorm.
// resmask	- mask resin areas.
// nthr		- thread count.
//
// Caller must dispose of scape with ImageIO::RasterFree().
//
uint8* CTileSet::Scape(
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
	int					nthr ) const
{
	if( !vid.size() ) {
		fprintf( flog, "Scape: Empty tile list.\n" );
		return NULL;
	}

	vector<TAffine>	vTadj;

	Scape_AdjustBounds( ws, hs, x0, y0, vTadj, vid, scale, szmult );

	int		ns		= ws * hs;
	uint8	*scp	= (uint8*)RasterAlloc( ns );

	if( scp ) {

		if( sdnorm > 0 )
			bkval = 127;

		memset( scp, bkval, ns );

		ME = this;
		GP = new CPaintPrms( scp, ws, hs,
					vTadj, vid, int(1/scale),
					bkval, lgord, sdnorm, resmask );
		Scape_PaintTH( nthr );
		delete GP;
	}
	else
		fprintf( flog, "Scape: Alloc failed (%d x %d).\n", ws, hs );

	return scp;
}


