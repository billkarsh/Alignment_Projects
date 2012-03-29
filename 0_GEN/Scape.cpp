

#include	"Scape.h"
#include	"ImageIO.h"
#include	"Maths.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const double	BIGD = 1.0e30;

/* --------------------------------------------------------------- */
/* AdjustBounds -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void AdjustBounds(
	uint32			&ws,
	uint32			&hs,
	double			&w0,
	double			&h0,
	vector<ScpTile>	&vTile,
	int				wi,
	int				hi,
	double			scale,
	int				szmult )
{
// maximum extents over all image corners

	double	xmin, xmax, ymin, ymax;
	int		nt = vTile.size();

	xmin =  BIGD;
	xmax = -BIGD;
	ymin =  BIGD;
	ymax = -BIGD;

	for( int i = 0; i < nt; ++i ) {

		vector<Point>	cnr( 4 );

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( wi-1, 0.0 );
		cnr[2] = Point( wi-1, hi-1 );
		cnr[3] = Point(  0.0, hi-1 );

		vTile[i].t2g.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			xmin = fmin( xmin, cnr[k].x );
			xmax = fmax( xmax, cnr[k].x );
			ymin = fmin( ymin, cnr[k].y );
			ymax = fmax( ymax, cnr[k].y );
		}
	}

// scale, and expand out to integer bounds

	w0 = floor( xmin * scale );
	h0 = floor( ymin * scale );
	ws = (int)ceil( (xmax + 1) * scale - w0 );
	hs = (int)ceil( (ymax + 1) * scale - h0 );

// ensure dims divisible by szmult

	int	rem;

	if( rem = ws % szmult )
		ws += szmult - rem;

	if( rem = hs % szmult )
		hs += szmult - rem;

// propagate new dims

	TForm	A, t;

	A.NUSetScl( scale );

	for( int i = 0; i < nt; ++i ) {

		TForm	&T = vTile[i].t2g;

		MultiplyTrans( T, A, t = T );
		T.AddXY( -w0, -h0 );
	}
}

/* --------------------------------------------------------------- */
/* ScanLims ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ScanLims(
	int			&x0,
	int			&xL,
	int			&y0,
	int			&yL,
	int			ws,
	int			hs,
	const TForm	&T,
	int			wi,
	int			hi )
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
/* Paint --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Paint(
	uint8			*scp,
	uint32			ws,
	uint32			hs,
	vector<ScpTile>	&vTile,
	int				rmvedges,
	FILE*			flog )
{
	int		nt = vTile.size();

	for( int i = 0; i < nt; ++i ) {

		TForm	inv;
		int		x0, xL, y0, yL,
				wL, hL;
		uint32	w,  h;
		uint8*	src = Raster8FromAny(
						vTile[i].name.c_str(), w, h, flog );

		ScanLims( x0, xL, y0, yL, ws, hs, vTile[i].t2g, w, h );

		wL = w - (rmvedges != 0);
		hL = h - (rmvedges != 0);

		InvertTrans( inv, vTile[i].t2g );

		for( int iy = y0; iy < yL; ++iy ) {

			for( int ix = x0; ix < xL; ++ix ) {

				Point	p( ix, iy );

				inv.Transform( p );

				if( p.x >= 0 && p.x < wL &&
					p.y >= 0 && p.y < hL ) {

					scp[ix+ws*iy] =
					(int)SafeInterp( p.x, p.y, src, w, h );
				}
			}
		}

		RasterFree( src );
	}
}

/* --------------------------------------------------------------- */
/* Scape --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Allocate and return pointer to new montage built by painting
// the listed tiles, and return its dims (ws, hs), and the top-left
// of the bounding box (w0, h0).
//
// Return NULL if unsuccessful.
//
// (wi, hi)	- specify input tile size (same for all).
// scale	- for example, 0.25 reduces by 4X.
// szmult	- scape dims made divisible by szmult.
// bkval	- default scape value where no data.
// rmvedges	- don't paint tile's edgemost pixels (smoother result).
//
// Caller must dispose of scape with ImageIO::RasterFree().
//
uint8* Scape(
	uint32			&ws,
	uint32			&hs,
	double			&w0,
	double			&h0,
	vector<ScpTile>	&vTile,
	int				wi,
	int				hi,
	double			scale,
	int				szmult,
	int				bkval,
	int				rmvedges,
	FILE*			flog )
{
	AdjustBounds( ws, hs, w0, h0, vTile, wi, hi, scale, szmult );

	int		ns		= ws * hs;
	uint8	*scp	= (uint8*)RasterAlloc( ns );

	if( scp ) {
		memset( scp, bkval, ns );
		Paint( scp, ws, hs, vTile, rmvedges, flog );
	}
	else
		fprintf( flog, "Scape: Alloc failed (%d x %d).\n", ws, hs );

	return scp;
}


