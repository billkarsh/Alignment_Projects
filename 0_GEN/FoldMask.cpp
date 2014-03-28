

#include	"Disk.h"
#include	"FoldMask.h"
#include	"ImageIO.h"

#include	<stdlib.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* PrintFoldmapHisto --------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintFoldmapHisto( const uint8* fm, int w, int h )
{
	vector<int>	cts( 256, 0 );
	int			i, n = w * h;

	for( i = 0; i < n; ++i )
		++cts[fm[i]];

	for( i = 0; i < 256; ++i ) {

		if( cts[i] )
			printf( "Foldmask: value=%3d, count=%8d\n", i, cts[i] );
	}
}

/* --------------------------------------------------------------- */
/* GetFoldMask --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Load or create a foldmask (always full size).
//
uint8* GetFoldMask(
	const string		&idb,
	const PicSpecs		&P,
	const char			*forcepath,
	const vector<uint8>	&resmsk,
	CCropMask			*CM,
	int					wf,
	int					hf,
	bool				nofile,
	bool				transpose,
	bool				force1rgn )
{
	uint8*	mask;
	int		np = wf * hf;

	if( nofile ) {
		mask = (uint8*)RasterAlloc( np );
		memset( mask, 1, np );
	}
	else {

		Til2FM	t2f;

		if( !forcepath ) {
			IDBTil2FM( t2f, idb, P.z, P.id );
			forcepath = t2f.path.c_str();
		}

		uint32	_w, _h;

		mask = Raster8FromAny( forcepath, _w, _h, stdout, transpose );

		if( _w != wf || _h != hf ) {

			printf(
			"GetFoldMask: Maps different size than input.\n" );
			exit( 42 );
		}

		PrintFoldmapHisto( mask, wf, hf );

		// force one (non-fold) region

		if( force1rgn ) {

			for( int i = 0; i < np; ++i ) {

				if( mask[i] )
					mask[i] = 1;
			}
		}
	}

// optionally remove resin

	if( resmsk.size() == np ) {

		for( int i = 0; i < np; ++i ) {

			if( !resmsk[i] )
				mask[i] = 0;
		}
	}

// optionally crop mask borders

	if( CM && CM->IsFile( idb ) ) {

		IBox	B;
		CM->GetBox( B, P.t2i.cam );

		for( int i = 0; i < np; ++i ) {

			int	y = i / wf,
				x = i - wf * y;

			if( y < B.B || y >= B.T ||
				x < B.L || x >= B.R ) {

				mask[i] = 0;
			}
		}

		printf( "Crop z %d id %d cam %d to x[%d %d) y[%d %d)\n",
		P.z, P.id, P.t2i.cam, B.L, B.R, B.B, B.T );
	}

	return mask;
}

/* --------------------------------------------------------------- */
/* SetWithinSectionBorders --------------------------------------- */
/* --------------------------------------------------------------- */

void SetWithinSectionBorders( uint8* foldMask, int wf, int hf )
{
	const int border = 100;

	for( int y = 0; y < hf; ++y ) {

		for( int x = 0; x < wf; ++x ) {

			int dt	= hf-1-y;	// distance from top
			int dr	= wf-1-x;	// distance from right
			int pix	= 0;

			if( x < border && y > x && dt > x && dr > x )
				pix = 1;	// left side

			if( y < border && x > y && dt > y && dr > y )
				pix = 2;	// bottom

			if( dt < border && x > dt && y > dt && dr > dt )
				pix = 3;	// top

			if( dr < border && x > dr && y > dr && dt > dr )
				pix = 4;	// right side

			foldMask[x + wf*y] = pix;
		}
	}
}

/* --------------------------------------------------------------- */
/* SetBoundsAndColors -------------------------------------------- */
/* --------------------------------------------------------------- */

// Given vector of ConnRegion whose points are determined...
// (1) Calculate each region's bounds
// (2) Color the folmask raster with region ids.
//
void SetBoundsAndColors(
	vector<ConnRegion>	&cr,
	uint8*				foldMask,
	int					wf,
	int					hf )
{
	int	np = wf * hf,
		nr = cr.size();

// zero foldmask (the fold value)
	memset( foldMask, 0, np );

// print region header
	printf( "SetBnd&Clr: Big enough regions = %d:\n", nr );

// for each region
	for( int ir = 0; ir < nr; ++ir ) {

		ConnRegion&	C = cr[ir];

		np = C.pts.size();

		// init bounds
		C.B.L	= BIG;
		C.B.B	= BIG;
		C.B.R	= -BIG;
		C.B.T	= -BIG;

		// assign region id
		C.id = ir + 1;

		// for each point
		for( int ip = 0; ip < np; ++ip ) {

			int	x = (int)C.pts[ip].x,
				y = (int)C.pts[ip].y;

			// update bounds
			if( x < C.B.L )
				C.B.L = x;
			else if( x > C.B.R )
				C.B.R = x;

			if( y < C.B.B )
				C.B.B = y;
			else if( y > C.B.T )
				C.B.T = y;

			// color mask pixel (these are non-zero)
			foldMask[x + wf*y] = ir + 1;
		}

		// print this region
		printf( "\tid=%2d, pts=%8d, x=[%4d %4d], y=[%4d %4d].\n",
			C.id, np, C.B.L, C.B.R, C.B.B, C.B.T );
	}
}


// Given vector of ConnRegion whose points are determined...
// (1) Calculate each region's bounds
// (2) Color the folmask raster with region ids...
//
// In this case, each point to be colored in the foldmask is
// enlarged to a (D+1)x(D+1) box and we color all the points
// of each such box that are also above thresh in image valid.
//
void SetBoundsAndColors(
	vector<ConnRegion>		&cr,
	uint8*					foldMask,
	const vector<double>	&valid,
	int						wf,
	int						hf,
	double					thresh,
	int						D )
{
	int	np = wf * hf,
		nr = cr.size();

// zero foldmask (the fold value)
	memset( foldMask, 0, np );

// print region header
	printf( "SetBnd&Clr: Big enough regions = %d:\n", nr );

// for each region
	for( int ir = 0; ir < nr; ++ir ) {

		ConnRegion&	C = cr[ir];

		np = C.pts.size();

		// init bounds
		C.B.L	= BIG;
		C.B.B	= BIG;
		C.B.R	= -BIG;
		C.B.T	= -BIG;

		// assign region id
		C.id = ir + 1;

		// for each point
		for( int ip = 0; ip < np; ++ip ) {

			int	x   = (int)C.pts[ip].x,
				y   = (int)C.pts[ip].y,
				xlo = max(    0, x - D ),
				xhi = min( wf-1, x + D ),
				ylo = max(    0, y - D ),
				yhi = min( hf-1, y + D ),
				clr = ir + 1,
				idx;

			// update bounds
			if( xlo < C.B.L )
				C.B.L = xlo;
			else if( xhi > C.B.R )
				C.B.R = xhi;

			if( ylo < C.B.B )
				C.B.B = ylo;
			else if( yhi > C.B.T )
				C.B.T = yhi;

			// color mask pixels

			for( y = ylo; y <= yhi; ++y ) {

				for( x = xlo; x <= xhi; ++x ) {

					if( valid[idx = x + wf*y] > thresh )
						foldMask[idx] = clr;
				}
			}
		}

		// print this region
		printf( "\tid=%2d, pts=%8d, x=[%4d %4d], y=[%4d %4d].\n",
			C.id, np, C.B.L, C.B.R, C.B.B, C.B.T );
	}
}

/* --------------------------------------------------------------- */
/* ConnRgnsFromFoldMask ------------------------------------------ */
/* --------------------------------------------------------------- */

// Scan given foldMask having pixels {0=fold, 1=rgn1, ...} and
// make an entry for each connected region in (cr). For each cr
// we fill in fields: {pts, B, id}.
//
// We will not create entries for the fold or for any region
// whose point count is below (minpts).
//
void ConnRgnsFromFoldMask(
	vector<ConnRegion>	&cr,
	const uint8*		foldMask,
	int					wf,
	int					hf,
	int					scale,
	uint32				minpts,
	FILE				*flog )
{
// Pass 1: gather ids and size info for cr and pts vectors

	vector<uint32>	szrgn( 256, 0 );
	int				N = wf * hf;
	int				max_id = 0;

	for( int i = 0; i < N; ++i )
		++szrgn[foldMask[i]];

// Find the highest occurring region id

	for( int i = 255; i > 0; --i ) {

		if( szrgn[i] ) {
			max_id = i;
			break;
		}
	}

// Report what will be included and excluded

	uint32	n_inc = 0, t_inc = 0, t_exc = 0;

	for( int i = 1; i <= max_id; ++i ) {

		if( szrgn[i] >= minpts ) {
			++n_inc;
			t_inc += szrgn[i];
		}
		else {
			t_exc += szrgn[i];
			szrgn[i] = 0;	// mark for skip
		}
	}

	fprintf( flog,
	"ConnRegion: FoldMask w=%d, h=%d, area=%d, scale=%d\n",
	wf, hf, N, scale );

	fprintf( flog,
	"ConnRegion: Included rgns=%3d, area=%8d, %%=%6.2f\n",
	n_inc, t_inc, 100.0*t_inc / N );

	fprintf( flog,
	"ConnRegion: Excluded rgns=%3d, area=%8d, %%=%6.2f\n",
	max_id - n_inc, t_exc, 100.0*t_exc / N );

	fprintf( flog,
	"ConnRegion:    Folds rgns=%3d, area=%8d, %%=%6.2f\n",
	(szrgn[0] != 0), szrgn[0], 100.0*szrgn[0] / N );

	szrgn[0] = 0;	// mark for skip

// Size the vectors (pts cpacity ample for any scaling)
// Also set the cr.id fields here.

	cr.clear();

	if( !max_id )
		return;

	cr.resize( max_id );

	for( int id = 1; id <= max_id; ++id ) {

		ConnRegion&	C = cr[id - 1];

		C.pts.resize( szrgn[id] );
		C.id = id;
	}

// Pass 2: gather unique scaled points and update BBoxes

	vector<uint32>	pcnt( max_id + 1, 0 );
	vector<uint8>	seen( N / (scale*scale), 0 );
	int				ws = wf / scale;

	for( int i = 0; i < N; ++i ) {

		int	id = foldMask[i];

		if( !szrgn[id] )
			continue;

		int	Y =  i / wf,
			y =  Y / scale,
			x = (i - wf * Y) / scale;

		if( seen[x + ws*y] )
			continue;

		ConnRegion&	C = cr[id - 1];

		seen[x + ws*y] = 1;

		if( x < C.B.L )
			C.B.L = x;
		else if( x > C.B.R )
			C.B.R = x;

		if( y < C.B.B )
			C.B.B = y;
		else if( y > C.B.T )
			C.B.T = y;

		C.pts[pcnt[id]++] = Point( x, y );
	}

// Lastly, remove empty cr, or if keeping, set pts actual size

	for( int i = 0; i < cr.size(); ++i ) {

		int	np = pcnt[cr[i].id];

		if( np )
			cr[i].pts.resize( np );
		else
			cr.erase( cr.begin() + i );
	}
}

/* --------------------------------------------------------------- */
/* ConnRgnForce1 ------------------------------------------------- */
/* --------------------------------------------------------------- */

void ConnRgnForce1( vector<ConnRegion> &cr, int ws, int hs )
{
	cr.clear();
	cr.resize( 1 );

	ConnRegion&	C = cr[0];

	MakeZeroBasedPoints( C.pts, ws, hs );

	C.B.L	= 0;
	C.B.R	= ws - 1;
	C.B.B	= 0;
	C.B.T	= hs - 1;

	C.id	= 1;
}


