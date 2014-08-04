

#include	"Cmdline.h"
#include	"Disk.h"
#include	"FoldMask.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Geometry.h"
#include	"CPicBase.h"
#include	"Memory.h"

#include	"numerical_recipes.h"


/* --------------------------------------------------------------- */
/* CArgs_tiny ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_tiny {

public:
	double		fmTOverride;
	const char	*infile,
				*nmrc,
				*fm,
				*fmd;
	int			Z,
				ID,
				D,
				minarea;
	bool		nomasks,
				oneregion,
				transpose,
				dumphist;

public:
	CArgs_tiny()
	{
		fmTOverride		= 0.0;
		infile			= NULL;
		nmrc			= NULL;
		fm				= NULL;
		fmd				= NULL;
		Z				= -1;
		ID				= -1;
		D				= -1;
		minarea			= 90000;
		nomasks			= false;
		oneregion		= false;
		transpose		= false;
		dumphist		= false;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_tiny	gArgs;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_tiny::SetCmdLine( int argc, char* argv[] )
{
// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: tiny <z> <id> <tif-file> <fold-mask-file> [option]\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( Z == -1 )
				Z = atoi( argv[i] );
			else if( ID == -1 )
				ID = atoi( argv[i] );
			else
				infile = argv[i];
		}
		else if( GetArgStr( nmrc, "-nmrc=", argv[i] ) )
			;
		else if( GetArgStr( fm, "-fm=", argv[i] ) )
			;
		else if( GetArgStr( fmd, "-fmd=", argv[i] ) )
			;
		else if( GetArg( &fmTOverride, "-fmto=%lf", argv[i] ) ) {

			printf( "Fold Mask Threshold overridden, now %f.\n",
			fmTOverride );
		}
		else if( GetArg( &D, "-D=%d", argv[i] ) ) {

			printf( "Fold radius overridden, now %d.\n", D );
		}
		else if( GetArg( &minarea, "-minarea=%d", argv[i] ) )
			;
		else if( IsArg( "-nf", argv[i] ) ) {

			nomasks = true;
			printf( "Not generating mask files.\n" );
		}
		else if( IsArg( "-one", argv[i] ) ) {

			oneregion = true;
			printf( "Using just one region.\n" );
		}
		else if( IsArg( "-h", argv[i] ) ) {

			dumphist = true;
			printf( "Dumping a histogram.\n" );
		}
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	if( nomasks ) {
		fm	= NULL;
		fmd	= NULL;
	}
	else if( !fm && !fmd )
		nomasks = true;

// For backwards compatibility, also accept the env. variable

	const char	*pe = getenv( "FoldMaskThreshold" );

	if( pe ) {

		fmTOverride = atof( pe );

		printf(
		"Environment variable over-ride of threshold to %f\n",
		fmTOverride );

		printf(
		"--- This is obsolete ---  Use the -fmto option instead\n" );
	}
}

/* --------------------------------------------------------------- */
/* WriteFOLDMAP2Entry -------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFOLDMAP2Entry( int ncr )
{
	CMutex	M;
	char	name[256];

	sprintf( name, "fm_%d", gArgs.Z );

	if( M.Get( name ) ) {

		FILE *f = fopen( "fm.same", "a" );

		if( f ) {

			fprintf( f,
			"FOLDMAP2 %d.%d %d\n", gArgs.Z, gArgs.ID, ncr );
			fflush( f );
			fclose( f );
		}
	}

	M.Release();
}

/* --------------------------------------------------------------- */
/* SingleValueDominates ------------------------------------------ */
/* --------------------------------------------------------------- */

// Tells if a single pixel value dominates the histogram.
//
static bool SingleValueDominates( const vector<int> &histo )
{
	int	total = 0, nh = histo.size();

	for( int i = 0; i < nh; ++i )
		total += histo[i];

	for( int i = 0; i < nh; ++i ) {

		if( histo[i] > total/2 ) {

			printf(
			"Dominant value %d, %d occurences in %d pixels\n",
			i, histo[i], total );

			return true;
		}
	}

	return false;
}

/* --------------------------------------------------------------- */
/* IsTooGaussian ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Tells whether a histogram is too well explained by a gaussian.
//
static bool IsTooGaussian( const vector<int> &histo )
{
	int	nhist = histo.size();

// Mean and std deviation.

	MeanStd	mh;
	double	mean, std;
	int		ndv = 0;	// number of different values

	for( int i = 0; i < nhist; ++i ) {

		mh.Run( i, histo[i] );
		ndv += (histo[i] > 0);
	}

	mh.Stats( mean, std );

	printf( "\nMean %f, std %f, # different values %d\n",
		mean, std, ndv );

// Unusual case; print them

	if( ndv <= 20 ) {

		for( int i = 0; i < nhist; ++i ) {

			if( histo[i] > 0 )
				printf( "[%d %d] ", i, histo[i] );
		}

		printf( "\n" );
	}

// Now do a least squares fit. Should yield little diff, in theory.
// Create y array (same data as double); s array (const std dev).

	VecDoub	x( nhist ), y( nhist ), s( nhist );
	MeanStd orig;

	for( int i = 0; i < nhist; ++i ) {

		x[i] = i;
		y[i] = histo[i];
		s[i] = 1.0;
		orig.Element( histo[i] );
	}

	double	omean, ostd;

	orig.Stats( omean, ostd );

	printf( "Histogram values: mean %f rms %f\n", omean, ostd );

// Now try a 1 gaussian fit

	VecDoub	a( 3 );
	double	sum = 0.0;
	double	cnt = 0;
	int		ilo = max(   0, int(mean-std) ),
			ihi = min( 255, int(mean+std) );

	for( int i = ilo; i <= ihi; ++i ) {

		sum += histo[i];
		++cnt;

		printf( "%d ", histo[i] );
	}

	printf( "\n" );
	printf( "Boxcar average is %f\n", sum/cnt );

	a[0] = 1.5 * sum/cnt;	// estimated peak value
	a[1] = mean;			// where it is
	a[2] = std * sqrt( 2 );	// and the divisor term

	printf( "Before fit fit: height %f, loc %f, width %f\n",
		a[0], mean, std );

	Fitmrq f( x, y, s, a, fgauss );

	try {
		f.fit();
	}
	catch(int) {		// if the fitting procedure blows up
		return false;	// it's safe to assume it's not very gaussian
	}

	printf( "After fit: height %f, loc %f, width %f\n",
		f.a[0], f.a[1], f.a[2]/sqrt(2) );

// Now look at residuals

	MeanStd	m;

	for( int i = 0; i < nhist; ++i )
		m.Element( y[i] - ygauss( x[i], f.a ) );

	m.Stats( mean, std );

	printf( "Residuals: mean %f, RMS about mean %f\n", mean, std );

	return std < ostd/8;	// a guess
}

/* --------------------------------------------------------------- */
/* RemoveTooGaussian --------------------------------------------- */
/* --------------------------------------------------------------- */

// March a sampling window of size SxS across raster in steps of
// size D < S, so windows overlap. Within each window get a value
// histogram and remove all pixels (in that window) from goodp if
// they look homogeneous in intensity. Homogeneous regions either
// are dominated by one value, or have a gaussian-like intensity
// distribution (narrow).
//
static void RemoveTooGaussian(
	vector<uint8>	&goodp,
	const uint8		*raster,
	int				w,
	int				h )
{
	vector<int>	histo( 256 );

	const int	D = 128;
	const int	S = 256;

	int	ymax = (int)ceil( (double(h) - S) / D ),
		xmax = (int)ceil( (double(w) - S) / D );

	for( int iy = 0; iy <= ymax; ++iy ) {

		int	yhi = min( h, iy*D + S ),
			ylo = yhi - S;

		for( int ix = 0; ix <= xmax; ++ix ) {

			int	xhi = min( w, ix*D + S ),
				xlo = xhi - S;

			memset( &histo[0], 0, 256*sizeof(int) );

			for( int jy = ylo; jy < yhi; ++jy ) {

				for( int jx = xlo; jx < xhi; ++jx )
					++histo[raster[jx + w*jy]];
			}

			if( SingleValueDominates( histo ) ||
				IsTooGaussian( histo ) ) {

				for( int jy = ylo; jy < yhi; ++jy ) {

					for( int jx = xlo; jx < xhi; ++jx )
						goodp[jx + w*jy] = 0;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* AccumNonSaturated --------------------------------------------- */
/* --------------------------------------------------------------- */

static double AccumNonSaturated(
	MeanStd				&m,
	const vector<uint8>	&goodp,
	const uint8			*raster,
	int					np,
	int					SAT )
{
	double	fracReal;
	int		ngood = 0;

	m.Reset();

	for( int i = 0; i < np; ++i ) {

		if( goodp[i] ) {

			int	pix = raster[i];

			++ngood;

			if( pix >= SAT && pix <= 255-SAT )
				m.Element( pix );
		}
	}

	fracReal = (double)m.HowMany() / ngood;

	printf( "SAT=%d: %d non-sat, %f percent non-sat\n",
		SAT, m.HowMany(), fracReal * 100.0 );

	return fracReal;
}

/* --------------------------------------------------------------- */
/* StatsForNonSaturated ------------------------------------------ */
/* --------------------------------------------------------------- */

// Get stats for pixels within SAT of 0, 255.
//
static int StatsForNonSaturated(
	double				&mean,
	double				&std,
	const vector<uint8>	&goodp,
	const uint8			*raster,
	int					w,
	int					h )
{
	MeanStd	m;
	double	fracReal;
	int		np = w * h, SAT;

	fracReal = AccumNonSaturated( m, goodp, raster, np, SAT = 3 );

	if( fracReal < 0.9 ) {
		printf( "Saturated image!  Retrying with SAT = 1\n");
		AccumNonSaturated( m, goodp, raster, np, SAT = 1 );
	}

	m.Stats( mean, std );
	printf( "Non-sat pixels: mean = %f stddev = %f\n", mean, std );

	return SAT;
}

/* --------------------------------------------------------------- */
/* StatsForRealPixels -------------------------------------------- */
/* --------------------------------------------------------------- */

// Calculate mean, std for those pixels that look real.
//
// 'Unreal' data suffer one or more of these faults:
// - Value histogram dominated by one value (very narrow).
// - Value histogram described by single gaussian (too narrow).
// - Value too close to saturation (within SAT of 0 or 255).
//
static int StatsForRealPixels (
	double		&mean,
	double		&std,
	const uint8	*raster,
	int			w,
	int			h )
{
	vector<uint8>	goodp( w * h, 1 );

	RemoveTooGaussian( goodp, raster, w, h );

	return StatsForNonSaturated( mean, std, goodp, raster, w, h );
}

/* --------------------------------------------------------------- */
/* ZeroWhitePixels ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Zero white pixels in raster iff white (near SAT) value is
// outside useful range, that is, if brighter than 2.5 sigma
// from mean.
//
static void ZeroWhitePixels(
	uint8	*raster,
	int		np,
	double	mean,
	double	std,
	int		SAT )
{
	if( mean + 2.5 * std < 255 ) {

		printf( "Removing white pixels\n" );

		for( int i = 0; i < np; ++i ) {

			if( raster[i] > 255 - SAT )
				raster[i] = 0;
		}
	}
}

/* --------------------------------------------------------------- */
/* SelectThreshAndD ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Choose value thresh such that normalized values > -thresh are
// real rather than folds. That's the test used in compiling cr's.
// The stats calculation calls 'good' pixels in range [SAT,255-SAT].
// So the highest fold value is SAT-1.
//
// Also select radius D such that all  pixels within D of a fold
// are ignored.
//
static void SelectThreshAndD(
	double	&thresh,
	int		&D,
	double	mean,
	double	std,
	int		SAT )
{
	thresh	= 4.0;	// try 4 stddevs and shrink as needed
	D		= 10;

	if( (SAT-1 - mean) / std > -thresh ) {

		// 4 stddevs too wide...
		// Set 95% of exact range.

		thresh = 0.95 * (mean - (SAT-1)) / std;

		printf( "Thresh forced down to %g\n", thresh );

		if( thresh < 2.0 ) {

			// Black area too expansive -> fragmentation...
			// Set as wide as possible, out to value = 0.5

			thresh	= (mean - 0.5) / std;
			D		= 0;

			printf(
			"Desparate action: Thresh = %g excluding only v = 0...\n"
			"...And setting D = 0\n", thresh );
		}
	}

	if( gArgs.fmTOverride != 0 ) {
		thresh = gArgs.fmTOverride;
		printf( "Explicit override of thresh %g\n", thresh );
	}

	if( gArgs.D > -1 ) {
		D = gArgs.D;
		printf( "Explicit override of D %d\n", D );
	}
}

/* --------------------------------------------------------------- */
/* RemoveLowContrast --------------------------------------------- */
/* --------------------------------------------------------------- */

// March a sampling window of size SxS across image in steps of
// size D < S, so windows overlap. Within each window evaluate
// contrast and if low, force those v-pixels below threshold.
//
static void RemoveLowContrast(
	vector<double>	&v,
	int				w,
	int				h,
	double			std,
	double			thresh )
{
	const int	D = 32;
	const int	S = 128;

	vector<double>	window( S*S );
	vector<int>		remove;

// Evaluate

	int	ymax = (int)ceil( (double(h) - S) / D ),
		xmax = (int)ceil( (double(w) - S) / D );

	for( int iy = 0; iy <= ymax; ++iy ) {

		int	yhi = min( h, iy*D + S ),
			ylo = yhi - S;

		for( int ix = 0; ix <= xmax; ++ix ) {

			int	xhi = min( w, ix*D + S ),
				xlo = xhi - S;

			for( int jy = ylo; jy < yhi; ++jy ) {

				memcpy(
					&window[S*(jy-ylo)],
					&v[xlo + w*jy],
					S*sizeof(double) );
			}

			if( IsLowContrast( window, std ) ) {

				for( int jy = ylo; jy < yhi; ++jy ) {

					for( int jx = xlo; jx < xhi; ++jx )
						remove.push_back( jx + w*jy );
				}
			}
		}
	}

// Remove

	int	nr = remove.size();

	printf( "Low contrast pixels = %d\n", nr );

	if( nr ) {

		thresh = -thresh - 1.0;

		for( int i = 0; i < nr; ++i )
			v[remove[i]] = thresh;
	}
}

/* --------------------------------------------------------------- */
/* Widen --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Remove pixels within radius D of any fold pixel.
//
static void Widen(
	vector<double>	&v,
	int				w,
	int				h,
	double			thresh,
	int				D )
{
	int	nr, np = w * h;

	vector<int>	remove;

	thresh = -thresh;

	for( int i = 0; i < np; ++i ) {

		if( v[i] <= thresh )
			remove.push_back( i );
	}

	if( nr = remove.size() ) {

		thresh -= 1.0;

		for( int ir = 0; ir < nr; ++ir ) {

			int	i = remove[ir],
				y = i / w,
				x = i - w * y,
				ylo = max(   0, y - D ),
				yhi = min( h-1, y + D ),
				xlo = max(   0, x - D ),
				xhi = min( w-1, x + D );

			for( y = ylo; y <= yhi; ++y ) {

				for( x = xlo; x <= xhi; ++x )
					v[x + w*y] = thresh;
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* ImageToFoldMap ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Here we convert an input image to a map where 0 = on fold,
// 1 = region 1, 2 = region 2, etc.
//
// Return maximal region id.
//
static int ImageToFoldMap(
	uint8*			FoldMask,
	const PicBase	&pic,
	bool			remove_low_contrast = false )
{
	double	mean, std, thresh;
	uint8	*raster = pic.raster;
	int		w = pic.w,
			h = pic.h,
			npixels = w * h,
			SAT, D;

// Calc fold thresh and width extension D

	SAT = StatsForRealPixels( mean, std, raster, w, h );

	ZeroWhitePixels( raster, npixels, mean, std, SAT );

	SelectThreshAndD( thresh, D, mean, std, SAT );

// Make normalized images

	vector<double>	v(npixels);

	for( int i = 0; i < npixels; ++i )
		v[i] = (raster[i] - mean) / std;

	vector<double>	vorig = v;

// Remove low contrast

	if( remove_low_contrast )
		RemoveLowContrast( v, w, h, std, thresh );

// Widen folds

	Widen( v, w, h, thresh, D );

// Propagate connected regions

	vector<ConnRegion>	cr;

	for( int i = 0; i < npixels; ++i ) {

		if( v[i] > -thresh ) {

			ConnRegion	c;
			int			npts;

			npts = Propagate( c.pts, v, w, h, i,
					-thresh, -thresh - 1.0 );

			printf(
			"ImageToFoldMap: ConnRegion with %d pixels\n", npts );

			if( npts > gArgs.minarea )
				cr.push_back( c );
		}
	}

// Final accounting

   SetBoundsAndColors( cr, FoldMask, vorig, w, h, -thresh, D );

   return cr.size();
}

/* --------------------------------------------------------------- */
/* MakeDrawingMask ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Create optimized drawing mask.
//
// Notes:
// For alignment, we remove/ignore low contrast areas to avoid poor
// alignment results. However, for drawing, once the transforms are
// already determined, we want to assign as many image pixels (hence
// map pixels) as possible to some appropriate transform ID, that is,
// to a map pixel value.
//
// The drawing mask typically has more non-zero pixels in it than
// the alignment mask and that's sufficient to cause the raster
// scanner that finds connected regions to find them in a different
// order, and so give them different ID's. Moreover, the extra
// pixels in the drawing map may bridge regions that are separate
// in the alignment map giving the two maps different region counts
// AND different labels.
//
// Getting the labels adjusted correctly has an easy part and a
// hard part. First, for a given pixel coordinate, if the alignment
// map has a non-zero value there, that's the value we want in the
// drawing map--easy. But if there is a (non-zero) pixel in the
// drawing map that has no (non-zero) pixel in the alignment map,
// what ID shall we use? There are several possible procedures, but
// what we elect to do here is to examine the drawing map (before
// it is relabeled) and for every non-zero region in it, determine
// which non-zero alignment region it has the largest area overlap
// with. This gives a mapping from any drawing label to an alignment
// label. This is simple to implement and probably good enough since
// low contrast cases are few.
//
static uint8* MakeDrawingMask(
	PicBase			&pic,
	const uint8		*FoldMaskAlign,
	int				np )
{
// First we generate the drawing mask in pretty much the same way
// as the alignment mask, except that low contrast regions are
// not removed for the drawing version.

	uint8	*FoldMaskDraw = (uint8*)malloc( np );

	ImageToFoldMap( FoldMaskDraw, pic, false );

// Next create the mapping from labels in the drawing map to
// that label/region in the alignment map with which there is
// greatest overlap.

	vector<int>	map( 256, 0 );

// Get a mapping for each label value used in the drawing map.
// Since region labels are assigned consecutively, as soon as
// we advance to a label that isn't used we are done.

	int used = true;

	for( int drawlabel = 1; drawlabel < 256 && used; ++drawlabel ) {

		vector<int>	alnarea( 256, 0 );

		used = false;

		for( int i = 0; i < np; ++i ) {

			if( FoldMaskDraw[i] == drawlabel ) {
				++alnarea[FoldMaskAlign[i]];
				used = true;
			}
		}

		if( used ) {

			// Find largest overlapped area, if any.

			int	max_A = 0;
			int	max_j = 0;

			for( int j = 1; j < 256; ++j ) {

				if( alnarea[j] > max_A ) {
					max_A = alnarea[j];
					max_j = j;
				}
			}

			printf(
			"Mapping drawing region %d to alignment region %d\n",
			drawlabel, max_j );

			map[drawlabel] = max_j;
		}
	}

// Now relabel drawing map pixels

	for( int i = 0; i < np; ++i ) {

		if( FoldMaskAlign[i] )
			FoldMaskDraw[i] = FoldMaskAlign[i];
		else if( FoldMaskDraw[i] )
			FoldMaskDraw[i] = map[FoldMaskDraw[i]];
	}

	return FoldMaskDraw;
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char **argv )
{
/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* -------- */
/* Load src */
/* -------- */

	PicBase	p;

	p.LoadOriginal( gArgs.infile, stdout, gArgs.transpose );

/* ---------------- */
/* Save src as nmrc */
/* ---------------- */

	if( strstr( gArgs.infile, ".mrc" ) && gArgs.nmrc )
		Raster8ToPng8( gArgs.nmrc, p.original, p.w, p.h );

/* ------------ */
/* Create masks */
/* ------------ */

	if( gArgs.nomasks ) {
		WriteFOLDMAP2Entry( 1 );
		return 0;
	}

// Create two different fold masks - one optimized for alignment,
// with all low contrast regions removed, and one optimized for
// drawing, with regions as big as possible.

	int		ncr, np = p.w * p.h;
	uint8	*FoldMaskAlign = (uint8*)malloc( np );
	uint8	*FoldMaskDraw;

	if( gArgs.oneregion ) {

		// Efficiently handle the '-one' option. In this case
		// the drawing and alignment masks are identical.

		ncr				= 1;
		FoldMaskDraw	= FoldMaskAlign;
		memset( FoldMaskAlign, 1, np );
	}
	else {

		// Otherwise, calculate real foldmasks

		ncr = ImageToFoldMap( FoldMaskAlign, p, true );

		if( gArgs.fmd )
			FoldMaskDraw = MakeDrawingMask( p, FoldMaskAlign, np );
	}

/* ------------------------------- */
/* Write masks and FOLDMAP entries */
/* ------------------------------- */

	if( gArgs.fm )
		Raster8ToPng8( gArgs.fm, FoldMaskAlign, p.w, p.h );

	if( gArgs.fmd )
		Raster8ToPng8( gArgs.fmd, FoldMaskDraw, p.w, p.h );

	WriteFOLDMAP2Entry( ncr );

	VMStats( stdout );
	return 0;
}


