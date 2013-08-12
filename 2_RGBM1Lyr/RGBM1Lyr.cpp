//
// For 1 layer, create new folder TIF_tag, and
// create scaled RGB images there according to:
//
// tag
// z layer
// -R=chn,pct
// -G=chn,pct
// -B=chn,pct
// -spanRGB=LLT
// -lrbt
//
// -spanRGB option needs three character string. Each character is
// either {L=whole layer, T=per tile} setting span of tiles used
// to scale each channel.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"TAffine.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;
	TAffine	T;

public:
	Picture( const TiXmlElement* ptch );
};


Picture::Picture( const TiXmlElement* ptch )
{
	fname = ptch->Attribute( "file_path" );

	T.ScanTrackEM2( ptch->Attribute( "transform" ) );
}

/* --------------------------------------------------------------- */
/* CArgs_rgbm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_rgbm {

public:
	IBox		roi;
	double		pct[3];
	const char	*infile,
				*tag,
				*span;
	int			z, RGB[3];

public:
	CArgs_rgbm()
	{
		roi.L	= roi.R = 0;
		pct[0]	= 99.5;
		pct[1]	= 99.5;
		pct[2]	= 99.5;
		infile	= NULL;
		tag		= NULL;
		span	= "LLL";
		z		= 0;
		RGB[0]	= -1;
		RGB[1]	= -1;
		RGB[2]	= -1;
	};

	bool ScanChan( int chn, const char *pat, char *argv );
	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_rgbm	gArgs;
static FILE*		flog = NULL;
static uint32		gW = 0,	gH = 0;		// universal pic dims






/* --------------------------------------------------------------- */
/* ScanChan ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Scan channel/pct arguments like:
// -R=0[,95.5]
//
bool CArgs_rgbm::ScanChan( int chn, const char *pat, char *argv )
{
	int	c;

	if( 1 == sscanf( argv, pat, &c ) ) {

		double	p;

		RGB[chn] = c;

		if( argv[4] == ',' && 1 == sscanf( argv + 5, "%lf", &p ) )
			pct[chn] = p;

		return true;
	}

	return false;
}

/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_rgbm::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "RGBM1Lyr.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf(
		"Usage: RGBM1Lyr <xml-file> <tag>"
		" -z=i <-[R,G,B]=i,pct> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		vector<int>	vi;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile )
				infile = argv[i];
			else
				tag = argv[i];
		}
		else if( GetArg( &z, "-z=%d", argv[i] ) )
			;
		else if( ScanChan( 0, "-R=%d", argv[i] ) )
			;
		else if( ScanChan( 1, "-G=%d", argv[i] ) )
			;
		else if( ScanChan( 2, "-B=%d", argv[i] ) )
			;
		else if( GetArgStr( span, "-spanRGB=", argv[i] )
			&& (span[0] == 'L' || span[0] == 'T')
			&& (span[1] == 'L' || span[1] == 'T')
			&& (span[2] == 'L' || span[2] == 'T') )
			;
		else if( GetArgList( vi, "-lrbt=", argv[i] ) && vi.size() == 4 )
			memcpy( &roi, &vi[0], 4*sizeof(int) );
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* GetTiles ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetTiles( vector<Picture> &vp, TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	if( ptch && !gW ) {
		gW = atoi( ptch->Attribute( "width" ) );
		gH = atoi( ptch->Attribute( "height" ) );
	}

	for( ; ptch; ptch = ptch->NextSiblingElement() )
		vp.push_back( Picture( ptch ) );
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<Picture> &vp )
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z < gArgs.z )
			continue;

		GetTiles( vp, layer );
		break;
	}
}

/* --------------------------------------------------------------- */
/* InROI --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool InROI( const Picture &p )
{
	if( gArgs.roi.L == gArgs.roi.R )
		return true;

	const double xolap = gW * 0.8;
	const double yolap = gH * 0.8;

	Point	c1, c2( gW, gH );
	double	t;

	p.T.Transform( c1 );
	p.T.Transform( c2 );

	if( c2.x < c1.x ) {
		t    = c1.x;
		c1.x = c2.x;
		c2.x = t;
	}

	if( c2.y < c1.y ) {
		t    = c1.y;
		c1.y = c2.y;
		c2.y = t;
	}

	if( c2.x < gArgs.roi.L + xolap )
		return false;

	if( c1.x > gArgs.roi.R - xolap )
		return false;

	if( c2.y < gArgs.roi.B + yolap )
		return false;

	if( c1.y > gArgs.roi.T - yolap )
		return false;

	return true;
}

/* --------------------------------------------------------------- */
/* ChanName ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static char *ChanName( char *buf, const Picture &p, int chn )
{
	const char	*src = p.fname.c_str();
	const char	*s, *u;
	char		*plt;

	sprintf( buf, "%s", src );
	plt	= strstr( buf, "Plate1_0" );
	s	= strrchr( src, '/' );
	u	= strrchr( src, '_' );

	sprintf( plt + 8, "/TIF%.*s_%c.tif",
		u - s, s, '0' + chn );

	return buf;
}

/* --------------------------------------------------------------- */
/* RGBName ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char *RGBName( char *buf, const Picture &p )
{
	const char	*src = p.fname.c_str();
	const char	*s, *u;
	char		*plt;

	sprintf( buf, "%s", src );
	plt	= strstr( buf, "Plate1_0" );
	s	= strrchr( src, '/' );
	u	= strrchr( src, '_' );

	sprintf( plt + 8, "/TIF_%s%.*s_%s.tif",
		gArgs.tag, u - s, s, gArgs.tag );

	return buf;
}

/* --------------------------------------------------------------- */
/* GetRas16 ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool GetRas16( uint16* &ras, const Picture &p, int chn )
{
	char	buf[2048];
	uint32	w, h;

	if( !DskExists( ChanName( buf, p, chn ) ) )
		return false;

	ras = Raster16FromTif16( buf, w, h, flog );

	return true;
}

/* --------------------------------------------------------------- */
/* Scale1Clr ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Scale1Clr(
	int						&mn,
	int						&mx,
	int						rgb,
	int						ip,		// -1 = whole layer
	const vector<Picture>	&vp )
{
	const int		nbins = 65536;	// also max val
	vector<double>	bins( nbins, 0.0 );
	double			uflo = 0.0,
					oflo = 0.0;
	int				np   = vp.size(),
					T, imin;

// collect histogram

	if( ip >= 0 ) {

		// tile ip

		uint16*	ras;

		if( !GetRas16( ras, vp[ip], gArgs.RGB[rgb] ) ) {

			mn	= 0;
			mx	= 65536;
			return;
		}

		Histogram( uflo, oflo, &bins[0], nbins,
			0.0, nbins, ras, gW * gH, false );

		RasterFree( ras );

	}
	else {

		// whole layer

		for( int i = 0; i < np; ++i ) {

			uint16*	ras;

			if( !InROI( vp[i] ) )
				continue;

			if( !GetRas16( ras, vp[i], gArgs.RGB[rgb] ) )
				continue;

			Histogram( uflo, oflo, &bins[0], nbins,
				0.0, nbins, ras, gW * gH, false );

			RasterFree( ras );
		}
	}

// mn is between lowest val and 2 sdev below mode
// Omit highest bins to avoid detector saturation

	imin = IndexOfMaxVal( &bins[0], nbins - 100 );
	T	 = int(2.0 * sqrt( imin ));
	mx   = imin - T;
	mn   = FirstNonzero( &bins[0], nbins );

	if( mx < 0 )
		mx = 0;

	if( mn < 0 )
		mn = 0;

	mn = (mn + mx) / 2;

// Get threshold for bkg-frg segmentation using Otsu method,
// but exclude values lower than imin + 2 sdev or higher than
// 98% of the object intensities.

	imin += T;
	mx = imin + PercentileBin( &bins[imin], nbins - imin, 0.98 );

	T = int(OtsuThresh( &bins[imin], mx - imin, imin, mx ));

// Now, we aren't going to find objects with the Otsu threshold,
// however, we are going to set the scale maximum as pct% of the
// foreground counts.

	if( gArgs.pct[rgb] < 100.0 ) {

		mx = T +
		PercentileBin( &bins[T], nbins - T, gArgs.pct[rgb]/100.0 );
	}
	else
		mx = T + int((nbins - T) * gArgs.pct[rgb]/100.0);
}

/* --------------------------------------------------------------- */
/* MakeFolder ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeFolder( const Picture &p )
{
	char	buf[2048];

	sprintf( buf, "%s", p.fname.c_str() );
	sprintf( strstr( buf, "Plate1_0" ) + 8, "/TIF_%s", gArgs.tag );
	DskCreateDir( buf, flog );
}

/* --------------------------------------------------------------- */
/* AddChn -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void AddChn(
	vector<uint32>	&RGB,
	int				mn,
	int				mx,
	int				shf,
	int				chn,
	const Picture	&p )
{
	uint16*	ras;

	if( GetRas16( ras, p, chn ) ) {

		int		pix, np = gW * gH;

		for( int i = 0; i < np; ++i ) {

			if( (pix = ras[i]) <= mn )
				pix = 0;
			else if( pix >= mx )
				pix = 255;
			else
				pix = (pix - mn)*255/(mx - mn);

			RGB[i] |= (pix << shf);
		}

		RasterFree( ras );
	}
}

/* --------------------------------------------------------------- */
/* ScaleLayer ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScaleLayer( const vector<Picture> &vp )
{
	int	rmn, rmx, gmn, gmx, bmn, bmx, np = vp.size();

// scan once to get scale for the layer/color

	if( gArgs.RGB[0] >= 0 && gArgs.span[0] == 'L' )
		Scale1Clr( rmn, rmx, 0, -1, vp );

	if( gArgs.RGB[1] >= 0 && gArgs.span[1] == 'L' )
		Scale1Clr( gmn, gmx, 1, -1, vp );

	if( gArgs.RGB[2] >= 0 && gArgs.span[2] == 'L' )
		Scale1Clr( bmn, bmx, 2, -1, vp );

// write each RGB image for this layer

	for( int i = 0; i < np; ++i ) {

		const Picture	&p = vp[i];
		vector<uint32>	RGB( gW * gH, 0xFF000000 );
		char			buf[2048];

		MakeFolder( p );

		if( gArgs.RGB[0] >= 0 ) {

			if( gArgs.span[0] != 'L' )
				Scale1Clr( rmn, rmx, 0, i, vp );

			AddChn( RGB, rmn, rmx, 0, gArgs.RGB[0], p );
		}

		if( gArgs.RGB[1] >= 0 ) {

			if( gArgs.span[1] != 'L' )
				Scale1Clr( gmn, gmx, 1, i, vp );

			AddChn( RGB, gmn, gmx, 8, gArgs.RGB[1], p );
		}

		if( gArgs.RGB[2] >= 0 ) {

			if( gArgs.span[2] != 'L' )
				Scale1Clr( bmn, bmx, 2, i, vp );

			AddChn( RGB, bmn, bmx, 16, gArgs.RGB[2], p );
		}

		Raster32ToTifRGBA( RGBName( buf, p ), &RGB[0], gW, gH );
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Picture>	vp;
	clock_t			T0 = StartTiming();

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ---------------- */
/* Read source file */
/* ---------------- */

	ParseTrakEM2( vp );

	fprintf( flog, "Got %d images.\n", vp.size() );

	if( !vp.size() )
		goto exit;

/* ----------------- */
/* Calculate scaling */
/* ----------------- */

	ScaleLayer( vp );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );
	StopTiming( stdout, "RGBM1Lyr", T0 );

	return 0;
}


