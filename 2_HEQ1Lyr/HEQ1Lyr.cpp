//
// Write scaled 8-bit HEQ images to 'tag' folder using:
//
// tag
// z layer
// pct,
// lrbt
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

/* --------------------------------------------------------------- */
/* CArgs_heq ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_heq {

public:
	IBox	roi;
	double	pct;
	char	*infile,
			*tag;
	int		z;

public:
	CArgs_heq()
	{
		roi.L	= roi.R = 0;
		pct		= 99.5;
		infile	= NULL;
		tag		= NULL;
		z		= 0;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_heq	gArgs;
static FILE*		flog = NULL;
static uint32		gW = 0,	gH = 0;		// universal pic dims






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_heq::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "HEQ1Lyr.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf( "Usage: HEQ1Lyr <xml-file> <tag> -z=i [options].\n" );
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
		else if( GetArg( &pct, "-pct=%lf", argv[i] ) )
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
/* Picture::Picture ---------------------------------------------- */
/* --------------------------------------------------------------- */

Picture::Picture( const TiXmlElement* ptch )
{
	fname = ptch->Attribute( "file_path" );
	T.ScanTrackEM2( ptch->Attribute( "transform" ) );
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

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

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
/* MakeFolder ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeFolder( const Picture &p )
{
	char		buf[2048];
	const char	*p1, *p2;

// folder path

	p1 = p.fname.c_str();
	p2 = strrchr( p1, '/' );
	sprintf( buf, "%.*s_%s", p2 - p1, p1, gArgs.tag );

// make dir

	DskCreateDir( buf, flog );
}

/* --------------------------------------------------------------- */
/* OutName ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char *OutName( char *buf, const Picture &p )
{
	const char	*p1 = p.fname.c_str();
	const char	*p2 = strrchr( p1, '/' );
	const char	*p3 = FileDotPtr( p2 );

// full path

	sprintf( buf,
		"%.*s_%s"
		"%.*s.%s.tif",
		p2 - p1, p1, gArgs.tag,
		p3 - p2, p2, gArgs.tag );

	return buf;
}

/* --------------------------------------------------------------- */
/* WriteImages --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteImages(
	const vector<Picture>	&vp,
	int						smin,
	int						smax )
{
	double	scale	= 255.0 / log((double)smax - smin);
	int		np		= vp.size(),
			npx		= gW * gH;

	for( int i = 0; i < np; ++i ) {

		const Picture	&p = vp[i];
		char			buf[2048];

		if( !DskExists( p.fname.c_str() ) )
			continue;

		MakeFolder( p );

		uint32	w, h;
		uint16*	ras = Raster16FromTif16(
						p.fname.c_str(),
						w, h, flog );

		vector<uint8>	i8( npx );

		for( int i = 0; i < npx; ++i ) {

			if( ras[i] <= smin )
				i8[i] = 0;
			else {

				double	pix = log(ras[i] - (double)smin) * scale;

				if( pix > 255 )
					pix = 255;

				i8[i] = (uint8)pix;
			}
		}

		RasterFree( ras );

		Raster8ToTif8( OutName( buf, p ), &i8[0], gW, gH );
	}
}

/* --------------------------------------------------------------- */
/* ScaleLayer ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ScaleLayer( const vector<Picture> &vp )
{
	const int		nbins = 65536;	// also max val
	vector<double>	bins( nbins, 0.0 );
	double			uflo = 0.0,
					oflo = 0.0;
	int				np   = vp.size(),
					T, imin, smin, smax;

// histogram whole layer

	for( int i = 0; i < np; ++i ) {

		if( !InROI( vp[i] ) )
			continue;

		if( !DskExists( vp[i].fname.c_str() ) )
			continue;

		uint32	w, h;
		uint16*	ras = Raster16FromTif16(
						vp[i].fname.c_str(),
						w, h, flog );

		Histogram( uflo, oflo, &bins[0], nbins,
			0.0, nbins, ras, w * h, false );

		RasterFree( ras );
	}

// smin is between lowest val and 2 sdev below mode
// Omit highest bins to avoid detector saturation

	imin = IndexOfMaxVal( &bins[0], nbins - 100 );
	T	 = int(2.0 * sqrt( imin ));
	smax = imin - T;
	smin = FirstNonzero( &bins[0], nbins );

	if( smax < 0 )
		smax = 0;

	if( smin < 0 )
		smin = 0;

	smin = (smin + smax) / 2;

// Get threshold for bkg-frg segmentation using Otsu method,
// but exclude values lower than imin + 2 sdev or higher than
// 98% of the object intensities.

	imin += T;
	smax  = imin + PercentileBin( &bins[imin], nbins - imin, 0.98 );

	T = int(OtsuThresh( &bins[imin], smax - imin, imin, smax ));

// Now, we aren't going to find objects with the Otsu threshold,
// however, we are going to set the scale maximum as pct% of the
// foreground counts.

	if( gArgs.pct < 100.0 ) {

		smax = T +
		PercentileBin( &bins[T], nbins - T, gArgs.pct/100.0 );
	}
	else
		smax = T + int((nbins - T) * gArgs.pct/100.0);

	WriteImages( vp, smin, smax );
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
	StopTiming( stdout, "HEQ1Lyr", T0 );

	return 0;
}


