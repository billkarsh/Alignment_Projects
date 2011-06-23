//
// Write entry in outfile GRScales.txt
// having scale min,max for given:
//
// z layer
// chn,
// pct,
// lrbt
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"CTForm.h"
#include	"Timer.h"

#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;
	TForm	T;

public:
	Picture( const TiXmlElement* ptch );
};

/* --------------------------------------------------------------- */
/* CArgs_gray ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_gray {

public:
	IBox	roi;
	double	pct;
	char	*infile;
	int		z, chn;

public:
	CArgs_gray()
	{
		roi.L	= roi.R = 0;
		pct		= 99.5;
		infile	= NULL;
		z		= 0;
		chn		= -1;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_gray	gArgs;
static FILE*		flog = NULL;
static uint32		gW = 0,	gH = 0;		// universal pic dims






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_gray::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "GraRan1Lyr.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: GraRan1Lyr <xml-file> -z=i [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		vector<int>	vi;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &z, "-z=%d", argv[i] ) )
			;
		else if( GetArg( &chn, "-chn=%d", argv[i] ) )
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
	if( gArgs.chn < 0 )
		fname = ptch->Attribute( "file_path" );
	else {

		char	buf[2048];
		int		len;

		len = sprintf( buf, "%s", ptch->Attribute( "file_path" ) );

		buf[len - 5] = '0' + gArgs.chn;

		fname = buf;
	}

	T.ScanTrackEM2( ptch->Attribute( "transform" ) );
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<Picture> &vp )
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );
	bool			loadOK = doc.LoadFile();

	if( !loadOK ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	layer = hDoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

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

		/* ------------------------------ */
		/* For each patch (tile) in layer */
		/* ------------------------------ */

		TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

		for( ; ptch; ptch = ptch->NextSiblingElement() ) {

			/* ---- */
			/* Dims */
			/* ---- */

			if( !gW ) {
				gW = atoi( ptch->Attribute( "width" ) );
				gH = atoi( ptch->Attribute( "height" ) );
			}

			vp.push_back( Picture( ptch ) );
		}

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
/* WriteScaleEntry ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScaleEntry( int smin, int smax )
{
	char	buf[2048];

	sprintf( buf, "GRTemp/z_%d.txt", gArgs.z );
	FILE *f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "%d\t%d\n", smin, smax );
	fclose( f );
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

	imin = IndexOfMaxVal( &bins[0], nbins );
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

	WriteScaleEntry( smin, smax );
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
	StopTiming( stdout, "GraRan1Lyr", T0 );

	return 0;
}


