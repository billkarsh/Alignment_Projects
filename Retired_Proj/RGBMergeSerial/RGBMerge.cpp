

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"ImageIO.h"
#include	"Maths.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;	// file name
	int		z;		// Z layer

public:
	Picture( const char *_n, int _z )	{fname=_n; z=_z;};

	bool operator < (const Picture &rhs) const
		{return z < rhs.z;};
};


class Fix {

public:
	int	z;

public:
	Fix( int _z )
		{z = _z;};
};

/* --------------------------------------------------------------- */
/* CArgs_rgbm ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_rgbm {

public:
	char	*infile;
	int		R, G, B,
			zmin, zmax;

public:
	CArgs_rgbm()
	{
		infile	= NULL;
		R		= -1;
		G		= -1;
		B		= -1;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_rgbm	gArgs;
static FILE*		flog = NULL;
static uint32		gW = 0,	gH = 0;		// universal pic dims






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_rgbm::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "RGBMerge.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 2 ) {
		printf( "Usage: RGBMerge <xml-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &R, "-R=%d", argv[i] ) )
			;
		else if( GetArg( &G, "-G=%d", argv[i] ) )
			;
		else if( GetArg( &B, "-B=%d", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
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

static void GetTiles(
	vector<Picture>	&vp,
	TiXmlElement*	layer,
	int				z )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	if( ptch && !gW ) {
		gW = atoi( ptch->Attribute( "width" ) );
		gH = atoi( ptch->Attribute( "height" ) );
	}

	for( ; ptch; ptch = ptch->NextSiblingElement() )
		vp.push_back( Picture( ptch->Attribute( "file_path" ), z ) );
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

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		GetTiles( vp, layer, z );
	}
}

/* --------------------------------------------------------------- */
/* GetLayerLimits ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Given starting index i0, which selects a layer z, set iN to be
// one beyond the highest index of a picture having same z. This
// makes loop limits [i0,iN) exclusive.
//
// If i0 or iN are out of bounds, both are set to -1.
//
static void GetLayerLimits(
	const vector<Picture>	&vp,
	int						&i0,
	int						&iN )
{
	int	np = vp.size();

	if( i0 < 0 || i0 >= np ) {

		i0 = -1;
		iN = -1;
	}
	else {

		int	Z = vp[i0].z;

		for( iN = i0 + 1; iN < np && vp[iN].z == Z; ++iN )
			;
	}
}

/* --------------------------------------------------------------- */
/* ChanName ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static char *ChanName( char *buf, const Picture &p, int chn )
{
	int	len = sprintf( buf, "%s", p.fname.c_str() );

	buf[len - 5] = '0' + chn;

	return buf;
}

/* --------------------------------------------------------------- */
/* RGBName ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static char *RGBName( char *buf, const Picture &p )
{
	char		name[256];
	const char	*p1, *p2;
	int			len;

// name part

	p1 = p.fname.c_str();
	p2 = strrchr( p1, '/' );
	len = sprintf( name, "%s", p2 + 1 );
	strcpy( name + len - 5, "RGB.tif" );

// full path

	sprintf( buf, "%.*s_RGB/%s", p2 - p1, p1, name );

	return buf;
}

/* --------------------------------------------------------------- */
/* GetRas16 ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static uint16* GetRas16( const Picture &p, int chn )
{
	char	buf[2048];
	uint32	w, h;

	return Raster16FromTif16( ChanName( buf, p, chn ), w, h, flog );
}

/* --------------------------------------------------------------- */
/* Scale1Lyr1Clr ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Scale1Lyr1Clr(
	int						&mn,
	int						&mx,
	int						chn,
	const vector<Picture>	&vp,
	int						i0,
	int						iN )
{
	const int		nbins = 65536;	// also max val
	vector<double>	bins( nbins, 0.0 );
	double			uflo = 0.0,
					oflo = 0.0;
	int				T,
					imin;

// histogram whole layer

	for( int i = i0; i < iN; ++i ) {

		uint16*	ras = GetRas16( vp[i], chn );

		Histogram( uflo, oflo, &bins[0], nbins,
			0.0, nbins, ras, gW * gH, false );

		RasterFree( ras );
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
// but exclude values lower than imin + 2 sdev.

	imin += T;
	T = int(OtsuThresh( &bins[imin], nbins - imin, imin, nbins ));

// Now, we aren't going to find objects with the Otsu threshold,
// however, we are going to set the scale maximum as 99% of the
// foreground counts.

	mx = PercentileBin( &bins[T], nbins - T, 0.99 );
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
	sprintf( buf, "%.*s_RGB", p2 - p1, p1 );

// make dir

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
	uint16*	ras = GetRas16( p, chn );
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

/* --------------------------------------------------------------- */
/* ScaleThisLayer ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ScaleThisLayer(
	vector<Fix>				&fix,
	const vector<Picture>	&vp,
	int						i0,
	int						iN )
{
	int	rmn, rmx, gmn, gmx, bmn, bmx;

// scan once to get scale for the layer/color

	if( gArgs.R >= 0 )
		Scale1Lyr1Clr( rmn, rmx, gArgs.R, vp, i0, iN );

	if( gArgs.G >= 0 )
		Scale1Lyr1Clr( gmn, gmx, gArgs.G, vp, i0, iN );

	if( gArgs.B >= 0 )
		Scale1Lyr1Clr( bmn, bmx, gArgs.B, vp, i0, iN );

// write each RGB image for this layer

	for( int i = i0; i < iN; ++i ) {

		vector<uint32>	RGB( gW * gH, 0xFF000000 );
		char			buf[2048];

		MakeFolder( vp[i] );

		if( gArgs.R >= 0 )
			AddChn( RGB, rmn, rmx, 0, gArgs.R, vp[i] );

		if( gArgs.G >= 0 )
			AddChn( RGB, gmn, gmx, 8, gArgs.G, vp[i] );

		if( gArgs.B >= 0 )
			AddChn( RGB, bmn, bmx, 16, gArgs.B, vp[i] );

		Raster32ToTifRGBA( RGBName( buf, vp[i] ), &RGB[0], gW, gH );
	}

	fix.push_back( Fix( vp[i0].z ) );
}

/* --------------------------------------------------------------- */
/* ScaleAllLayers ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Loop over layers, process images.
//
static void ScaleAllLayers(
	vector<Fix>				&fix,
	const vector<Picture>	&vp )
{
	int	i0, iN;

	GetLayerLimits( vp, i0 = 0, iN );

	while( iN != -1 ) {

		fprintf( flog, "\n---- Layer %d ----\n", vp[i0].z );

		ScaleThisLayer( fix, vp, i0, iN );

		printf( "Done layer %d\n", vp[i0].z );

		GetLayerLimits( vp, i0 = iN, iN );
	}
}

/* --------------------------------------------------------------- */
/* EditPath ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Change pattern
// from: .../dir/name_0.tif			// e.g. chan #0
// to:   .../dir_tag/name_tag.tif	// e.g. tag = RGB
//
static void EditPath( TiXmlElement* ptch )
{
	const char	*tag = "RGB";

	char		buf[2048], name[128];
	const char	*n = ptch->Attribute( "file_path" ),
				*s = strrchr( n, '/' ),
				*u = strrchr( s, '_' );

// get the '/name' part

	sprintf( name, "%.*s", u - s, s );

// rebuild path

	sprintf( buf, "%.*s_%s%s_%s.tif", s - n, n,
		tag, name, tag );

	ptch->SetAttribute( "file_path", buf );
}

/* --------------------------------------------------------------- */
/* FixTiles ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void FixTiles( TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		EditPath( ptch );
		ptch->SetAttribute( "type", 4 );
		ptch->SetAttribute( "min", 0 );
		ptch->SetAttribute( "max", 255 );
	}
}

/* --------------------------------------------------------------- */
/* WriteXML ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteXML( const vector<Fix> &fix )
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ---------- */
/* Fix layers */
/* ---------- */

	int	nz = fix.size();

	// for each layer we wish to fix...
	for( int iz = 0; iz < nz; ++iz ) {

		// advance xml to that layer
		while( atoi( layer->Attribute( "z" ) ) < fix[iz].z )
			layer = layer->NextSiblingElement();

		if( !layer )
			break;

		FixTiles( layer );
	}

/* ---- */
/* Save */
/* ---- */

	xml.Save( "xmltmp.txt", true );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Picture>	vp;
	vector<Fix>		fix;

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

/* ------------------ */
/* Sort by Z and tile */
/* ------------------ */

	sort( vp.begin(), vp.end() );

/* ----------------- */
/* Calculate scaling */
/* ----------------- */

	ScaleAllLayers( fix, vp );

/* ------------- */
/* Write new xml */
/* ------------- */

	WriteXML( fix );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


