

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Geometry.h"
#include	"CTForm.h"

#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	minolap	0.025

// Special override for Davi: allow tiny overlap
//#define	minolap	0.0003

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;	// file name
	double	r;		// inlayer radius from center
	int		z;		// Z layer
	int		id;		// inlayer id
	TForm	tr;		// local to global
	TForm	inv;	// global to local

public:
	Picture()	{id = -1;};
};


class Pair {

public:
	int	a, b;

public:
	Pair( int _a, int _b )	{a = _a; b = _b;};
};

/* --------------------------------------------------------------- */
/* CArgs_scr ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_scr {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images.
	CRegexID	re_id;

public:
	char	*infile,
			*outdir;
	int		zmin,
			zmax;
	bool	Connect,		// just connect two layers
			Simple,
			NoFolds,
			NoDirs;

public:
	CArgs_scr()
	{
		infile		=
		outdir		= "NoSuch";	// prevent overwriting real dir
		zmin		= 0;
		zmax		= 32768;
		Connect		= false;
		Simple		= false;
		NoFolds		= false;
		NoDirs		= false;
	};

	void SetCmdLine( int argc, char* argv[] );

	int DecodeID( const char *name );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const vector<Picture>	*_vpsort;

static char			gtopdir[2048];
static CArgs_scr	gArgs;
static FILE*		flog	= NULL;
static uint32		gW		= 0,	// universal pic dims
					gH		= 0;
static int			gZMax	= 0;
static int			ismrc	= false;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_scr::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "scr.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Make scripts: %s ", atime );

// parse command line args

	char	*pat;

	re_id.Set( "/N" );

	if( argc < 2 ) {
		printf( "Usage: scr <source-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( IsArg( "-connect", argv[i] ) )
			Connect = true;
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( IsArg( "-simple", argv[i] ) )
			Simple = true;
		else if( IsArg( "-nf", argv[i] ) )
			NoFolds = true;
		else if( IsArg( "-nd", argv[i] ) )
			NoFolds = NoDirs = true;
		else if( GetArgStr( outdir, "-d", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n" );

	re_id.Compile( flog );

	fflush( flog );
}

/* --------------------------------------------------------------- */
/* DecodeID ------------------------------------------------------ */
/* --------------------------------------------------------------- */

int CArgs_scr::DecodeID( const char *name )
{
	const char	*s = strrchr( name, '/' );
	int			id;

	if( !s ) {
		fprintf( flog, "No '/' in [%s].\n", name );
		exit( 42 );
	}

	if( !re_id.Decode( id, ++s ) ) {
		printf( "No tile-id found in [%s].\n", s );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* ParseSimple --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseSimple( vector<Picture> &vp )
{
	FILE	*fp = FileOpenOrDie( gArgs.infile, "r", flog );

/* ---------- */
/* Scan lines */
/* ---------- */

	for( ;; ) {

		Picture	p;
		char	name[2048];
		int		x, y, z;

		/* ---------- */
		/* Get a line */
		/* ---------- */

		if( fscanf( fp, "%s %d %d %d", name, &x, &y, &z ) != 4 )
			break;

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		if( z > gZMax )
			gZMax = z;

		/* ---------------------- */
		/* Read actual dimensions */
		/* ---------------------- */

		if( !gW ) {

			uint8 *ras = Raster8FromAny( name, gW, gH, flog );

			if( !ras || !gW ) {
				fprintf( flog, "Error loading [%s].\n", name );
				exit( 42 );
			}

			RasterFree( ras );
		}

		/* ----------------- */
		/* Set picture entry */
		/* ----------------- */

		p.fname	= name;
		p.z		= z;
		p.id	= gArgs.DecodeID( name );
		p.tr.SetXY( x, y );
		InvertTrans( p.inv, p.tr );

		vp.push_back( p );
	}

/* ----- */
/* Close */
/* ----- */

	fclose( fp );
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

	//fprintf( flog, "Child element value %s.\n", layer->Value() );

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

		//fprintf( flog, "Got a <t2_layer>.\n" );

		const char	*sz = layer->Attribute( "z" );
		int			z	= int(atof(sz) + 0.5);

		//fprintf( flog, "z = %s.\n", sz );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		if( z > gZMax )
			gZMax = z;

		/* ------------------------------ */
		/* For each patch (tile) in layer */
		/* ------------------------------ */

		TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

		for( ; ptch; ptch = ptch->NextSiblingElement() ) {

			//fprintf( flog, "Got a <t2_patch>.\n" );

			Picture		p;
			const char	*name = ptch->Attribute( "file_path" );

			/* ---- */
			/* Dims */
			/* ---- */

			if( !gW ) {

				uint8 *ras = Raster8FromAny( name, gW, gH, flog );

				if( !ras || !gW ) {
					fprintf( flog, "Error loading [%s].\n", name );
					exit( 42 );
				}

				RasterFree( ras );
			}

			/* ----------------- */
			/* Set picture entry */
			/* ----------------- */

			p.fname	= name;
			p.z		= z;
			p.id	= gArgs.DecodeID( name );

			p.tr.ScanTrackEM2( ptch->Attribute( "transform" ) );
			InvertTrans( p.inv, p.tr );

			vp.push_back( p );
		}
	}
}

/* --------------------------------------------------------------- */
/* Sort_z_inc ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_z_inc( const Picture &A, const Picture &B )
{
	return A.z < B.z;
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
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Bounds(
	DBox					&B,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	const double	BIGD = 1.0e30;

	B.L = BIGD, B.R = -BIGD,
	B.B = BIGD, B.T = -BIGD;

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];
		vector<Point>	cnr( 4 );

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		P.tr.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			B.L = fmin( B.L, cnr[k].x );
			B.R = fmax( B.R, cnr[k].x );
			B.B = fmin( B.B, cnr[k].y );
			B.T = fmax( B.T, cnr[k].y );
		}
	}
}

/* --------------------------------------------------------------- */
/* AssignR ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void AssignR(
	vector<Picture>			&vp,
	int						is0,
	int						isN,
	const DBox				&B )
{
	Point	C( (B.L + B.R) / 2, (B.B + B.T) / 2 );

	for( int i = is0; i < isN; ++i ) {

		Picture&	P = vp[i];
		Point		cnr;	// (0,0)

		P.tr.Transform( cnr );
		P.r = cnr.DistSqr( C );
	}
}

/* --------------------------------------------------------------- */
/* Sort_r_inc ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_r_inc( const Picture &A, const Picture &B )
{
	return A.r < B.r || (A.r == B.r && A.id < B.id);
}

/* --------------------------------------------------------------- */
/* SortTilesRadially --------------------------------------------- */
/* --------------------------------------------------------------- */

// Within each layer/montage, sort the tiles into ascending
// order of distance from center. We want jobs to be issued
// from the center out to get more reliable matches earlier.
//
static void SortTilesRadially( vector<Picture> &vp )
{
// First, just sort into layers

	sort( vp.begin(), vp.end(), Sort_z_inc );

// For each layer...

	int		is0, isN;

	GetLayerLimits( vp, is0 = 0, isN );

	while( isN != -1 ) {

		// Get montage bounds
		// Assign radii from montage center
		// Sort this layer by r

		DBox	B;

		Bounds( B, vp, is0, isN );
		AssignR( vp, is0, isN, B );
		sort( vp.begin() + is0, vp.begin() + isN, Sort_r_inc );

		GetLayerLimits( vp, is0 = isN, isN );
	}
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
	char	name[2048];

// gtopdir gets the full path to the top directory
	getcwd( gtopdir, sizeof(gtopdir) );
	sprintf( name, "/%s", gArgs.outdir );
	strcat( gtopdir, name );

// create the top dir
	DskCreateDir( gArgs.outdir, flog );

// create stack subdir
	sprintf( name, "%s/stack", gArgs.outdir );
	DskCreateDir( name, flog );
}

/* --------------------------------------------------------------- */
/* WriteImageparamsFile ------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteImageparamsFile()
{
	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/imageparams.txt", gArgs.outdir );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "IMAGESIZE %d %d\n", gW, gH );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* CreateLayerDir ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Each layer gets a directory named by its z-index. All content
// pertains to this layer, or to this layer acting as a source
// onto itself or other layers.
//
// For example, make.down will contain ptest jobs aligning this
// layer onto that below (z-1).
//
static void CreateLayerDir( char *lyrdir, int L )
{
	fprintf( flog, "\n\nCreateLayerDir: layer %d\n", L );

	sprintf( lyrdir, "%s/%d", gArgs.outdir, L );
	DskCreateDir( lyrdir, flog );
}

/* --------------------------------------------------------------- */
/* CreateTileSubdirs --------------------------------------------- */
/* --------------------------------------------------------------- */

// Each tile gets a directory named by its picture id. All content
// pertains to this tile, or this tile acting as a source onto
// other tiles.
//
// For example, folder 8/10 contains the foldmask fm.png for tile
// 10 in layer 8. If this folder contains file 7.11.tf.txt it
// lists the transforms mapping tile 8/10 onto tile 7/11.
//
static void CreateTileSubdirs(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	fprintf( flog, "--CreateTileSubdirs: layer %d\n", vp[is0].z );

	for( int i = is0; i < isN; ++i ) {

		char	subdir[2048];

		sprintf( subdir, "%s/%d", lyrdir, vp[i].id );
		DskCreateDir( subdir, flog );
	}
}

/* --------------------------------------------------------------- */
/* Sort_tile_inc ------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool Sort_tile_inc( int a, int b )
{
	return (*_vpsort)[a].id < (*_vpsort)[b].id;
}

/* --------------------------------------------------------------- */
/* Make_TileToImage ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Write TileToImage.txt for this layer.
// Each row gives {tile, global-Tr, image path}.
//
static void Make_TileToImage(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
// Locally sort entries in this file by tile-id

	isN -= is0;

	vector<int>	order( isN );

	for( int i = 0; i < isN; ++i )
		order[i] = is0 + i;

	_vpsort = &vp;

	sort( order.begin(), order.end(), Sort_tile_inc );

// Open file

	char	name[2048];
	FILE	*f;

	sprintf( name, "%s/TileToImage.txt", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// Header

	fprintf( f, "Tile\tT0\tT1\tX\tT3\tT4\tY\tPath\n" );

// Write sorted entries
// Use given name, unless mrc images.

	if( !ismrc ) {

		for( int i = 0; i < isN; ++i ) {

			const Picture&	P = vp[order[i]];
			const double*	T = P.tr.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
				P.id, T[0], T[1], T[2], T[3], T[4], T[5],
				P.fname.c_str() );
		}
	}
	else {

		char	basepath[2048];

		sprintf( basepath, "%s/%d/", gtopdir, vp[order[0]].z );

		for( int i = 0; i < isN; ++i ) {

			const Picture&	P = vp[order[i]];
			const double*	T = P.tr.t;

			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s%d/nmrc_%d_%d.png\n",
				P.id, T[0], T[1], T[2], T[3], T[4], T[5],
				basepath, P.id, P.z, P.id );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_ThmPairFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void OneTprFile( const char *lyrdir, int a, int b )
{
    char	name[2048];
    FILE	*f;

// approximate transform file

	sprintf( name, "%s/ThmPair_%d_@_%d.txt", lyrdir, a, b );

	f = FileOpenOrDie( name, "w", flog );

	fprintf( f, "Atl\tBtl\tAcr\tBcr\tErr\tDeg\tQ\tR"
	"\tT0\tT1\tX\tT3\tT4\tY\n" );

	fclose( f );
}


// For a layer, make ThmPair files for this and adjacent layers.
//
static void Make_ThmPairFile(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						id0,
	int						iu0 )
{
	fprintf( flog, "--Make_ThmPairFile: layer %d\n", vp[is0].z );

	OneTprFile( lyrdir, vp[is0].z, vp[is0].z );

	if( id0 != -1 )
		OneTprFile( lyrdir, vp[is0].z, vp[id0].z );

	//if( iu0 != -1 )
	//	OneTprFile( lyrdir, vp[is0].z, vp[iu0].z );
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
	OrdVert( const vertex& _v )
		{v = _v;};

	bool operator <  (const OrdVert &rhs) const
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
static double ABOlap( const Picture &a, const Picture &b )
{
	vector<vertex>	va( 4 ), vb( 4 );

// Set b vertices (b-system) CCW ordered for LeftSide calcs

	vb[0] = vertex( 0   , 0 );
	vb[1] = vertex( gW-1, 0 );
	vb[2] = vertex( gW-1, gH-1 );
	vb[3] = vertex( 0   , gH-1 );

// Set a vertices (b-system) CCW ordered for LeftSide calcs

	{
		TForm			T;	// map a->b
		vector<Point>	p( 4 );

		MultiplyTrans( T, b.inv, a.tr );

		p[0] = Point( 0.0 , 0.0 );
		p[1] = Point( gW-1, 0.0 );
		p[2] = Point( gW-1, gH-1 );
		p[3] = Point( 0.0 , gH-1 );

		T.Transform( p );

		for( int i = 0; i < 4; ++i ) {

			va[i].x = int(p[i].x);
			va[i].y = int(p[i].y);

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

	double	A = AreaOfPolygon( pgon ) / (gW * gH);

	fprintf( flog, "----ABOlap: Tile %3d - %3d; area frac %f.\n",
	a.id, b.id, A );

	return A;
}

/* --------------------------------------------------------------- */
/* ConvertSpaces ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Make files have targets, dependencies and rules. Make can not
// scan dependency strings that contain embedded spaces. However,
// it works if " " is substituted by "\ ".
//
// Note too, that make does not like dependency strings to be
// enclosed in any quotes, so that will not solve this issue.
//
static void ConvertSpaces( char *out, const char *in )
{
	while( *out++ = *in++ ) {

		if( in[-1] == ' ' ) {

			out[-1]	= '\\';
			*out++	= ' ';
		}
	}
}

/* --------------------------------------------------------------- */
/* WriteThumbMakeFile -------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell thumbs to process the pairs
// of images described by (P). Argument mkname is a string from
// {"same", "up", "down"}.
//
static void WriteThumbMakeFile(
	const char				*lyrdir,
	const char				*mkname,
	const vector<Picture>	&vp,
	const vector<Pair>		&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/thumbs.%s", lyrdir, mkname );

	f = FileOpenOrDie( name, "w", flog );

// write 'all' targets line

	fprintf( f, "all:\n" );

// rule lines

	for( int i = 0; i < np; ++i ) {

		const Picture&	A = vp[P[i].a];
		const Picture&	B = vp[P[i].b];

		fprintf( f,
		"\tthumbs %d/%d@%d/%d ${EXTRA}\n",
		A.z, A.id, B.z, B.id );
	}

	fprintf( f, "\n" );

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_ThumbsSame ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting thumbs jobs for pairs
// of intersecting images within this layer.
//
// (a, b) = (source, target).
//
static void Make_ThumbsSame(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_ThumbsSame: layer %d\n", vp[is0].z );

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = a + 1; b < isN; ++b ) {

			if( ABOlap( vp[a], vp[b] ) > minolap )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

	WriteThumbMakeFile( lyrdir, "same", vp, P );
}

/* --------------------------------------------------------------- */
/* Make_ThumbsDown ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting thumbs jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[below layer]).
//
static void Make_ThumbsDown(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN,
	int						id0,
	int						idN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_ThumbsDown: layer %d @ %d\n",
		vp[is0].z, (id0 != -1 ? vp[id0].z : -1) );

// write dummy file even if no targets

	if( id0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = id0; b < idN; ++b ) {

			//if( vp[a].id != vp[b].id )
			//	continue;

			if( ABOlap( vp[a], vp[b] ) > minolap )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteThumbMakeFile( lyrdir, "down", vp, P );
}

/* --------------------------------------------------------------- */
/* Make_MakeFM --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write script to tell program tiny to calculate a foldmask for
// each tile.
//
static void Make_MakeFM(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	char	name[2048];
	FILE	*f;

	fprintf( flog, "--Make_MakeFM: layer %d\n", vp[is0].z );

	sprintf( name, "%s/make.fm", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// master target depends on all others

	fprintf( f, "all: " );

	for( int i = is0; i < isN; ++i )
		fprintf( f, "%d/fm.png ", vp[i].id );

	fprintf( f, "\n\n" );

// subtargets and rules

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];
		char dep[2048];

		ConvertSpaces( dep, P.fname.c_str() );

		fprintf( f, "%d/fm.png: %s\n", P.id, dep );

		if( gArgs.NoFolds ) {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=%d/nmrc_%d_%d.png'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str(),
				P.id, P.z, P.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str() );
			}
		}
		else {
			if( ismrc ) {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-nmrc=%d/nmrc_%d_%d.png'"
				" '-fm=%d/fm.png'"
				" '-fmd=%d/fmd.png'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str(),
				P.id, P.z, P.id,
				P.id,
				P.id );
			}
			else {
				fprintf( f,
				"\ttiny %d %d '%s'"
				" '-fm=%d/fm.png'"
				" '-fmd=%d/fmd.png'"
				" ${EXTRA}\n",
				P.z, P.id, P.fname.c_str(),
				P.id,
				P.id );
			}
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_fmsame --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write fm.same file with FOLDMAP2 entry for each tile.
//
// This is used only for '-nf' option because these entries
// all get connected-region count = 1.
//
static void Make_fmsame(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	char	name[2048];
	FILE	*f;

	fprintf( flog, "--Make_fmsame: layer %d\n", vp[is0].z );

	sprintf( name, "%s/fm.same", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// FOLDMAP2 entries

	for( int i = is0; i < isN; ++i ) {

		const Picture&	P = vp[i];

		fprintf( f, "FOLDMAP2 %d.%d 1\n", P.z, P.id );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteMakeFile ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Actually write the script to tell ptest to process the pairs
// of images described by (P). Argument mkname is a string from
// {"same", "up", "down"}.
//
static void WriteMakeFile(
	const char				*lyrdir,
	const char				*mkname,
	const vector<Picture>	&vp,
	const vector<Pair>		&P )
{
    char	name[2048];
	FILE	*f;
	int		np = P.size();

// open the file

	sprintf( name, "%s/make.%s", lyrdir, mkname );

	f = FileOpenOrDie( name, "w", flog );

// write 'all' targets line

	fprintf( f, "all: " );

	for( int i = 0; i < np; ++i ) {

		const Picture&	A = vp[P[i].a];
		const Picture&	B = vp[P[i].b];

		fprintf( f, "%d/%d.%d.map.tif ", A.id, B.z, B.id );
	}

	fprintf( f, "\n\n" );

// Write each 'target: dependencies' line
//		and each 'rule' line

	if( gArgs.NoFolds ) {

		for( int i = 0; i < np; ++i ) {

			const Picture&	A = vp[P[i].a];
			const Picture&	B = vp[P[i].b];

			fprintf( f,
			"%d/%d.%d.map.tif:\n",
			A.id, B.z, B.id );

			fprintf( f,
			"\tptest %d/%d@%d/%d -nf ${EXTRA}\n\n",
			A.z, A.id, B.z, B.id );
		}
	}
	else {

		for( int i = 0; i < np; ++i ) {

			const Picture&	A = vp[P[i].a];
			const Picture&	B = vp[P[i].b];
			char depA[2048], depB[2048];

			ConvertSpaces( depA, A.fname.c_str() );
			ConvertSpaces( depB, B.fname.c_str() );

			fprintf( f,
			"%d/%d.%d.map.tif:"
			" %s %s"
			" ../%d/%d/fm.png"
			" ../%d/%d/fm.png\n",
			A.id, B.z, B.id,
			depA, depB,
			A.z, A.id,
			B.z, B.id );

			fprintf( f,
			"\tptest %d/%d@%d/%d ${EXTRA}\n\n",
			A.z, A.id, B.z, B.id );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Make_MakeSame ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images within this layer.
//
// (a, b) = (source, target).
//
static void Make_MakeSame(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeSame: layer %d\n", vp[is0].z );

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = a + 1; b < isN; ++b ) {

			if( ABOlap( vp[a], vp[b] ) > minolap )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

	WriteMakeFile( lyrdir, "same", vp, P );
}

/* --------------------------------------------------------------- */
/* Make_MakeDown ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[below layer]).
//
static void Make_MakeDown(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN,
	int						id0,
	int						idN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeDown: layer %d @ %d\n",
		vp[is0].z, (id0 != -1 ? vp[id0].z : -1) );

// write dummy file even if no targets

	if( id0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = id0; b < idN; ++b ) {

			//if( vp[a].id != vp[b].id )
			//	continue;

			if( ABOlap( vp[a], vp[b] ) > minolap )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteMakeFile( lyrdir, "down", vp, P );
}

/* --------------------------------------------------------------- */
/* Make_MakeUp --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Write a make file submitting ptest jobs for pairs
// of intersecting images across layers.
//
// (a, b) = (source[this layer], target[above layer]).
//
static void Make_MakeUp(
	const char				*lyrdir,
	const vector<Picture>	&vp,
	int						is0,
	int						isN,
	int						iu0,
	int						iuN )
{
	vector<Pair>	P;

	fprintf( flog, "--Make_MakeUp: layer %d @ %d\n",
		vp[is0].z, (iu0 != -1 ? vp[iu0].z : -1) );

// write dummy file even if no targets

	if( iu0 == -1 )
		goto write;

// collect job indices

	for( int a = is0; a < isN; ++a ) {

		for( int b = iu0; b < iuN; ++b ) {

			if( vp[a].id != vp[b].id )
				continue;

			if( ABOlap( vp[a], vp[b] ) > minolap )
				P.push_back( Pair( a, b ) );
		}
	}

// write jobs

write:
	WriteMakeFile( lyrdir, "up", vp, P );
}

/* --------------------------------------------------------------- */
/* ForEachLayer -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Loop over layers, creating all: subdirs, scripts, work files.
//
static void ForEachLayer( const vector<Picture> &vp )
{
	int		id0, idN, is0, isN, iu0, iuN;

	id0 = -1;
	idN = -1;
	GetLayerLimits( vp, is0 = 0, isN );
	GetLayerLimits( vp, iu0 = isN, iuN );

	while( isN != -1 ) {

		char	lyrdir[2048];

		CreateLayerDir( lyrdir, vp[is0].z );

		if( !gArgs.NoDirs )
			CreateTileSubdirs( lyrdir, vp, is0, isN );

		Make_TileToImage( lyrdir, vp, is0, isN );

		Make_ThmPairFile( lyrdir, vp, is0, id0, iu0 );

		//Make_ThumbsSame( lyrdir, vp, is0, isN );
		//Make_ThumbsDown( lyrdir, vp, is0, isN, id0, idN );

		if( gArgs.NoFolds ) {

			if( ismrc )
				Make_MakeFM( lyrdir, vp, is0, isN );
			else
				Make_fmsame( lyrdir, vp, is0, isN );
		}
		else
			Make_MakeFM( lyrdir, vp, is0, isN );

		Make_MakeSame( lyrdir, vp, is0, isN );
		Make_MakeDown( lyrdir, vp, is0, isN, id0, idN );
		//Make_MakeUp( lyrdir, vp, is0, isN, iu0, iuN );

		id0 = is0;
		idN = isN;
		is0 = iu0;
		isN = iuN;
		GetLayerLimits( vp, iu0 = iuN, iuN );
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Picture>	vp;

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ---------------- */
/* Read source file */
/* ---------------- */

	if( gArgs.Simple )
		ParseSimple( vp );
	else
		ParseTrakEM2( vp );

	fprintf( flog, "Got %d images.\n", vp.size() );

	if( !vp.size() )
		goto exit;

	ismrc = strstr( vp[0].fname.c_str(), ".mrc" ) != NULL;

/* ------------------------------------------------- */
/* Within each layer, sort tiles by dist from center */
/* ------------------------------------------------- */

	SortTilesRadially( vp );

/* ------------------- */
/* Handle connect mode */
/* ------------------- */

// Connect mode expects exactly two adjacent layers
// and creates only make.up and make.down for them.

	if( gArgs.Connect ) {

		if( vp.size() < 2 || vp[0].z != vp[vp.size()-1].z - 1 ) {

			fprintf( flog,
			"Bogons! Expected two consecutive layers for"
			" connect mode.\n" );
			exit( 42 );
		}

		char	lyrdir[2048];
		int		is0, isN, iu0, iuN;

		// makeUp for layer 0 onto 1

		GetLayerLimits( vp, is0 = 0, isN );
		iu0 = isN;
		iuN = vp.size();

		//sprintf( lyrdir, "%s/0", gArgs.outdir );
		//Make_MakeUp( lyrdir, vp, is0, isN, iu0, iuN );

		// makeDown for layer 1 onto 0

		sprintf( lyrdir, "%s/1", gArgs.outdir );
		Make_MakeDown( lyrdir, vp, iu0, iuN, is0, isN );

		goto exit;
	}

/* --------------- */
/* Create dir tree */
/* --------------- */

	CreateTopDir();

	WriteImageparamsFile();

	ForEachLayer( vp );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


