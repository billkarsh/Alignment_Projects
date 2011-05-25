

// Based on Lou's 11/23/2010 copy of lsq.cpp


#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"LinEqu.h"
#include	"TrakEM2_DTD.h"
#include	"CPoint.h"
#include	"CTForm.h"

#include	<math.h>

#include	<map>
#include	<set>
#include	<stack>
#include	<string>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// ----------------------------------------

// Pair of matched connRgn

class CRPair {

public:
	int	a, b;   // always stored so a <= b

public:
	CRPair( int aa, int bb )
		{
			if( aa <= bb )
				{a = aa; b = bb;}
			else
				{a = bb; b = aa;}
		};

	bool operator < (const CRPair &rhs) const
		{return a < rhs.a || (a == rhs.a && b < rhs.b);};

	bool operator == (const CRPair &rhs) const
		{return a == rhs.a && b == rhs.b;};
};

// ----------------------------------------

// Any matched point pair (POINT entry)

class Constraint {

public:
	Point	p1, p2;	// the two points
	int		r1, r2;	// indexes to the two regions
	bool	used,	// is this constraint used?
			inlier;	// is this a RANSAC inlier or outlier?

public:
	Constraint( int rr1, Point pp1, int rr2, Point pp2 )
		{r1 = rr1; p1 = pp1; r2 = rr2; p2 = pp2;};
};

// ----------------------------------------

// Image identifier for FOLDMAP lookups

class ZID {

public:
	int	z, id;

public:
	ZID( int zz, int iid )
		{z = zz; id = iid;};

	bool operator < (const ZID &rhs) const
		{return z < rhs.z || (z == rhs.z && id < rhs.id);};

	bool operator == (const ZID &rhs) const
		{return z == rhs.z && id == rhs.id;};
};

/* --------------------------------------------------------------- */
/* DIR ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Maps name-tag to z-layer

// @@@ -----------------------------------------
// Possible for POINT entries to directly carry {z,tile,rgn}
// and get rid of DIR...at least here. mos does not use POINT
// entries.
// @@@ -----------------------------------------

class DIR {

private:
	typedef struct DirEntry {
		char	*tag;	// directory part of image name
		int		z;		// layer
		DirEntry( char *line )
			{
				z	= atoi( strtok( line + 3, " " ) );
				tag	= strdup( strtok( NULL, " \n" ) );
			};
	} DirEntry;

private:
	vector<DirEntry>	dirTbl;

public:
	void ReadDIRFile();
	int  ZFromName( const char *name ) const;
};

/* --------------------------------------------------------------- */
/* RGN ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Connected region

class RGN {

public:
	string	name;	// image filename
	int		z,		// layer
			id,		// index in layer
			rgn,	// connRgn index
			itr;	// global-space transform index

public:
	RGN( int zz, int iid );
	RGN( const char *c, const DIR *dir );
	RGN( const char *path, const char *key );
};

/* --------------------------------------------------------------- */
/* ZIDR ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Image identifier for POINT lookups

class ZIDR {

public:
	int	z, id, rgn;

public:
	ZIDR( const RGN& R )
		{z = R.z; id = R.id; rgn = R.rgn;};

	bool operator < (const ZIDR &rhs) const
		{
			if( z < rhs.z )
				return true;
			if( z > rhs.z )
				return false;
			if( id < rhs.id )
				return true;
			if( id > rhs.id )
				return false;

			return rgn < rhs.rgn;
		};

	bool operator == (const ZIDR &rhs) const
		{return z == rhs.z && id == rhs.id && rgn == rhs.rgn;};
};

/* --------------------------------------------------------------- */
/* zsort --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Orders regions by z, for writing files

class zsort{

public:
	int	z, id, rgn, i;

public:
	zsort()	{};

	zsort( const RGN& R, int ii )
		{z = R.z; id = R.id; rgn = R.rgn; i = ii;};

	bool operator < (const zsort &rhs) const
		{
			if( z < rhs.z )
				return true;
			if( z > rhs.z )
				return false;
			if( id < rhs.id )
				return true;
			if( id > rhs.id )
				return false;

			return rgn < rhs.rgn;
		};
};

/* --------------------------------------------------------------- */
/* CArgs_lsq ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_lsq {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images.
	//
	CRegexID	re_id;

public:
	// xml_type values: these are ImagePlus codes:
	// GRAY8		= 0
	// GRAY16		= 1
	// GRAY32		= 2
	// COLOR_256	= 3
	// COLOR_RGB	= 4
	//
	vector<RGN>	include_only;		// if given, include only these
	double		same_strength,
				square_strength,
				scale_strength,
				thresh,				// outlier if worse than this
				trim;				// trim this off XML images
	char		*pts_file,
				*dir_file,
				*tfm_file;
	int			unite_layer,
				ref_layer,
				max_pass,
				xml_type;
	bool		make_layer_square,
				use_all;			// align even if #pts < 3/tile

public:
	CArgs_lsq()
	{
		same_strength		= 10.0;
		square_strength		= 0.1;
		scale_strength		= 1.0;
		thresh				= 700.0;
		trim				= 0.0;
		pts_file			= NULL;
		dir_file			= NULL;
		tfm_file			= NULL;
		unite_layer			= -1;
		ref_layer			= -1;
		max_pass			= 1;
		xml_type			= 0;
		make_layer_square	= false;
		use_all				= false;
	};

	void SetCmdLine( int argc, char* argv[] );

	int TileIDFromName( const char *name );
};

/* --------------------------------------------------------------- */
/* RGD ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Paired points used for rigid alignment testing

class RGD {

private:
	typedef struct {
		// rgn-pair's pts
		vector<Point>	pa;
		vector<Point>	pb;
	} CorrPts;

private:
	vector<CorrPts>	r12Pts;

private:
	double CanAlign(
		const vector<Point>	&p1,
		const vector<Point>	&p2,
		bool				print );

public:
	void AddPOINTPair(
		int			r1,
		const Point	&p1,
		int			r2,
		const Point	&p2 );

	void TestPairAlignments();
};

/* --------------------------------------------------------------- */
/* CNX ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Connection table operations

class CNX {

private:
	static const int	cnxMinLinks = 3;

	typedef struct {
		vector<int>	linkto;	// which regions it connects to
		vector<int>	nlinks;	// # points in that connection
		int			seed;	// starting entry of conn graph
	} CnxEntry;

private:
	vector<CnxEntry>	cnxTbl;

private:
	void ListWeakConnections( set<CRPair> &r12Bad );
	void MaxConnectedSet( set<int> &ignore );
	void Set_itr_set_used( set<CRPair> &r12Bad, set<int> &ignore );

public:
	void AddCorrespondence( int r1, int r2 );
	void SelectIncludedImages();
};

/* --------------------------------------------------------------- */
/* EVL ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Error evaluation

class EVL {

private:
	typedef struct Error {

		double	amt;
		int		idx;	// which constraint

		Error()	{};
		Error( double err, int id )
			{amt = err; idx = id;};

		bool operator < (const Error &rhs) const
			{return amt < rhs.amt;};
	} Error;

	typedef struct SecErr {

		Point	loc;	// global error location
		double	err;	// worst err
		int		idx;	// which constraint

		SecErr()
			{err = 0; idx = -1;};

		SecErr( const Point	&p1,
				const Point	&p2,
				double		e,
				int			id )
			{
				loc = Point( (p1.x+p2.x)/2, (p1.y+p2.y)/2 );
				err = e;
				idx = id;
			};
	} SecErr;

private:
	vector<Error>	Epnt;	// each constraint
	vector<SecErr>	Ein,	// worst in- and between-layer
					Ebt;

private:
	void Tabulate(
		const vector<zsort>		&z,
		const vector<double>	&X );

	void Line(
		FILE	*f,
		double	xfrom,
		double	yfrom,
		double	xto,
		double	yto );

	void BoxOrCross( FILE *f, double x, double y, bool box );
	void Arrow( FILE *f, const Point &g1, const Point &g2 );

	void Print_be_and_se_files( const vector<zsort> &z );
	void Print_errs_by_layer( const vector<zsort> &z );

public:
	void Evaluate(
		const vector<zsort>		&z,
		const vector<double>	&X );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_lsq			gArgs;
static FILE					*FOUT;		// outfile: 'simple'
static vector<RGN>			vRgn;		// the regions
static map<CRPair, int>		r12Idx;		// idx from region-pair
static vector<Constraint>	vAllC;		// all POINT entries
static map<ZID, int>		nConRgn;	// # corr-rgn this image
static int					gW = 4056,	// default image dims
							gH = 4056,
							gNTr = 0;	// Set by Set_itr_set_used






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_lsq::SetCmdLine( int argc, char* argv[] )
{
// Parse command line args

	printf( "---- Read params ----\n" );

	vector<int>	vi;
	char		*pat;

	re_id.Set( "/N" );

	if( argc < 2 ) {
		printf(
		"Usage: lsq <file of points>"
		" [<file of directory-name - layer-number correspondences>]"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		char	*unite;

		if( argv[i][0] != '-' ) {

			if( !pts_file )
				pts_file = argv[i];
			else
				dir_file = argv[i];
		}
		else if( GetArgList( vi, "-only=", argv[i] ) ) {

			int ni = vi.size();

			if( ni & 1 )
				--ni;

			for( int i = 0; i < ni; i += 2 ) {

				include_only.push_back(
					RGN( vi[i], vi[i+1] ) );

				printf( "Include only %4d %4d.\n", vi[i], vi[i+1] );
			}
		}
		else if( GetArg( &same_strength, "-same=%lf", argv[i] ) )
			printf( "Setting same-layer strength to %f.\n", same_strength );
		else if( GetArg( &square_strength, "-square=%lf", argv[i] ) )
			printf( "Setting square strength to %f.\n", square_strength );
		else if( GetArg( &scale_strength, "-scale=%lf", argv[i] ) )
			printf( "Setting scale strength to %f.\n", scale_strength );
		else if( GetArg( &thresh, "-threshold=%lf", argv[i] ) )
			printf( "Setting threshold to %f.\n", thresh );
		else if( GetArg( &trim, "-trim=%lf", argv[i] ) )
			printf( "Setting trim amount to %f.\n", trim );
		else if( GetArgStr( unite, "-unite=", argv[i] ) ) {

			char	buf[2048];

			sscanf( unite, "%d,%s", &unite_layer, buf );
			tfm_file = strdup( buf );

			printf( "Uniting: layer %d of '%s'.\n",
			unite_layer, tfm_file );
		}
		else if( GetArg( &ref_layer, "-refz=%d", argv[i] ) )
			printf( "Reference layer %d.\n", ref_layer );
		else if( GetArg( &max_pass, "-pass=%d", argv[i] ) )
			printf( "Setting maximum passes to %d.\n", max_pass );
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			printf( "Setting xml image type to %d.\n", xml_type );
		else if( IsArg( "-mls", argv[i] ) ) {
			make_layer_square = true;
			printf( "Making reference layer square.\n" );
		}
		else if( IsArg( "-all", argv[i] ) ) {
			use_all = true;
			printf( "Using all correspondences.\n" );
		}
		else if( GetArgStr( pat, "-p", argv[i] ) ) {
			re_id.Set( pat );
			printf( "Setting pattern '%s'.\n", pat );
		}
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	re_id.Compile( stdout );
}

/* --------------------------------------------------------------- */
/* TileIDFromName ------------------------------------------------ */
/* --------------------------------------------------------------- */

int CArgs_lsq::TileIDFromName( const char *name )
{
	const char	*s = strrchr( name, '/' );
	int			id;

	if( !s ) {
		printf( "No '/' in '%s'.\n", name );
		exit( 42 );
	}

	if( !re_id.Decode( id, ++s ) ) {
		printf( "No tile-id found in '%s'.\n", s );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* ReadDIRFile --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read file who's entries associate:
//
//	[layer-#] <-> [imgname-tag].
//
// Example entry:
//
//	DIR 0 A1_
//
void DIR::ReadDIRFile()
{
	printf( "---- Read DIR ----\n" );

// Assume there is only layer z=0 if no file

	if( !gArgs.dir_file )
		return;

// Read file

	FILE	*f = FileOpenOrDie( gArgs.dir_file, "r" );

	for(;;) {

		char	line[2048];

		if( !fgets( line, sizeof(line), f ) )
			break;

		if( !strncmp( line, "DIR", 3 ) ) {
			fprintf( FOUT, line );
			dirTbl.push_back( DirEntry( line ) );
		}
		else {
			printf( "Bad line '%s' in DIR file.\n", line );
			exit( 42 );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ZFromName ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Find layer number by searching DIR entries.
// Inefficient, but should be OK here.
//
int DIR::ZFromName( const char *name ) const
{
	int	nd = dirTbl.size();

	for( int i = 0; i < nd; ++i ) {

		if( strstr( name, dirTbl[i].tag ) )
			return dirTbl[i].z;
	}

	return 0;	// not found, return layer 0
}

/* --------------------------------------------------------------- */
/* RGN::RGN ------------------------------------------------------ */
/* --------------------------------------------------------------- */

RGN::RGN( int zz, int iid )
{
	name	= "";
	z		= zz;
	id		= iid;
	rgn		= 1;
	itr		= -1;
}


RGN::RGN( const char *c, const DIR *dir )
{
	char	*s = strrchr( c, ':' );

	s[-1]	= 0;
	name	= c;
	z		= dir->ZFromName( c );
	id		= gArgs.TileIDFromName( c );
	rgn		= atoi( s + 1 );
	itr		= -1;
}


RGN::RGN( const char *path, const char *key )
{
	name	= path;
	itr		= -1;

	sscanf( key, "%d.%d:%d", &z, &id, &rgn );
}

/* --------------------------------------------------------------- */
/* CanAlign ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return RMS error from assuming that a rigid transformation
// (rot + trans) takes points {p1} into corr. points {p2}:
//
//		a x  -  b y  +  c  =  x'
//		b x  +  a y  +  d  =  y'
//
// Here there are 4 free params: a, b, c, d, and each pair of
// points (p1, p2) gives us two equations:
//
//      | p1x -p1y  1  0 |   | a |   | p2x |
//      | p1y  p1x  0  1 | x | b | = | p2y |
//                           | c |
//                           | d |
//
double RGD::CanAlign(
	const vector<Point>	&p1,
	const vector<Point>	&p2,
	bool				print )
{
// Create system of normal equations

	vector<double>	X( 4 );
	vector<double>	RHS( 4, 0.0 );
	vector<LHSCol>	LHS( 4 );
	int				np    = p1.size(),
					i1[3] = { 0, 1, 2 },
					i2[3] = { 0, 1, 3 };

	for( int i = 0; i < np; ++i ) {

		double	v1[3] = { p1[i].x, -p1[i].y, 1.0 };
		double	v2[3] = { p1[i].y,  p1[i].x, 1.0 };

		AddConstraint( LHS, RHS, 3, i1, v1, p2[i].x );
		AddConstraint( LHS, RHS, 3, i2, v2, p2[i].y );
	}

// Solve

	WriteSolveRead( X, LHS, RHS, false );

	if( print ) {

		double	mag = sqrt( X[0]*X[0] + X[1]*X[1] );

		printf( "  a=%f b=%f x0=%f y0=%f mag=%f.\n",
		X[0], X[1], X[2], X[3], mag );
	}

// Calc error

	double	sum = 0.0;

	for( int i = 0; i < np; ++i ) {

		double px = p1[i].x * X[0] - p1[i].y * X[1] + X[2];
		double py = p1[i].x * X[1] + p1[i].y * X[0] + X[3];
		double dx = px - p2[i].x;
		double dy = py - p2[i].y;

		sum += dx*dx + dy*dy;

		if( print ) {
			printf(
			"  (%8.2f %8.2f) -> (%8.2f %8.2f) =?"
			" (%8.2f %8.2f) d=%8.2f.\n",
			p1[i].x, p1[i].y, px, py,
			p2[i].x, p2[i].y, sqrt(dx*dx + dy*dy) );
		}
	}

	double	rms = sqrt( sum / np );

	if( print )
		printf( "  RMS error %f pixels.\n", rms );

	return rms;
}

/* --------------------------------------------------------------- */
/* AddPOINTPair -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Add new pair of correspondence points to r12Pts table.
//
void RGD::AddPOINTPair(
	int			r1,
	const Point	&p1,
	int			r2,
	const Point	&p2 )
{
	CRPair						r12( r1, r2 );
	map<CRPair, int>::iterator	mi = r12Idx.find( r12 );
	int							loc;

	if( mi == r12Idx.end() ) {

		// did not exist; add it
		r12Pts.push_back( CorrPts() );
		r12Idx[r12] = loc = r12Pts.size() - 1;
	}
	else
		loc = mi->second;

// Add points in correct order

	if( r1 < r2 ) {
		r12Pts[loc].pa.push_back( p1 );
		r12Pts[loc].pb.push_back( p2 );
	}
	else {
		r12Pts[loc].pa.push_back( p2 );
		r12Pts[loc].pb.push_back( p1 );
	}
}

/* --------------------------------------------------------------- */
/* TestPairAlignments -------------------------------------------- */
/* --------------------------------------------------------------- */

// For each pair of regions that we can (5 or more corr points)
// solve for a rigid transform. If the rms error exceeds tol.
// write a nasty report.
//
// There are no other material consequences to this step.
//
void RGD::TestPairAlignments()
{
	printf(
	"---- Lyr.til:rgn pairs with rigid align err > 70 pix ----\n" );

	map<CRPair, int>::iterator	pi;	// iterator over pairs

	for( pi = r12Idx.begin(); pi != r12Idx.end(); ++pi ) {

		if( r12Pts[pi->second].pa.size() >= 5 ) {

			double	RMS = CanAlign(
							r12Pts[pi->second].pa,
							r12Pts[pi->second].pb, false );

			if( RMS > 70.0 ) {	// rough limit

				const RGN	&A = vRgn[(pi->first).a];
				const RGN	&B = vRgn[(pi->first).b];

				printf( "%d.%d:%d - %d.%d:%d\n",
				A.z, A.id, A.rgn, B.z, B.id, B.rgn );

				// repeat for report
				CanAlign(
					r12Pts[pi->second].pa,
					r12Pts[pi->second].pb, true );
			}
		}
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* ListWeakConnections ------------------------------------------- */
/* --------------------------------------------------------------- */

// Add to r12Bad those region pairs having fewer than three
// corr. points between them. The real purpose here is to
// log suspicious pairs for later human debugging.
//
void CNX::ListWeakConnections( set<CRPair> &r12Bad )
{
	if( gArgs.use_all )
		return;

// Scan for weak pairs

	printf( "---- Scan for weakly connected pairs ----\n" );

	int	nct = cnxTbl.size();

	for( int i = 0; i < nct; ++i ) {

		const CnxEntry&	Ci = cnxTbl[i];
		int				np = Ci.nlinks.size();

		for( int j = 0; j < np; ++j ) {

			if( Ci.nlinks[j] < cnxMinLinks )
				r12Bad.insert( CRPair( i, Ci.linkto[j] ) );
		}
	}

// Print weak pairs

	if( r12Bad.size() ) {

		printf(
		"%d weak pair entries reported in 'FewCorr'.\n",
		r12Bad.size() );

		FILE	*f = FileOpenOrDie( "FewCorr", "w" );

		fprintf( f, "Lyr.til:rgn pairs with few ctl pts between\n" );

		set<CRPair>::iterator	pi;

		for( pi = r12Bad.begin(); pi != r12Bad.end(); ++pi ) {

			const RGN	&A = vRgn[pi->a];
			const RGN	&B = vRgn[pi->b];

			fprintf( f, "%d.%d:%d - %d.%d:%d\n",
			A.z, A.id, A.rgn, B.z, B.id, B.rgn );
		}

		fclose ( f );
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* MaxConnectedSet ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Find the longest interconnected set of regions, and
// fill list (ignore) with those that do not belong.
//
void CNX::MaxConnectedSet( set<int> &ignore )
{
	printf( "---- Find maximum set ----\n" );

	int	nct			= cnxTbl.size(),
		bestSeed	= -1,
		bestLen		= -1;

/* --------------------------------- */
/* Init all seeds to -1 'unassigned' */
/* --------------------------------- */

	for( int i = 0; i < nct; ++i )
		cnxTbl[i].seed = -1;

/* ------------------------------------------ */
/* Consider each entry as a possible new seed */
/* ------------------------------------------ */

	for( int i = 0; i < nct; ++i ) {

		/* ----------------- */
		/* Already assigned? */
		/* ----------------- */

		if( cnxTbl[i].seed >= 0 )
			continue;

		/* -------------- */
		/* Got a new seed */
		/* -------------- */

		printf( "Next seed %d...\n", i );

		// Push the new seed...
		// ...and, as we go, its connected
		// neighbors, for further processing.

		stack<int>	s;
		int			len = 0;

		s.push( i );

		/* --------------------------------- */
		/* Pop and process connected entries */
		/* --------------------------------- */

		while( !s.empty() ) {

			/* ------- */
			/* Get one */
			/* ------- */

			int			j   = s.top();
			CnxEntry	&Cj = cnxTbl[j];

			s.pop();

			/* ------------------ */
			/* Already processed? */
			/* ------------------ */

			if( Cj.seed >= 0 )
				continue;

			/* ----------------- */
			/* Belongs to seed i */
			/* ----------------- */

			Cj.seed = i;
			++len;

			/* --------------------------------- */
			/* Push all well-connected neighbors */
			/* --------------------------------- */

			int	nneib = Cj.linkto.size();

			for( int k = 0; k < nneib; ++k ) {

				// require cnxMinLinks corr pts
				if( Cj.nlinks[k] >= cnxMinLinks || gArgs.use_all )
					s.push( Cj.linkto[k] );
			}
		}

		/* ------------------- */
		/* This seed completed */
		/* ------------------- */

		printf( "...Got %d connections.\n", len );

		if( len > bestLen ) {
			bestSeed	= i;
			bestLen		= len;
		}

		/* --------------------------------- */
		/* If > nct/2 can't get a longer one */
		/* --------------------------------- */

		if( len * 2 >= nct )
			break;
	}

	printf( "Maximum length %d of %d candidates.\n", bestLen, nct );

/* -------------------------------- */
/* Create and print the ignore list */
/* -------------------------------- */

	if( bestLen < nct ) {

		printf( "\nIgnoring:\n" );

		for( int i = 0; i < nct; ++i ) {

			if( cnxTbl[i].seed != bestSeed ) {

				const RGN	&A = vRgn[i];

				ignore.insert( i );
				printf( "\t%d.%d:%d\n", A.z, A.id, A.rgn );
			}
		}
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* Set_itr_set_used ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Assign each regions's global-space transform index (itr).
// The itrs will be: {-1=invalid, 0+=valid} and a region's
// six transform values start at location X[itr * 6].
//
// Set the Constraint.used flags selecting those that will be fit.
//
// Also set global valid transform count gNTr.
//
void CNX::Set_itr_set_used( set<CRPair> &r12Bad, set<int> &ignore )
{
	printf( "---- Count valid regions ----\n" );

	int	nGoodC = 0, nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		int	r1 = vAllC[i].r1,
			r2 = vAllC[i].r2;

		vAllC[i].used = false;

		if( ignore.find( r1 ) != ignore.end() ||
			ignore.find( r2 ) != ignore.end() ||
			r12Bad.find( CRPair( r1, r2 ) ) != r12Bad.end() ) {

			continue;
		}

		if( vRgn[r1].itr == -1 )
			vRgn[r1].itr = gNTr++;

		if( vRgn[r2].itr == -1 )
			vRgn[r2].itr = gNTr++;

		vAllC[i].used = true;
		++nGoodC;
	}

	printf(
	"%d transforms to be determined, %d point correspondences.\n",
	gNTr, nGoodC );

	if( gNTr == 0 || nGoodC < 4 ) {
		printf( "Too few transforms or constraints.\n" );
		exit( 42 );
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* AddCorrespondence --------------------------------------------- */
/* --------------------------------------------------------------- */

// Update global cnxTbl whose entries track which regions connect
// to this one and how many points compose that connection.
//
// Note cnxTbl.size() kept same as vRgn.size().
//
void CNX::AddCorrespondence( int r1, int r2 )
{
// Make room

	if( r1 >= cnxTbl.size() || r2 >= cnxTbl.size() )
		cnxTbl.resize( max( r1, r2 ) + 1 );

	CnxEntry	&C1 = cnxTbl[r1],
				&C2 = cnxTbl[r2];
	int			j, n;

// Counts for r1

	n = C1.linkto.size();

	for( j = 0; j < n; ++j ) {

		if( C1.linkto[j] == r2 ) {
			++C1.nlinks[j];
			goto do_n2;
		}
	}

	C1.linkto.push_back( r2 );
	C1.nlinks.push_back( 1 );

// Counts for r2

do_n2:
	n = C2.linkto.size();

	for( j = 0; j < n; ++j ) {

		if( C2.linkto[j] == r1 ) {
			++C2.nlinks[j];
			return;
		}
	}

	C2.linkto.push_back( r1 );
	C2.nlinks.push_back( 1 );
}

/* --------------------------------------------------------------- */
/* SelectIncludedImages ------------------------------------------ */
/* --------------------------------------------------------------- */

// Decide which regions to include based on connection strength.
//
// The set<> container holds sorted, unique key-values.
//
void CNX::SelectIncludedImages()
{
	set<CRPair>	r12Bad;		// suspicious region pairs
	set<int>	ignore;		// unconnected regions

	ListWeakConnections( r12Bad );
	MaxConnectedSet( ignore );
	Set_itr_set_used( r12Bad, ignore );
}

/* --------------------------------------------------------------- */
/* FindOrAdd ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return zero-based index for given RGN.
//
// If already stored, that index is returned. Else, new
// entries are created with index (nr); nr is incremented.
//
static int FindOrAdd( map<ZIDR, int> &m, int &nr, const RGN &R )
{
// Already mapped?

	ZIDR						key( R );
	map<ZIDR, int>::iterator	it = m.find( key );

	if( it != m.end() )
		return it->second;

// No - add to map

	m[key] = nr;

// Add to vRgn

	vRgn.push_back( R );

	return nr++;
}

/* --------------------------------------------------------------- */
/* ReadPtsFile --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadPtsFile( CNX *cnx, RGD *rgd, const DIR *dir )
{
	printf( "---- Read pts ----\n" );

	FILE	*f = FileOpenOrDie( gArgs.pts_file, "r" );

	map<ZIDR, int>	getRGN;
	int				nr = 0, nlines = 0;

	for(;;) {

		char	line[4096], name1[2048], name2[2048];

		if( !fgets( line, sizeof(line), f ) )
			break;

		++nlines;

		if( !strncmp( line, "POINT", 5 ) ) {

			Point	p1, p2;

			if( 6 != sscanf( line + 6,
						"%s %lf %lf %s %lf %lf",
						name1, &p1.x, &p1.y,
						name2, &p2.x, &p2.y ) ) {

				printf(
				"WARNING: 'POINT' format error; line %d.\n",
				nlines );

				continue;
			}

			RGN	R1( name1, dir );
			RGN	R2( name2, dir );
			int r1 = FindOrAdd( getRGN, nr, R1 );
			int r2 = FindOrAdd( getRGN, nr, R2 );

			cnx->AddCorrespondence( r1, r2 );
			rgd->AddPOINTPair( r1, p1, r2, p2 );

			vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
		}
		else if( !strncmp( line, "CPOINT", 6 ) ) {

			char	key1[32], key2[32];
			Point	p1, p2;

			if( 8 != sscanf( line + 7,
						"'%[^']' %s %lf %lf '%[^']' %s %lf %lf",
						name1, key1, &p1.x, &p1.y,
						name2, key2, &p2.x, &p2.y ) ) {

				printf(
				"WARNING: 'CPOINT' format error; line %d.\n",
				nlines );

				continue;
			}

			RGN	R1( name1, key1 );
			RGN	R2( name2, key2 );
			int r1 = FindOrAdd( getRGN, nr, R1 );
			int r2 = FindOrAdd( getRGN, nr, R2 );

			cnx->AddCorrespondence( r1, r2 );
			rgd->AddPOINTPair( r1, p1, r2, p2 );

			vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
		}
		else if( !strncmp( line, "IMAGESIZE", 9 ) ) {

			if( 2 != sscanf( line + 10, "%d %d", &gW, &gH ) ) {
				printf( "Bad IMAGESIZE line '%s'.\n", line );
				exit( 42 );
			}

			fprintf( FOUT, line );
			printf( line );
		}
		else if( !strncmp( line, "FOLDMAP", 7 ) ) {

			int	z, id, nrgn = -1;

			sscanf( line + 8,
			"'%*[^']' ../%d/%d%*s %d", &z, &id, &nrgn );

			nConRgn[ZID( z, id )] = nrgn;

			fprintf( FOUT, line );
		}
		else {

			printf(
			"WARNING: Unknown entry type; line %d.\n",
			nlines );
		}
	}

	fclose( f );

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* SetPointPairs ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SetPointPairs(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
{
	int	nc	= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		double	fz =
		(vRgn[C.r1].z == vRgn[C.r2].z ? gArgs.same_strength : 1);

		double	x1 = C.p1.x * fz / sc,
				y1 = C.p1.y * fz / sc,
				x2 = C.p2.x * fz / sc,
				y2 = C.p2.y * fz / sc;
		int		j  = vRgn[C.r1].itr * 6,
				k  = vRgn[C.r2].itr * 6;

		// T1(p1) - T2(p2) = 0

		double	v[6]  = {x1,  y1,   fz, -x2, -y2,  -fz};
		int		i1[6] = {j,   j+1, j+2, k,   k+1,  k+2};
		int		i2[6] = {j+3, j+4, j+5, k+3, k+4,  k+5};

		AddConstraint( LHS, RHS, 6, i1, v, 0.0 );
		AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
	}
}

/* --------------------------------------------------------------- */
/* SetIdentityTForm ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Explicitly set some TForm to Identity.
// @@@ Does it matter which one we use?
//
static void SetIdentityTForm(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
	double	stiff	= 1.0;

	double	one	= stiff;
	int		j	= (gNTr / 2) * 6;

	AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;

// Report which tile we set

#if 1
	int	nr = vRgn.size();

	for( int k = 0; k < nr; ++k ) {

		if( vRgn[k].itr == gNTr / 2 ) {

			printf( "Ref region z=%d, id=%d.\n",
			vRgn[k].z, vRgn[k].id );
			break;
		}
	}
#endif
}

/* --------------------------------------------------------------- */
/* SetUniteLayer ------------------------------------------------- */
/* --------------------------------------------------------------- */

// ----------------------------------------

// Map {z,id,rgn} <-> TForm

class MZIDR {

public:
	int	z, id, rgn;

public:
	MZIDR() {};

	MZIDR( const RGN &R )
		{z = R.z; id = R.id; rgn = R.rgn;};

	bool operator < (const MZIDR &rhs) const
		{
			if( z < rhs.z )
				return true;
			if( z > rhs.z )
				return false;
			if( id < rhs.id )
				return true;
			if( id > rhs.id )
				return false;

			return rgn < rhs.rgn;
		};

	bool operator == (const MZIDR &rhs) const
		{return z == rhs.z && id == rhs.id && rgn == rhs.rgn;};
};

// ----------------------------------------


// Set one layer-full of TForms to those from a previous
// solution output file gArgs.tfm_file.
//
static void SetUniteLayer(
	vector<LHSCol>	&LHS,
	vector<double>	&RHS,
	double			sc )
{
/* ------------------------------- */
/* Load TForms for requested layer */
/* ------------------------------- */

	map<MZIDR, TForm>	M;

	FILE	*f = FileOpenOrDie( gArgs.tfm_file, "r" );

	for(;;) {

		char	line[256];

		if( !fgets( line, sizeof(line), f ) )
			break;

		MZIDR	R;
		TForm	T;

		R.z = atoi( line );

		if( R.z < gArgs.unite_layer )
			continue;

		if( R.z > gArgs.unite_layer )
			break;

		sscanf( line, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&R.z, &R.id, &R.rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5] );

		M[R] = T;
	}

	fclose( f );

/* ----------------------------- */
/* Set each TForm in given layer */
/* ----------------------------- */

	double	stiff	= 10.0;

	int	nr = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		if( vRgn[i].z != gArgs.unite_layer || vRgn[i].itr < 0 )
			continue;

		map<MZIDR, TForm>::iterator		it;

		it = M.find( MZIDR( vRgn[i] ) );

		if( it == M.end() )
			continue;

		double	one	= stiff,
				*t	= it->second.t;
		int		j	= vRgn[i].itr * 6;

		AddConstraint( LHS, RHS, 1, &j, &one, one*t[0] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[1] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[2] / sc );	j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[3] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[4] );		j++;
		AddConstraint( LHS, RHS, 1, &j, &one, one*t[5] / sc );	j++;
	}
}

/* --------------------------------------------------------------- */
/* PrintMagnitude ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void PrintMagnitude( const vector<double> &X, int nvars )
{
	int	k = X.size() - nvars;

	if( k >= 0 ) {

		double	mag	= sqrt( X[k]*X[k] + X[k+1]*X[k+1] );

		printf( "Final magnitude is %f = %.6e.\n", mag, mag );
	}
}

/* --------------------------------------------------------------- */
/* SolveWithSquareness ------------------------------------------- */
/* --------------------------------------------------------------- */

static void SolveWithSquareness(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
/* -------------------------- */
/* Add squareness constraints */
/* -------------------------- */

	double	stiff = gArgs.square_strength;

	for( int i = 0; i < gNTr; ++i ) {

		int	j = i * 6;

		// equal cosines
		{
			double	V[2] = {stiff, -stiff};
			int		I[2] = {j, j+4};

			AddConstraint( LHS, RHS, 2, I, V, 0.0 );
		}

		// opposite sines
		{
			double	V[2] = {stiff, stiff};
			int		I[2] = {j+1, j+3};

			AddConstraint( LHS, RHS, 2, I, V, 0.0 );
		}
	}

/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

	printf( "Solve with [transform squareness].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X, 6 );
}

/* --------------------------------------------------------------- */
/* SolveWithMontageSqr ------------------------------------------- */
/* --------------------------------------------------------------- */

// Add some constraints so the left edge of the array
// needs to be the same length as the right edge, and
// repeat for top and bottom. Of course, this assumes
// a symmetric montage.
//
static void SolveWithMontageSqr(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
	const int	MAXINT = 0x7FFFFFFF;

	int	nr = vRgn.size();

/* ------------ */
/* Which layer? */
/* ------------ */

	if( gArgs.ref_layer < 0 ) {

		gArgs.ref_layer = MAXINT;

		for( int i = 0; i < nr; ++i ) {

			if( vRgn[i].itr >= 0 ) {

				gArgs.ref_layer =
					min( gArgs.ref_layer, vRgn[i].z );
			}
		}

		printf(
		"\nNo reference layer specified,"
		" using lowest layer %d.\n", gArgs.ref_layer );
	}

/* ---------------------------------------- */
/* Search for greatest outward translations */
/* ---------------------------------------- */

	double	vNW = -MAXINT,
			vNE = -MAXINT,
			vSE = -MAXINT,
			vSW = -MAXINT;
	int		iNW, iNE, iSE, iSW,
			jNW, jNE, jSE, jSW;

	for( int i = 0; i < nr; ++i ) {

		if( vRgn[i].z != gArgs.ref_layer )
			continue;

		if( vRgn[i].itr < 0 )
			continue;

		double	v;
		int		j = vRgn[i].itr * 6;

		if( (v =  X[j+2] + X[j+5]) > vNE )
			{iNE = i; jNE = j; vNE = v;}

		if( (v =  X[j+2] - X[j+5]) > vSE )
			{iSE = i; jSE = j; vSE = v;}

		if( (v = -X[j+2] + X[j+5]) > vNW )
			{iNW = i; jNW = j; vNW = v;}

		if( (v = -X[j+2] - X[j+5]) > vSW )
			{iSW = i; jSW = j; vSW = v;}
	}

	printf(
	"Corner tiles are:"
	" se %d (%f %f),"
	" ne %d (%f %f),"
	" nw %d (%f %f),"
	" sw %d (%f %f).\n",
	vRgn[iSE].id, X[jSE+2], X[jSE+5],
	vRgn[iNE].id, X[jNE+2], X[jNE+5],
	vRgn[iNW].id, X[jNW+2], X[jNW+5],
	vRgn[iSW].id, X[jSW+2], X[jSW+5] );

/* ------------------------------------------- */
/* Use these corner tiles to impose squareness */
/* ------------------------------------------- */

	double	stiff = gArgs.square_strength;

// Top = bottom (DX = DX)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+2,  jSW+2,  jNE+2, jNW+2};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

// Left = right (DY = DY)

	{
		double	V[4] = {stiff, -stiff, -stiff, stiff};
		int		I[4] = {jSE+5,  jSW+5,  jNE+5, jNW+5};

		AddConstraint( LHS, RHS, 4, I, V, 0.0 );
	}

/* --------------- */
/* Update solution */
/* --------------- */

	printf( "Solve with [montage squareness].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	PrintMagnitude( X, 6 );
}

/* --------------------------------------------------------------- */
/* SolveWithUnitMag ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Effectively, we want to constrain the cosines and sines
// so that c^2 + s^2 = 1. We can't make constraints that are
// non-linear in the variables X[], but we can construct an
// approximation using the {c,s = X[]} of the previous fit:
// c*x + s*y = 1. To reduce sensitivity to the sizes of the
// previous fit c,s, we normalize them by m = sqrt(c^2 + s^2).
//
static void SolveWithUnitMag(
	vector<double>	&X,
	vector<LHSCol>	&LHS,
	vector<double>	&RHS )
{
	double	stiff = gArgs.scale_strength;

	for( int i = 0; i < gNTr; ++i ) {

		int		j = i * 6;
		double	c = X[j];
		double	s = X[j+3];
		double	m = sqrt( c*c + s*s );

		// c/m*x + s/m*y = 1

		double	V[2] = {c * stiff, s * stiff};
		int		I[2] = {j, j+3};

		AddConstraint( LHS, RHS, 2, I, V, m * stiff );
	}

	printf( "Solve with [unit magnitude].\n" );
	WriteSolveRead( X, LHS, RHS, false );
	printf( "\t\t\t\t" );
	PrintMagnitude( X, 6 );
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Build and solve system of linear equations.
//
// Note:
// All translational variables are scaled down by 'scale' so they
// are sized similarly to the sine/cosine variables. This is only
// to stabilize solver algorithm. We undo the scaling on exit.
//
static void SolveSystem( vector<double> &X )
{
	double	scale	= 2 * max( gW, gH );
	int		nvars	= gNTr * 6;

	printf( "%d unknowns; %d constraints.\n", nvars, vAllC.size() );

	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

	X.resize( nvars );

/* ------------------ */
/* Get rough solution */
/* ------------------ */

	SetPointPairs( LHS, RHS, scale );

	if( gArgs.unite_layer < 0 )
		SetIdentityTForm( LHS, RHS );
	else
		SetUniteLayer( LHS, RHS, scale );

	SolveWithSquareness( X, LHS, RHS );

/* ----------------------------------------- */
/* Use solution to add torsional constraints */
/* ----------------------------------------- */

	if( gArgs.make_layer_square )
		SolveWithMontageSqr( X, LHS, RHS );

	SolveWithUnitMag( X, LHS, RHS );

/* --------------------------- */
/* Rescale translational terms */
/* --------------------------- */

	int	nr	= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		if( vRgn[i].itr >= 0 ) {

			int j = vRgn[i].itr * 6;

			X[j+2] *= scale;
			X[j+5] *= scale;
		}
	}
}

/* --------------------------------------------------------------- */
/* SolveSystemRigid ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Like CanAlign, solve for 4-var transforms Y, but return 6-var
// transforms X to caller.
//
// Although this version functions correctly, cross-layer matching
// is not as good as the full affine (as expected). This version is
// retained only as an example.
//
static void SolveSystemRigid( vector<double> &X )
{
	int	nvars	= gNTr * 4,
		nc		= vAllC.size(),
		nr		= vRgn.size();

	printf( "%d unknowns; %d constraints.\n", nvars, nc );

	vector<double> Y( nvars );
	vector<double> RHS( nvars, 0.0 );
	vector<LHSCol> LHS( nvars );

/* ------------------------- */
/* Add point-pair transforms */
/* ------------------------- */

// All translational variables are scaled down by 'scale' so they
// are sized similarly to the sine/cosine variables. This is only
// to stabilize solver algorithm. We undo the scaling on exit.

	double	scale = 2 * max( gW, gH );

	for( int i = 0; i < nc; ++i ) {

		const Constraint &C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		double	x1 = C.p1.x / scale,
				y1 = C.p1.y / scale,
				x2 = C.p2.x / scale,
				y2 = C.p2.y / scale;
		int		j  = vRgn[C.r1].itr * 4,
				k  = vRgn[C.r2].itr * 4;

		// T1(p1) - T2(p2) = 0

		double	v[6]  = {x1,  -y1, 1.0, -x2,  y2, -1.0};
		int		i1[6] = {j,   j+1, j+2, k,   k+1,  k+2};
		int		i2[6] = {j+1,   j, j+3, k+1,   k,  k+3};

		AddConstraint( LHS, RHS, 6, i1, v, 0.0 );

		v[1] = y1;
		v[4] = -y2;

		AddConstraint( LHS, RHS, 6, i2, v, 0.0 );
	}

/* ------------------------------- */
/* Make one the identity transform */
/* ------------------------------- */

// Explicitly set each element of some transform.
// @@@ Does it matter which one we use?

	double	stiff	= 1.0;

	double	one	= stiff;
	int		j	= (gNTr / 2) * 6;

	AddConstraint( LHS, RHS, 1, &j, &one, one );	j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;
	AddConstraint( LHS, RHS, 1, &j, &one, 0 );		j++;

/* ----------------- */
/* 1st pass solution */
/* ----------------- */

// We have enough info for first estimate of the global
// transforms. We will need these to formulate further
// constraints on the global shape and scale.

	printf( "Solve with [rigid only].\n" );
	WriteSolveRead( Y, LHS, RHS, false );
	PrintMagnitude( Y, 4 );

/* ------------------- */
/* Make montage square */
/* ------------------- */

// Add some constraints so the left edge of the array
// needs to be the same length as the right edge, and
// repeat for top and bottom.

	const int	MAXINT = 0x7FFFFFFF;

	if( gArgs.make_layer_square ) {

		/* ------------ */
		/* Which layer? */
		/* ------------ */

		if( gArgs.ref_layer < 0 ) {

			gArgs.ref_layer = MAXINT;

			for( int i = 0; i < nr; ++i ) {

				if( vRgn[i].itr >= 0 ) {

					gArgs.ref_layer =
						min( gArgs.ref_layer, vRgn[i].z );
				}
			}

			printf(
			"\nNo reference layer specified,"
			" using lowest layer %d.\n", gArgs.ref_layer );
		}

		/* ---------------------------------------- */
		/* Search for greatest outward translations */
		/* ---------------------------------------- */

		double	vNW = -MAXINT,
				vNE = -MAXINT,
				vSE = -MAXINT,
				vSW = -MAXINT;
		int		iNW, iNE, iSE, iSW,
				jNW, jNE, jSE, jSW;

		for( int i = 0; i < nr; ++i ) {

			if( vRgn[i].z != gArgs.ref_layer )
				continue;

			if( vRgn[i].itr < 0 )
				continue;

			double	v;
			int		j = vRgn[i].itr * 4;

			if( (v =  Y[j+2] + Y[j+3]) > vNE )
				{iNE = i; jNE = j; vNE = v;}

			if( (v =  Y[j+2] - Y[j+3]) > vSE )
				{iSE = i; jSE = j; vSE = v;}

			if( (v = -Y[j+2] + Y[j+3]) > vNW )
				{iNW = i; jNW = j; vNW = v;}

			if( (v = -Y[j+2] - Y[j+3]) > vSW )
				{iSW = i; jSW = j; vSW = v;}
		}

		printf(
		"Corner tiles are:"
		" se %d (%f %f),"
		" ne %d (%f %f),"
		" nw %d (%f %f),"
		" sw %d (%f %f).\n",
		vRgn[iSE].id, Y[jSE+2], Y[jSE+3],
		vRgn[iNE].id, Y[jNE+2], Y[jNE+3],
		vRgn[iNW].id, Y[jNW+2], Y[jNW+3],
		vRgn[iSW].id, Y[jSW+2], Y[jSW+3] );

		/* ------------------------------------------- */
		/* Use these corner tiles to impose squareness */
		/* ------------------------------------------- */

		stiff = gArgs.square_strength;

		// Top = bottom (DX = DX)
		{
			double	V[4] = {stiff, -stiff, -stiff, stiff};
			int		I[4] = {jSE+2,  jSW+2,  jNE+2, jNW+2};

			AddConstraint( LHS, RHS, 4, I, V, 0.0 );
		}

		// Left = right (DY = DY)
		{
			double	V[4] = {stiff, -stiff, -stiff, stiff};
			int		I[4] = {jSE+3,  jSW+3,  jNE+3, jNW+3};

			AddConstraint( LHS, RHS, 4, I, V, 0.0 );
		}

		// Update solution

		printf( "Solve with [montage squareness].\n" );
		WriteSolveRead( Y, LHS, RHS, false );
		PrintMagnitude( Y, 4 );
	}

/* ------------------- */
/* Make unit magnitude */
/* ------------------- */

// Effectively, we want to constrain the cosines and sines
// so that c^2 + s^2 = 1. We can't make constraints that are
// non-linear in the variables Y[], but we can construct an
// approximation using the {c,s = Y[]} of the previous fit:
// c*x + s*y = 1. To reduce sensitivity to the sizes of the
// previous fit c,s, we normalize them by m = sqrt(c^2 + s^2).

	stiff = gArgs.scale_strength;

	for( int i = 0; i < gNTr; ++i ) {

		int		j = i * 4;
		double	c = Y[j];
		double	s = Y[j+1];
		double	m = sqrt( c*c + s*s );

		// c/m*x + s/m*y = 1

		double	V[2] = {c * stiff, s * stiff};
		int		I[2] = {j, j+1};

		AddConstraint( LHS, RHS, 2, I, V, m * stiff );
	}

	printf( "Solve with [unit magnitude].\n" );
	WriteSolveRead( Y, LHS, RHS, false );
	printf( "\t\t\t\t" );
	PrintMagnitude( Y, 4 );

/* --------------------------- */
/* Rescale translational terms */
/* --------------------------- */

	for( int i = 0; i < nr; ++i ) {

		if( vRgn[i].itr >= 0 ) {

			int j = vRgn[i].itr * 4;

			Y[j+2] *= scale;
			Y[j+3] *= scale;
		}
	}

/* ---------------------------- */
/* Return expanded 6-var copies */
/* ---------------------------- */

	X.resize( gNTr * 6 );

	for( int i = 0; i < gNTr; ++i ) {

		double	*dst = &X[i * 6],
				*src = &Y[i * 4];

		dst[0] = src[0];
		dst[4] = src[0];
		dst[1] = -src[1];
		dst[3] = src[1];
		dst[2] = src[2];
		dst[5] = src[3];
	}
}

/* --------------------------------------------------------------- */
/* IterateInliers ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void IterateInliers( vector<double> &X )
{
	printf( "---- Iterative solver ----\n" );

/* -------------------------- */
/* Init and count constraints */
/* -------------------------- */

	int	nsame	= 0,
		ndiff	= 0,
		nc		= vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		Constraint	&C = vAllC[i];

		if( C.used ) {

			C.inlier = true;

			if( vRgn[C.r1].z == vRgn[C.r2].z )
				++nsame;
			else
				++ndiff;
		}
		else
			C.inlier = false;
	}

	printf(
	"Constraints: %d same-layer; %d diff-layer.\n", nsame, ndiff );

/* ---------------------------------------- */
/* Repeat while any new inliers or outliers */
/* ---------------------------------------- */

	int NewInlier = 1;
	int NewOutlier = 0;

	for( int pass = 1;
		pass <= gArgs.max_pass && (NewInlier || NewOutlier);
		++pass ) {

		printf( "\nPASS %d >>>>>>>>\n", pass );

		/* ----- */
		/* Solve */
		/* ----- */

		SolveSystem( X );
//		SolveSystemRigid( X );

		/* -------------------------- */
		/* Count inliers and outliers */
		/* -------------------------- */

		NewInlier = 0;
		NewOutlier = 0;

		double	sum		= 0.0,
				big_in	= 0.0,
				big_out	= 0.0;
		int		in		= 0,
				out		= 0;

		for( int i = 0; i < nc; ++i ) {

			Constraint	&C = vAllC[i];

			if( !C.used )
				continue;

			/* ----------------------------- */
			/* Global space points and error */
			/* ----------------------------- */

			TForm	T1( &X[vRgn[C.r1].itr * 6] ),
					T2( &X[vRgn[C.r2].itr * 6] );
			Point	g1 = C.p1,
					g2 = C.p2;

			T1.Transform( g1 );
			T2.Transform( g2 );

			double	err = g2.DistSqr( g1 );
			bool	old = C.inlier;

			if( C.inlier = (sqrt( err ) <= gArgs.thresh) ) {

				sum   += err;
				big_in = max( big_in, err );

				++in;
				NewInlier += !old;
			}
			else {
				big_out = max( big_out, err );

				++out;
				NewOutlier += old;
			}
		}

		/* ------- */
		/* Reports */
		/* ------- */

		printf( "\n\t\t\t\t"
		"%d new inliers; %d new outliers.\n",
		NewInlier, NewOutlier );

		printf( "\t\t\t\t"
		"%d active constraints;"
		" %d inliers (%.2f%%),"
		" %d outliers (%.2f%%).\n",
		in + out,
		in,  double(in )/(in+out)*100.0,
		out, double(out)/(in+out)*100.0 );

		// Overall error

		double		rms	= sqrt( sum / in ),
					big	= sqrt( big_in );
		const char	*flag;

		if( rms > 20.0 )
			flag = "<---------- rms!";
		else if( big > 75.0 )
			flag = "<---------- big!";
		else
			flag = "";

		printf( "\t\t\t\t"
		"RMS error %.2f, max error inlier %.2f, max outlier %.2f %s\n",
		rms, big, sqrt( big_out ), flag );
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* Bounds -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Bounds(
	double			&xbnd,
	double			&ybnd,
	vector<double>	&X )
{
	printf( "---- Global bounds ----\n" );

// For plotting, we'd like to know the global XY-bounds.
// So transform each included regions's rectangle and find
// bounds over whole set.

	const double	BIGD = 1.0e30;

	double	xmin = BIGD, xmax = -BIGD,
			ymin = BIGD, ymax = -BIGD;
	int		nr   = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		TForm			T( &X[itr * 6] );
		vector<Point>	cnr( 4 );

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			xmin = fmin( xmin, cnr[k].x );
			xmax = fmax( xmax, cnr[k].x );
			ymin = fmin( ymin, cnr[k].y );
			ymax = fmax( ymax, cnr[k].y );
		}
	}

// Translate all transforms to put global origin at ~(0,0).
// An integer change makes layer-to-layer alignment easier.

	int	xfl = int(floor( xmin )),
		yfl = int(floor( ymin ));

	xmin -= xfl;
	xmax -= xfl;
	ymin -= yfl;
	ymax -= yfl;

	for( int i = 0; i < nr; ++i ) {

		int	j = vRgn[i].itr;

		if( j >= 0 ) {
			j		*= 6;
			X[j+2]	-= xfl;
			X[j+5]	-= yfl;
		}
	}

// Open GNUPLOT files for debugging

	FILE	*fEven		= FileOpenOrDie( "pf.even", "w" ),
			*fOdd		= FileOpenOrDie( "pf.odd", "w" ),
			*fLabEven	= FileOpenOrDie( "pf.labels.even", "w" ),
			*fLabOdd	= FileOpenOrDie( "pf.labels.odd", "w" );

// Write rects and labels

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		TForm			T( &X[itr * 6] );
		vector<Point>	cnr( 4 );
		double			xmid = 0.0, ymid = 0.0;

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {
			xmid += cnr[k].x;
			ymid += cnr[k].y;
		}

		xmid /= 4.0;
		ymid /= 4.0;

		// select even or odd reporting

		FILE	*f, *flab;
		int		color;

		if( vRgn[i].z & 1 ) {
			f		= fOdd;
			flab	= fLabOdd;
			color	= 1;
		}
		else {
			f		= fEven;
			flab	= fLabEven;
			color	= 2;
		}

		// transformed rect

		for( int k = 0; k < 5; ++k )
			fprintf( f, "%f %f\n", cnr[k%4].x, cnr[k%4].y );

		fprintf( f, "\n" );

		// label

		fprintf( flab, "set label \"%d:%d.%d \" at %f,%f tc lt %d\n",
		vRgn[i].z, vRgn[i].id, vRgn[i].rgn, xmid, ymid, color );
	}

// Close files

	fclose( fLabOdd );
	fclose( fLabEven );
	fclose( fOdd );
	fclose( fEven );

// Report

	fprintf( FOUT, "BBOX %f %f %f %f\n", xmin, ymin, xmax, ymax );

	printf( "Bounds of global image are x=[%f %f] y=[%f %f].\n\n",
	xmin, xmax, ymin, ymax );

	xbnd = xmax;
	ybnd = ymax;
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTransforms(
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	printf( "---- Write transforms ----\n" );

	FILE	*f   = FileOpenOrDie( "TFormTable.txt", "w" );
	double	smin = 100.0,
			smax = 0.0,
			smag = 0.0;
	int		nr   = vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[z[i].i];

		if( I.itr < 0 )
			continue;

		int	j = I.itr * 6;

		fprintf( f, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
		I.z, I.id, I.rgn,
		X[j  ], X[j+1], X[j+2],
		X[j+3], X[j+4], X[j+5] );

		fprintf( FOUT, "TRANSFORM %s::%d %f %f %f %f %f %f\n",
		I.name.c_str(), I.rgn,
		X[j  ], X[j+1], X[j+2],
		X[j+3], X[j+4], X[j+5] );

		double	mag = sqrt( X[j]*X[j+4] - X[j+1]*X[j+3] );

		smag += mag;
		smin  = fmin( smin, mag );
		smax  = fmax( smax, mag );
	}

	fclose( f );

	printf(
	"Average magnitude=%f, min=%f, max=%f, max/min=%f.\n\n",
	smag/gNTr, smin, smax, smax/smin );
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTrakEM(
	double					xmax,
	double					ymax,
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	FILE	*f = FileOpenOrDie( "MultLayAff.xml", "w" );

	int	oid = 3;

	fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

	TrakEM2WriteDTD( f );

	fprintf( f, "<trakem2>\n" );

	fprintf( f,
	"\t<project\n"
	"\t\tid=\"0\"\n"
	"\t\ttitle=\"Project\"\n"
	"\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
	"\t/>\n" );

	fprintf( f,
	"\t<t2_layer_set\n"
	"\t\toid=\"%d\"\n"
	"\t\ttransform=\"matrix(1.0,0.0,0.0,1.0,0.0,0.0)\"\n"
	"\t\ttitle=\"Top level\"\n"
	"\t\tlayer_width=\"%.2f\"\n"
	"\t\tlayer_height=\"%.2f\"\n"
	"\t>\n",
	oid++, xmax, ymax );

	int	prev	= -1;	// will be previously written layer
	int	offset	= int(2 * gArgs.trim + 0.5);
	int	nr		= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[z[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		// changed layer
		if( z[i].z != prev ) {

			if( prev != -1 )
				fprintf( f, "\t\t</t2_layer>\n" );

			fprintf( f,
			"\t\t<t2_layer\n"
			"\t\t\toid=\"%d\"\n"
			"\t\t\tz=\"%d\"\n"
			"\t\t>\n",
			oid++, z[i].z );

			prev = z[i].z;
		}

		// trim trailing quotes and '::'
		// s = filename only
		char	buf[2048];
		strcpy( buf, I.name.c_str() );
		char	*p = strtok( buf, " ':\n" );
		char	*s = strrchr( p, '/' );
		s = (s ? s+1 : p);

		// fix origin : undo trimming
		int		j = I.itr * 6;
		double	x = -gArgs.trim;
		double	x_orig = X[j  ]*x + X[j+1]*x + X[j+2];
		double	y_orig = X[j+3]*x + X[j+4]*x + X[j+5];

		fprintf( f,
		"\t\t\t<t2_patch\n"
		"\t\t\t\toid=\"%d\"\n"
		"\t\t\t\twidth=\"%d\"\n"
		"\t\t\t\theight=\"%d\"\n"
		"\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
		"\t\t\t\ttitle=\"%s\"\n"
		"\t\t\t\ttype=\"%d\"\n"
		"\t\t\t\tfile_path=\"%s\"\n"
		"\t\t\t\to_width=\"%d\"\n"
		"\t\t\t\to_height=\"%d\"\n"
		"\t\t\t/>\n",
		oid++, gW - offset, gH - offset,
		X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
		s, gArgs.xml_type, p, gW - offset, gH - offset );
	}

	if( nr > 0 )
		fprintf( f,"\t\t</t2_layer>\n");

	fprintf( f, "\t</t2_layer_set>\n");
	fprintf( f, "</trakem2>\n");
	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteJython(
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	FILE	*f = FileOpenOrDie( "JythonTransforms.txt", "w" );

	fprintf( f, "transforms = {\n" );

	int	nr = vRgn.size();

	for( int i = 0, itrf = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[z[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		++itrf;

		// trim trailing quotes and '::'
		char	buf[2048];
		strcpy( buf, I.name.c_str() );
		char	*p = strtok( buf, " ':\n" );

		// fix origin : undo trimming
		int		j = I.itr * 6;
		double	x = -gArgs.trim;
		double	x_orig = X[j  ]*x + X[j+1]*x + X[j+2];
		double	y_orig = X[j+3]*x + X[j+4]*x + X[j+5];

		fprintf( f, "\"%s\" : [%f, %f, %f, %f, %f, %f]%s\n",
			p, X[j], X[j+3], X[j+1], X[j+4], x_orig, y_orig,
			(itrf == gNTr ? "" : ",") );
	}

	fprintf( f, "}\n" );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* AontoBOverlap ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return approximated fraction of a on b area overlap.
//
static double AontoBOverlap( TForm &a, TForm &b )
{
	TForm			T;	// T = a->b
	vector<Point>	corners( 4 );

	AToBTrans( T, a, b );

	corners[0] = Point(  0.0, 0.0 );
	corners[1] = Point( gW-1, 0.0 );
	corners[2] = Point( gW-1, gH-1 );
	corners[3] = Point(  0.0, gH-1 );

// bounding box.

	double xmin = 1E9, xmax = -1E9;
	double ymin = 1E9, ymax = -1E9;

	for( int k = 0; k < 4; ++k ) {

		T.Transform( corners[k] );

		xmin = fmin( xmin, corners[k].x );
		ymin = fmin( ymin, corners[k].y );
		xmax = fmax( xmax, corners[k].x );
		ymax = fmax( ymax, corners[k].y );
	}

// any overlap possibility?

	if( xmin > gW-1 || xmax < 0 || ymin > gH-1 || ymax < 0 )
		return 0.0;

// approximate area using sampling of random b-points.

	const int count = 4000;

	double	wf	= double(gW-1) / RAND_MAX;
	double	hf	= double(gH-1) / RAND_MAX;
	int		in	= 0;

	for( int i = 0; i < count; ++i ) {

		Point p( wf*rand(), hf*rand() );

		T.Transform( p );

		if( p.x >= 0 && p.x < gW && p.y >= 0.0 && p.y < gH )
			++in;
	}

	//printf( "----AontoBOverlap fraction %f.\n", double(in)/count );

	return double(in)/count;
}

/* --------------------------------------------------------------- */
/* NoCorrs ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Examine region pairs having ZERO corr points between,
// but that might be expected to have connections.
//
// Note that these do not appear in the r12Bad or ignore lists.
// To get into either of those you had to be in the cnxTbl, and
// entries in the cnxTbl come only from 'POINT' entries. So our
// interesting cases are not among those. Moreover, we can skip
// cases having (itr < 0) because, again, those reflect r12Bad
// and ignore listings.
//
static void NoCorrs(
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	printf( "---- Check NoCorrs ----\n" );

	FILE	*fscr = FileOpenOrDie( "NoCorr", "w" ),
			*flog = FileOpenOrDie( "NoCorrLog", "w" );

/* ------------------------ */
/* Look at each region i... */
/* ------------------------ */

	int	nr = vRgn.size(), nreports = 0;

	for( int i = 0; i < nr; ++i ) {

		int i1 = z[i].i,
			z1 = z[i].z;

		fprintf( fscr, "#Start region %d, layer %d\n", i1, z1 );

		const RGN	&A = vRgn[i1];

		if( A.itr < 0 )
			continue;

		/* ---------------------------------------------- */
		/* ...Against each region j in same or next layer */
		/* ---------------------------------------------- */

		for( int j = i+1; j < nr && z[j].z <= z1+1; ++j ) {

			int i2 = z[j].i;

			const RGN	&B = vRgn[i2];

			if( B.itr < 0 )
				continue;

			// diff only by rgn?
			if( z1 == z[j].z && A.id == B.id )
				continue;

			// mapped pairs not interesting here
			if( r12Idx.find( CRPair( i1, i2 ) ) != r12Idx.end() )
				continue;

			/* ------------------------- */
			/* OK, this was never a pair */
			/* ------------------------- */

			TForm	t1( &X[A.itr * 6] ),
					t2( &X[B.itr * 6] );
			double	olap = AontoBOverlap( t1, t2 );

			if( olap <= 0.25 )
				continue;

			/* ----------------------- */
			/* But there is overlap... */
			/* ----------------------- */

			// Count conRgns for each

			int	nr1 = nConRgn[ZID( A.z, A.id )],
				nr2	= nConRgn[ZID( B.z, B.id )];

			// only consider cases without folds
			if( nr1 != 1 || nr2 != 1 )
				continue;

			/* ---------------- */
			/* Report this case */
			/* ---------------- */

			++nreports;

			fprintf( flog, "No points in common -"
			" Lyr.til:rgn %d.%d:%d - %d.%d:%d, overlap %.1f%%\n"
			" - %s\n"
			" - %s\n",
			A.z, A.id, A.rgn, B.z, B.id, B.rgn, olap*100.0,
			A.name.c_str(), B.name.c_str() );

			/* ---------------- */
			/* Report in NoCorr */
			/* ---------------- */

			// Create:
			// forward = t1 -> t2
			// reverse = t2 -> t1

			TForm	forward, reverse;

			AToBTrans( forward, t1, t2 );
			InvertTrans( reverse, forward );

			/* ---------------------------------------- */
			/* Instructions to redo A onto B (1 onto 2) */
			/* ---------------------------------------- */

			fprintf( fscr,
			"cd %d\n"
			"rm %d/%d.%d.*\n", A.z, A.id, B.z, B.id );

			fprintf( fscr, "#Transform 1->2: " );
			forward.PrintTransformAsParam( fscr, true );

			fprintf( fscr, "make -f make.up EXTRA='-F=../param.redo" );
			forward.PrintTransformAsParam( fscr );
			fprintf( fscr, "'\n" );

			fprintf( fscr, "cd ..\n" );

			/* ---------------------------------------- */
			/* Instructions to redo B onto A (2 onto 1) */
			/* ---------------------------------------- */

			fprintf( fscr,
			"cd %d\n"
			"rm %d/%d.%d.*\n", B.z, B.id, A.z, A.id );

			fprintf( fscr, "#Transform 2->1: " );
			reverse.PrintTransformAsParam( fscr, true );

			fprintf( fscr, "make -f make.down EXTRA='-F=../param.redo" );
			reverse.PrintTransformAsParam( fscr );
			fprintf( fscr, "'\n" );

			fprintf( fscr, "cd ..\n");
		}
	}

	fclose( flog );
	fclose( fscr );

	printf( "%d NoCorr cases reported.\n\n", nreports );
}

/* --------------------------------------------------------------- */
/* Tabulate ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Loop over all constraints (point-pairs) and---
//
// - convert constraint's points to global space
// - compute err = (p2-p1)^2 in global space
// - Store errs in Epnt[]
// - Sum energies = err^2 in Ergn[]
// - Record worst error data for whole layers and layer-pairs
//
// Report some summary results in log.
//
void EVL::Tabulate(
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	int				nr		= vRgn.size(),
					nc		= vAllC.size();
	vector<Error>	Ergn( nr );	// total error, by region
	double			sum		= 0.0,
					biggest	= 0.0;
	int				ne		= 0;

// Init region energies

	for( int i = 0; i < nr; ++i ) {
		Ergn[i].amt	= 0.0;
		Ergn[i].idx	= i;
	}

// Init whole layer data (size: max_layer_id + 1)

	Ein.resize( z[z.size()-1].z + 1 );
	Ebt = Ein;

// Tabulate errors per constraint and per region

	for( int i = 0; i < nc; ++i ) {

		Constraint	&C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		/* ----------------------------- */
		/* Global space points and error */
		/* ----------------------------- */

		TForm	T1( &X[vRgn[C.r1].itr * 6] ),
				T2( &X[vRgn[C.r2].itr * 6] );
		Point	&g1 = C.p1,
				&g2 = C.p2;

		T1.Transform( g1 );
		T2.Transform( g2 );

		double	err = g2.DistSqr( g1 );

		/* --------- */
		/* Reporting */
		/* --------- */

		sum += err;
		biggest = max( biggest, err );

		fprintf( FOUT, "MPOINTS %d %f %f %d %f %f\n",
		vRgn[C.r1].z, g1.x, g1.y,
		vRgn[C.r2].z, g2.x, g2.y );

		/* ------ */
		/* Epnt[] */
		/* ------ */

		Epnt.push_back( Error( err, ne++ ) );

		/* ------ */
		/* Ergn[] */
		/* ------ */

		Ergn[C.r1].amt += err;
		Ergn[C.r2].amt += err;

		/* ----------- */
		/* Whole layer */
		/* ----------- */

		SecErr	*S;
		int		z1 = vRgn[C.r1].z,
				z2 = vRgn[C.r2].z;

		if( z1 == z2 )
			S = &Ein[z1];
		else
			S = &Ebt[min( z1, z2 )];

		if( err > S->err )
			*S = SecErr( g1, g2, err, i );
	}

// Print overall error

	double		rms	= sqrt( sum / ne ),
				big	= sqrt( biggest );
	const char	*flag;

	if( rms > 20.0 )
		flag = "<---------- rms!";
	else if( big > 75.0 )
		flag = "<---------- big!";
	else
		flag = "";

	printf( "%d transforms, RMS error %.2f, max error %.2f %s\n\n",
	gNTr, rms, big, flag );

// Print 10 biggest errors

	int	istart, iend;

	printf( "Ten largest constraint errors---\n" );
	printf( "Index\tError\n" );

	sort( Epnt.begin(), Epnt.end() );

	iend	= Epnt.size();
	istart	= max( 0, iend - 10 );

	for( int i = istart; i < iend; ++i ) {

		printf( "%4d\t%f\n",
		Epnt[i].idx, sqrt( Epnt[i].amt ) );
	}

	printf( "\n" );

// Print regions with largest strain energies

	printf( "Ten largest region energies---\n" );
	printf( "      Energy\tLayer\tTile\t Rgn\tName\n" );

	sort( Ergn.begin(), Ergn.end() );

	iend	= Ergn.size();
	istart	= max( 0, iend - 10 );

	for( int i = istart; i < iend; ++i ) {

		const Error	&E = Ergn[i];
		const RGN	&I = vRgn[E.idx];

		printf( "%12.1f\t%4d\t%4d\t%4d\t%s\n",
		E.amt, I.z, I.id, I.rgn, I.name.c_str() );
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* Line ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void EVL::Line(
	FILE	*f,
	double	xfrom,
	double	yfrom,
	double	xto,
	double	yto )
{
	fprintf( f, "\n%f %f\n%f %f\n", xfrom, yfrom, xto, yto );
}

/* --------------------------------------------------------------- */
/* BoxOrCross ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Draw a box or a cross at the specified location.
//
void EVL::BoxOrCross( FILE *f, double x, double y, bool box )
{
	if( box ) {
		fprintf( f, "\n%f %f\n%f %f\n", x-20, y-20, x-20, y+20 );
		fprintf( f,   "%f %f\n%f %f\n", x-20, y+20, x+20, y+20 );
		fprintf( f,   "%f %f\n%f %f\n", x+20, y+20, x+20, y-20 );
		fprintf( f,   "%f %f\n%f %f\n", x+20, y-20, x-20, y-20 );
	}
	else {	// otherwise draw a cross
		Line( f, x-20, y-20, x+20, y+20 );
		Line( f, x-20, y+20, x+20, y-20 );
	}
}

/* --------------------------------------------------------------- */
/* Arrow --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Draw arrow head '>' pointing at g2, from g1-direction.
//
void EVL::Arrow( FILE *f, const Point &g1, const Point &g2 )
{
	double	q, x, y, dx, dy, c, s;

// 30 pixel length
// 50 degree opening angle

	dx = g1.x - g2.x;
	dy = g1.y - g2.y;
	q  = 30 / sqrt( dx*dx + dy*dy );
	s  = 25 * PI / 180;
	c  = cos( s );
	s  = sin( s );

	x  = q * (c*dx - s*dy) + g2.x;
	y  = q * (s*dx + c*dy) + g2.y;

	Line( f, x, y, g2.x, g2.y );

	s  = -s;
	x  = q * (c*dx - s*dy) + g2.x;
	y  = q * (s*dx + c*dy) + g2.y;

	Line( f, x, y, g2.x, g2.y );
}

/* --------------------------------------------------------------- */
/* Print_be_and_se_files ----------------------------------------- */
/* --------------------------------------------------------------- */

// Log the NPRNT biggest errors.
//
// Plot the NPLOT biggest errors, and all those over 75 pixels.
// pf.se contains those that were OK, for contrast.
//
// Display plot files using:
// > gnuplot
// gnuplot> plot 'pf.be' with lines
// gnuplot> exit
//
void EVL::Print_be_and_se_files( const vector<zsort> &z )
{
	const int NPRNT = 10;
	const int NPLOT = 50;

	FILE	*fbe = FileOpenOrDie( "pf.be", "w" ),
			*fse = FileOpenOrDie( "pf.se", "w" );

	int		ne		= Epnt.size(),
			nc		= vAllC.size();
	double	bigpr	= (ne > NPRNT ? Epnt[ne - NPRNT].amt : 0.0),
			bigpl	= (ne > NPLOT ? Epnt[ne - NPLOT].amt : 0.0);

	printf( "Maximum layer number is %d.\n\n", z[z.size()-1].z );

	printf( "Ten largest constraint errors---\n" );
	printf( "     Error\tLayer\tTile\t Rgn\tLayer\tTile\t Rgn\n" );

	for( int i = 0; i < nc; ++i ) {

		const Constraint	&C = vAllC[i];

		if( !C.used || !C.inlier )
			continue;

		const Point	&g1 = C.p1,
					&g2 = C.p2;

		double	err = g2.DistSqr( g1 );
		int		z1  = vRgn[C.r1].z,
				z2  = vRgn[C.r2].z;

		// print out if big enough

		if( err >= bigpr ) {

			printf( "%10.1f\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",
			sqrt( err ),
			z1, vRgn[C.r1].id, vRgn[C.r1].rgn,
			z2, vRgn[C.r2].id, vRgn[C.r2].rgn );

			printf( "%s\n%s\n",
			vRgn[C.r1].name.c_str(), vRgn[C.r2].name.c_str() );
		}

		// and plot

		if( err >= bigpl || sqrt( err ) > 75.0 ) {

			Line( fbe, g1.x, g1.y, g2.x, g2.y );

			BoxOrCross( fbe, g1.x, g1.y, !(z1 & 1) );
			BoxOrCross( fbe, g2.x, g2.y, !(z2 & 1) );

			Arrow( fbe, g1, g2 );
		}
		else
			Line( fse, g1.x, g1.y, g2.x, g2.y );

	}

	fclose( fbe );
	fclose( fse );
}

/* --------------------------------------------------------------- */
/* Print_errs_by_layer ------------------------------------------- */
/* --------------------------------------------------------------- */

void EVL::Print_errs_by_layer( const vector<zsort> &z )
{
	FILE	*f = FileOpenOrDie( "errs_by_layer.txt", "w" );

	int	zmax = z[z.size()-1].z;

	for( int i = z[0].z; i <= zmax; i++ ) {

		const SecErr	&Ei = Ein[i];
		const SecErr	&Eb = Ebt[i];

		int	it1 = 0, it2 = 0,	// in-layer tiles
			bt1 = 0, bt2 = 0,	// tween tiles
			bz1 = 0, bz2 = 0;	// tween z's

		if( Ei.idx >= 0 ) {

			const Constraint	&C = vAllC[Ei.idx];

			it1  = vRgn[C.r1].id;
			it2  = vRgn[C.r2].id;
		}

		if( Eb.idx >= 0 ) {

			const Constraint	&C = vAllC[Eb.idx];

			bz1	= vRgn[C.r1].z;
			bt1	= vRgn[C.r1].id;
			bz2	= vRgn[C.r2].z;
			bt2	= vRgn[C.r2].id;
		}

		fprintf( f,
		"Layer %4d:"
		" %3d:%3d %8.1f at (%8.1f, %8.1f),"
		" %4d:%3d <-> %4d:%3d %8.1f at (%8.1f, %8.1f)\n",
		i,
		it1, it2, sqrt( Ei.err ), Ei.loc.x, Ei.loc.y,
		bz1, bt1, bz2, bt2, sqrt( Eb.err ), Eb.loc.x, Eb.loc.y );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Evaluate ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void EVL::Evaluate(
	const vector<zsort>		&z,
	const vector<double>	&X )
{
	printf( "---- Evaluate errors ----\n" );

	Tabulate( z, X );
	Print_be_and_se_files( z );
	Print_errs_by_layer( z );

	printf( "\n" );
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

/* ------------------ */
/* Create output file */
/* ------------------ */

	FOUT = FileOpenOrDie( "simple", "w" );

/* --------------------- */
/* Read input data files */
/* --------------------- */

	CNX	*cnx = new CNX;	// collect connection data
	RGD	*rgd = new RGD;	// collect rigid alignment data
	DIR	*dir = new DIR;	// map name strings to z layers

	dir->ReadDIRFile();
	ReadPtsFile( cnx, rgd, dir );

	delete dir;

/* ------------------------- */
/* Try aligning region pairs */
/* ------------------------- */

// This logs reports about suspicious pairs
// but has no other impact on the real solver.

	rgd->TestPairAlignments();

	delete rgd;

/* ------------------------------------------- */
/* Decide which regions have valid connections */
/* ------------------------------------------- */

// Results mark the global RGN.itr fields

	cnx->SelectIncludedImages();

	delete cnx;

/* ----- */
/* Solve */
/* ----- */

// X are the global transform elements;
// six packed doubles per valid region.

	vector<double>	X;

	IterateInliers( X );

/* ------------------ */
/* Calc global bounds */
/* ------------------ */

	double	xbnd, ybnd;

	Bounds( xbnd, ybnd, X );

/* ---------------------------------------------- */
/* Sort regions by z -- for writing ordered files */
/* ---------------------------------------------- */

	int				nr = vRgn.size();
	vector<zsort>	z( nr );

	for( int i = 0; i < nr; ++i )
		z[i] = zsort( vRgn[i], i );

	sort( z.begin(), z.end() );

/* ---------------- */
/* Write transforms */
/* ---------------- */

	WriteTransforms( z, X );
	WriteTrakEM( xbnd, ybnd, z, X );
	WriteJython( z, X );

/* ---------------------------------- */
/* Report any missing correspondences */
/* ---------------------------------- */

	NoCorrs( z, X );

/* ------------------------ */
/* Assess and report errors */
/* ------------------------ */

	EVL	evl;

	evl.Evaluate( z, X );

/* ---- */
/* Done */
/* ---- */

	fclose( FOUT );

	return 0;
}


