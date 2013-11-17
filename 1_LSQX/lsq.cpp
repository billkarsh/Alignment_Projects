

// Based on Lou's 11/23/2010 copy of lsq.cpp


#include	"lsq_ReadPts.h"
#include	"lsq_MTrans.h"
#include	"lsq_MSimlr.h"
#include	"lsq_MAffine.h"
#include	"lsq_MHmgphy.h"

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"CAffineLens.h"
#include	"LinEqu.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Memory.h"
#include	"Timer.h"

#include	<math.h>

#include	<algorithm>
using namespace std;


/* --------------------------------------------------------------- */
/* CArgs_lsq ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_lsq {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID	re_id;

public:
	// xml_type values: these are ImagePlus codes:
	// AUTO			= -1
	// GRAY8		= 0
	// GRAY16		= 1
	// GRAY32		= 2
	// COLOR_256	= 3
	// COLOR_RGB	= 4
	//
	vector<MZID>	include_only;	// if given, include only these
	vector<double>	lrbt;			// if not empty, forced bbox
	double	same_strength,
			square_strength,
			scale_strength,
			scaf_strength,
			tfm_tol,			// transform uniformity tol
			thresh,				// outlier if worse than this
			trim,				// trim this off XML images
			degcw;				// rotate clockwise degrees
	char	*pts_file,
			*dir_file,
			*unt_file,
			*priorafftbl;		// start affine model from these
	int		model,				// model type {T,S,A,H}
			minMtgLinks,		// min connected neib/tile in montage
			unite_layer,
			ref_layer,
			nproc,
			max_pass,
			xml_type,
			xml_min,
			xml_max,
			viserr;				// 0, or, error scale
	bool	strings,
			lens,
			use_all,			// align even if #pts < 3/tile
			davinocorn;			// no davi bock same lyr corners

public:
	CArgs_lsq()
	{
		same_strength		= 1.0;
		square_strength		= 0.1;
		scale_strength		= 0.1;
		scaf_strength		= 0.01;
		tfm_tol				= -1.0;
		thresh				= 700.0;
		trim				= 0.0;
		degcw				= 0.0;
		pts_file			= NULL;
		dir_file			= NULL;
		unt_file			= NULL;
		priorafftbl			= NULL;
		model				= 'A';
		minMtgLinks			= 1;
		unite_layer			= -1;
		ref_layer			= -1;
		nproc				= 1;
		max_pass			= 1;
		xml_type			= 0;
		xml_min				= 0;
		xml_max				= 0;
		viserr				= 0;
		strings				= false;
		lens				= false;
		use_all				= false;
		davinocorn			= false;
	};

	void SetCmdLine( int argc, char* argv[] );

	int TileIDFromName( const char *name );
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

	typedef struct VisErr {

		double	L, R, B, T, D;
	}  VisErr;

private:
	vector<Error>	Epnt;	// each constraint
	vector<SecErr>	Ein,	// worst in- and between-layer
					Ebt;

private:
	void Tabulate(
		const vector<zsort>		&zs,
		const vector<double>	&X );

	void Line(
		FILE	*f,
		double	xfrom,
		double	yfrom,
		double	xto,
		double	yto );

	void BoxOrCross( FILE *f, double x, double y, bool box );
	void Arrow( FILE *f, const Point &g1, const Point &g2 );

	void Print_be_and_se_files( const vector<zsort> &zs );
	void Print_errs_by_layer( const vector<zsort> &zs );

	void ViseEval1(
		vector<VisErr>			&ve,
		const vector<double>	&X );

	void ViseEval2(
		vector<VisErr>			&ve,
		const vector<double>	&X );

	void BuildVise(
		double					xmax,
		double					ymax,
		const vector<zsort>		&zs,
		const vector<double>	&X );

public:
	void Evaluate(
		double					xmax,
		double					ymax,
		const vector<zsort>		&zs,
		const vector<double>	&X );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_lsq	gArgs;
static FILE			*FOUT;			// outfile: 'simple'
static MDL			*M		= NULL;
static int			gNTr	= 0;	// Set by Set_itr_set_used






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_lsq::SetCmdLine( int argc, char* argv[] )
{
// Parse command line args

	printf( "---- Read params ----\n" );

	vector<int>	vi;
	const char	*pat;

	re_id.Set( "/N" );

	if( argc < 2 ) {
		printf(
		"Usage: lsq <file of points>"
		" [<file of directory-name - layer-number correspondences>]"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		const char	*instr;

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
					MZID( vi[i], vi[i+1] ) );

				printf( "Include only %4d %4d\n", vi[i], vi[i+1] );
			}

			printf( "<-only> option not implemented; IGNORED.\n" );
		}
		else if( GetArgList( lrbt, "-lrbt=", argv[i] ) ) {

			if( 4 != lrbt.size() ) {
				printf( "Bad format in -lrbt [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else if( GetArg( &same_strength, "-same=%lf", argv[i] ) )
			printf( "Setting same-layer strength to %f\n", same_strength );
		else if( GetArg( &square_strength, "-square=%lf", argv[i] ) )
			printf( "Setting square strength to %f\n", square_strength );
		else if( GetArg( &scale_strength, "-scale=%lf", argv[i] ) )
			printf( "Setting scale strength to %f\n", scale_strength );
		else if( GetArg( &scaf_strength, "-scaf=%lf", argv[i] ) )
			printf( "Setting scaffold strength to %f\n", scaf_strength );
		else if( GetArg( &tfm_tol, "-tformtol=%lf", argv[i] ) )
			printf( "Setting tform uniformity to %f\n", tfm_tol );
		else if( GetArg( &thresh, "-threshold=%lf", argv[i] ) )
			printf( "Setting threshold to %f\n", thresh );
		else if( GetArg( &trim, "-trim=%lf", argv[i] ) )
			printf( "Setting trim amount to %f\n", trim );
		else if( GetArg( &degcw, "-degcw=%lf", argv[i] ) )
			printf( "Setting deg-cw to %f\n", degcw );
		else if( GetArgStr( instr, "-model=", argv[i] ) ) {

			model = toupper( instr[0] );

			switch( model ) {
				case 'T':
					printf( "Setting model to translation\n" );
				break;
				case 'S':
					printf( "Setting model to similarity\n" );
				break;
				case 'A':
					printf( "Setting model to affine\n" );
				break;
				case 'H':
					printf( "Setting model to homography\n" );
				break;
				default:
					printf( "Bad model [%s].\n", argv[i] );
					exit( 42 );
				break;
			}
		}
		else if( GetArgStr( instr, "-unite=", argv[i] ) ) {

			char	buf[2048];

			sscanf( instr, "%d,%s", &unite_layer, buf );
			unt_file = strdup( buf );

			printf( "Uniting: layer %d of '%s'.\n",
			unite_layer, unt_file );
		}
		else if( GetArgStr( instr, "-prior=", argv[i] ) ) {

			priorafftbl = strdup( instr );
			printf( "Prior solutions: '%s'.\n", priorafftbl );
		}
		else if( GetArg( &minMtgLinks, "-minmtglinks=%d", argv[i] ) )
			printf( "Minimum montage neib/tile %d\n", minMtgLinks );
		else if( GetArg( &ref_layer, "-refz=%d", argv[i] ) )
			printf( "Reference layer %d\n", ref_layer );
		else if( GetArg( &nproc, "-nproc=%d", argv[i] ) )
			printf( "Processors %d\n", nproc );
		else if( GetArg( &max_pass, "-pass=%d", argv[i] ) )
			printf( "Setting maximum passes to %d\n", max_pass );
		else if( GetArg( &xml_type, "-xmltype=%d", argv[i] ) )
			printf( "Setting xml image type to %d\n", xml_type );
		else if( GetArg( &xml_min, "-xmlmin=%d", argv[i] ) )
			printf( "Setting xml image min to %d\n", xml_min );
		else if( GetArg( &xml_max, "-xmlmax=%d", argv[i] ) )
			printf( "Setting xml image max to %d\n", xml_max );
		else if( GetArg( &viserr, "-viserr=%d", argv[i] ) )
			printf( "Setting visual error scale to %d\n", viserr );
		else if( IsArg( "-strings", argv[i] ) )
			strings = true;
		else if( IsArg( "-lens", argv[i] ) )
			lens = true;
		else if( IsArg( "-all", argv[i] ) ) {
			use_all = true;
			printf( "Using all correspondences.\n" );
		}
		else if( IsArg( "-davinc", argv[i] ) ) {
			davinocorn = true;
			printf( "Davi no same layer corners.\n" );
		}
		else if( GetArgStr( pat, "-p=", argv[i] ) ) {
			re_id.Set( pat );
			printf( "Setting pattern '%s'.\n", pat );
		}
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

	if( strings )
		re_id.Compile( stdout );
	else
		printf( "\n" );
}

/* --------------------------------------------------------------- */
/* TileIDFromName ------------------------------------------------ */
/* --------------------------------------------------------------- */

int CArgs_lsq::TileIDFromName( const char *name )
{
	const char	*s = FileNamePtr( name );
	int			id;

	if( !re_id.Decode( id, s ) ) {
		printf( "No tile-id found in '%s'.\n", s );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* Callback IDFromName ------------------------------------------- */
/* --------------------------------------------------------------- */

static int IDFromName( const char *name )
{
	return gArgs.TileIDFromName( name );
}

/* --------------------------------------------------------------- */
/* OutlierTform classes ------------------------------------------ */
/* --------------------------------------------------------------- */

// Some set of values that are derived from TForm

class Tfmval {

private:
	double	t0, t1, t3, t4;

public:
	Tfmval()
		{};
	Tfmval( const double *X )
		{t0=X[0];t1=X[1];t3=X[3];t4=X[4];};

	bool IsOutlier(
		const Tfmval	&ref,
		double			tol,
		int				z,
		int				id,
		FILE			*f );

	static void PrintHdr( FILE *f )
		{fprintf( f, "Lyr\tTil\tDT0\tDT1\tDT3\tDT4\n" );};
};

bool Tfmval::IsOutlier(
	const Tfmval	&ref,
	double			tol,
	int				z,
	int				id,
	FILE			*f )
{
	double dt0 = fabs( t0 - ref.t0 );
	double dt1 = fabs( t1 - ref.t1 );
	double dt3 = fabs( t3 - ref.t3 );
	double dt4 = fabs( t4 - ref.t4 );

	if( dt0 > tol || dt1 > tol || dt3 > tol || dt4 > tol ) {

		fprintf( f, "%d\t%d\t%f\t%f\t%f\t%f\n",
		z, id, dt0, dt1, dt3, dt4 );

		return true;
	}

	return false;
}


// Implement layer-wise calc for one or more values
// derived from TForm elements

class CLayerTfmvalCalc {

private:
	vector<double>	t0, t1, t3, t4;

public:
	void Add( const double *X );

	int Size()
		{return t0.size();};

	Tfmval LayerVals();

	void Reset()
		{t0.clear();t1.clear();t3.clear();t4.clear();};
};

void CLayerTfmvalCalc::Add( const double *X )
{
	t0.push_back( X[0] );
	t1.push_back( X[1] );
	t3.push_back( X[3] );
	t4.push_back( X[4] );
}

Tfmval CLayerTfmvalCalc::LayerVals()
{
	double	X[6] =
		{MedianVal( t0 ), MedianVal( t1 ), 0,
		 MedianVal( t3 ), MedianVal( t4 ), 0};

	return Tfmval( X );
}

/* --------------------------------------------------------------- */
/* MapFromZtoMedianTfmval ---------------------------------------- */
/* --------------------------------------------------------------- */

static void MapFromZtoMedianTfmval(
	map<int,Tfmval>			&mzT,
	const vector<double>	&X,
	const vector<zsort>		&zs )
{
	CLayerTfmvalCalc	LC;
	int					nr		= vRgn.size(),
						zcur	= -1;

// Loop over all RGN in z-order and compute median values

	for( int i = 0; i < nr; ++i ) {

		// If new layer then finish prev layer

		if( zs[i].z != zcur ) {

			if( LC.Size() ) {
				mzT[zcur] = LC.LayerVals();
				LC.Reset();
			}

			zcur = zs[i].z;
		}

		// get values for this tile

		const RGN&	I = vRgn[zs[i].i];

		if( I.itr >= 0 )
			LC.Add( &X[0] + I.itr * 6 );
	}

	if( LC.Size() )
		mzT[zcur] = LC.LayerVals();
}

/* --------------------------------------------------------------- */
/* MarkWildItrsInvalid ------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if any regions marked invalid here.
//
static bool MarkWildItrsInvalid(
	map<int,Tfmval>			&mzT,
	const vector<double>	&X,
	const vector<zsort>		&zs )
{
	FILE	*f	= FileOpenOrDie( "WildTFormTiles.txt", "w" );
	int		nr	= vRgn.size();
	bool	marked = false;

	Tfmval::PrintHdr( f );

	for( int i = 0; i < nr; ++i ) {

		RGN&	I = vRgn[zs[i].i];

		if( I.itr < 0 )
			continue;

		Tfmval	T( &X[0] + I.itr * 6 );

		if( T.IsOutlier( mzT.find( I.z )->second,
				gArgs.tfm_tol, I.z, I.id, f ) ) {

			I.itr	= -1;
			marked	= true;
		}
	}

	fclose( f );

// Repack transform indices (itr)

	if( marked ) {

		gNTr = 0;

		for( int i = 0; i < nr; ++i ) {

			RGN&	I = vRgn[zs[i].i];

			if( I.itr >= 0 )
				I.itr = gNTr++;
		}
	}

	return marked;
}

/* --------------------------------------------------------------- */
/* KillOulierTForms ---------------------------------------------- */
/* --------------------------------------------------------------- */

// After solving for transforms, calculate the median tform values
// for each layer, and then set (used = false) for each constraint
// referencing a tile whose tform is greater than tfm_tol from the
// median.
//
// Return true if any changes made here.
//
static bool KillOulierTForms(
	const vector<double>	&X,
	const vector<zsort>		&zs )
{
	map<int,Tfmval>	mzT;	// z-layer maps to median Tfmval

	MapFromZtoMedianTfmval( mzT, X, zs );

	bool changed = MarkWildItrsInvalid( mzT, X, zs );

// Disable referring constraints

	if( changed ) {

		int	nc = vAllC.size();

		for( int i = 0; i < nc; ++i ) {

			Constraint	&C = vAllC[i];

			if( C.used ) {

				if( vRgn[C.r1].itr < 0 || vRgn[C.r2].itr < 0 )
					C.used = false;
			}
		}
	}

	return changed;
}

/* --------------------------------------------------------------- */
/* IterateInliers ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void IterateInliers(
	vector<double>		&X,
	const vector<zsort>	&zs,
	int					nignored )
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

	double	lastrms, lastin, lastout;
	int		NewInlier  = 1;
	int		NewOutlier = 0;

	for( int pass = 1;
		pass <= gArgs.max_pass && (NewInlier || NewOutlier);
		++pass ) {

		printf( "\nPASS %d >>>>>>>>\n", pass );

		/* ----- */
		/* Solve */
		/* ----- */

		M->SolveSystem( X, gNTr );

		/* -------------------------- */
		/* Apply transform uniformity */
		/* -------------------------- */

		if( pass == 1 && gArgs.tfm_tol > 0.0 ) {

			if( KillOulierTForms( X, zs ) ) {

				printf( "\nPASS %d (Wild TFs Rmvd) >>>>>>>>\n", pass );

				M->SolveSystem( X, gNTr );
			}
		}

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

			Point	g1 = C.p1,
					g2 = C.p2;

			M->L2GPoint( g1, X, vRgn[C.r1].itr );
			M->L2GPoint( g2, X, vRgn[C.r2].itr );

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

		lastrms = sqrt( sum / in );
		lastin  = sqrt( big_in );
		lastout = sqrt( big_out );

		const char	*flag;

		if( lastrms > 20.0 )
			flag = "<---------- rms!";
		else if( lastin > 75.0 )
			flag = "<---------- big!";
		else
			flag = "";

		printf( "\t\t\t\t"
		"RMS error %.2f, max error inlier %.2f, max outlier %.2f %s\n",
		lastrms, lastin, lastout, flag );
	}

	printf( "\nFINAL RMS %.2f MAXIN %.2f MAXOUT %.2f IGNORED %d\n\n",
	lastrms, lastin, lastout, nignored );
}

/* --------------------------------------------------------------- */
/* ApplyLens ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ApplyLens( vector<double> &X, bool inv )
{
	if( !gArgs.lens )
		return;

	CAffineLens	LN;
	int			nr = vRgn.size();

	if( !LN.ReadIDB( idb ) )
		exit( 42 );

	for( int i = 0; i < nr; ++i ) {

		int	itr = vRgn[i].itr;

		if( itr < 0 )
			continue;

		const Til2Img	*m;
		RGN::GetMeta( &m, NULL, vRgn[i], vRgn[i] );

		LN.UpdateDoublesRHS( &X[itr * 6], m->cam, inv );
	}

	IDBT2ICacheClear();
}

#if 0

// @@@ Not sure this diagnostic for missing solutions is that
// informative. Also, it builds scripts for retries that are
// no longer correct.
//

/* --------------------------------------------------------------- */
/* AontoBOverlap ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return approximated fraction of a on b area overlap.
//
static double AontoBOverlap( TAffine &a, TAffine &b )
{
	TAffine			T;	// T = a->b
	vector<Point>	corners( 4 );

	T.FromAToB( a, b );

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

	//printf( "----AontoBOverlap fraction %f\n", double(in)/count );

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
	const vector<zsort>		&zs,
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

		int i1 = zs[i].i,
			z1 = zs[i].z;

		fprintf( fscr, "#Start region %d, layer %d\n", i1, z1 );

		const RGN	&A = vRgn[i1];

		if( A.itr < 0 )
			continue;

		/* ---------------------------------------------- */
		/* ...Against each region j in same or next layer */
		/* ---------------------------------------------- */

// zs[j].z <= z1+1 tests same,down,up; older practice
//		for( int j = i+1; j < nr && zs[j].z <= z1+1; ++j ) {

// zs[j].z <= z1 tests same and down, as per current practice

		for( int j = i+1; j < nr && zs[j].z <= z1; ++j ) {

			int i2 = zs[j].i;

			const RGN	&B = vRgn[i2];

			if( B.itr < 0 )
				continue;

			// diff only by rgn?
			if( z1 == zs[j].z && A.id == B.id )
				continue;

			// mapped pairs not interesting here
			if( r12Idx.find( CRPair( i1, i2 ) ) != r12Idx.end() )
				continue;

			/* ------------------------- */
			/* OK, this was never a pair */
			/* ------------------------- */

			TAffine	t1( &X[A.itr * 6] ),
					t2( &X[B.itr * 6] );
			double	olap = AontoBOverlap( t1, t2 );

			if( olap <= 0.25 )
				continue;

			/* ----------------------- */
			/* But there is overlap... */
			/* ----------------------- */

			// Count conRgns for each

			int	nr1 = nConRgn[MZID( A.z, A.id )],
				nr2	= nConRgn[MZID( B.z, B.id )];

			// only consider cases without folds
			if( nr1 != 1 || nr2 != 1 )
				continue;

			/* ---------------- */
			/* Report this case */
			/* ---------------- */

			++nreports;

			const Til2Img	*ma, *mb;
			RGN::GetMeta( &ma, &mb, A, B );

			fprintf( flog, "No points in common -"
			" Lyr.til:rgn %d.%d:%d - %d.%d:%d, overlap %.1f%%\n"
			" - %s\n"
			" - %s\n",
			A.z, A.id, A.rgn, B.z, B.id, B.rgn, olap*100.0,
			ma->path.c_str(), mb->path.c_str() );

			/* ---------------- */
			/* Report in NoCorr */
			/* ---------------- */

			// Create:
			// forward = t1 -> t2
			// reverse = t2 -> t1

			TAffine	forward, reverse;

			forward.FromAToB( t1, t2 );
			reverse.InverseOf( forward );

			/* ---------------------------------------- */
			/* Instructions to redo A onto B (1 onto 2) */
			/* ---------------------------------------- */

			fprintf( fscr,
			"cd %d\n"
			"rm %d/%d.%d.*\n", A.z, A.id, B.z, B.id );

			fprintf( fscr, "#Transform 1->2: " );
			forward.TPrintAsParam( fscr, true );

			fprintf( fscr, "make -f make.up EXTRA='-F=../param.redo" );
			forward.TPrintAsParam( fscr );
			fprintf( fscr, "'\n" );

			fprintf( fscr, "cd ..\n" );

			/* ---------------------------------------- */
			/* Instructions to redo B onto A (2 onto 1) */
			/* ---------------------------------------- */

			fprintf( fscr,
			"cd %d\n"
			"rm %d/%d.%d.*\n", B.z, B.id, A.z, A.id );

			fprintf( fscr, "#Transform 2->1: " );
			reverse.TPrintAsParam( fscr, true );

			fprintf( fscr, "make -f make.down EXTRA='-F=../param.redo" );
			reverse.TPrintAsParam( fscr );
			fprintf( fscr, "'\n" );

			fprintf( fscr, "cd ..\n");
		}
	}

	fclose( flog );
	fclose( fscr );

	printf( "%d NoCorr cases reported.\n\n", nreports );

	IDBT2ICacheClear();
}

#endif

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
// - Enter error in 'PostFitErrs.txt'
//
// Report some summary results in log.
//
void EVL::Tabulate(
	const vector<zsort>		&zs,
	const vector<double>	&X )
{
	FILE			*f		= FileOpenOrDie( "PostFitErrs.txt", "w" );
	int				nr		= vRgn.size(),
					nc		= vAllC.size();
	vector<Error>	Ergn( nr );	// total error, by region
	double			sum		= 0.0,
					biggest	= 0.0;

// 'PostFitErrs.txt' headers

	fprintf( f, "LyrA\tTilA\tRgnA\tLyrB\tTilB\tRgnB\tSqrErr\n" );

// Init region energies

	for( int i = 0; i < nr; ++i ) {
		Ergn[i].amt	= 0.0;
		Ergn[i].idx	= i;
	}

// Init whole layer data (size: max_layer_id + 1)

	Ein.resize( zs[zs.size()-1].z + 1 );
	Ebt = Ein;

// Tabulate errors per constraint and per region

	for( int i = 0; i < nc; i += 4 ) {

		double		err;
		Constraint	&C = vAllC[i];
		const RGN&	R1 = vRgn[C.r1];
		const RGN&	R2 = vRgn[C.r2];

		if( !C.used )
			err = -2;
		else if( !C.inlier )
			err = -1;
		else {

			/* ----------------------------- */
			/* Global space points and error */
			/* ----------------------------- */

			Point	&g1 = C.p1,
					&g2 = C.p2;

			M->L2GPoint( g1, X, R1.itr );
			M->L2GPoint( g2, X, R2.itr );

			err = g2.DistSqr( g1 );

			/* --------- */
			/* Reporting */
			/* --------- */

			sum += err;
			biggest = max( biggest, err );

			fprintf( FOUT, "MPOINTS %d %f %f %d %f %f\n",
			R1.z, g1.x, g1.y,
			R2.z, g2.x, g2.y );

			/* ------ */
			/* Epnt[] */
			/* ------ */

			Epnt.push_back( Error( err, i ) );

			/* ------ */
			/* Ergn[] */
			/* ------ */

			Ergn[C.r1].amt += err;
			Ergn[C.r2].amt += err;

			/* ----------- */
			/* Whole layer */
			/* ----------- */

			SecErr	*S;
			int		z1 = R1.z,
					z2 = R2.z;

			if( z1 == z2 )
				S = &Ein[z1];
			else
				S = &Ebt[min( z1, z2 )];

			if( err > S->err )
				*S = SecErr( g1, g2, err, i );
		}

		/* ---- */
		/* File */
		/* ---- */

		fprintf( f, "%d\t%d\t%d\t%d\t%d\t%d\t%f\n",
		R1.z, R1.id, R1.rgn,
		R2.z, R2.id, R2.rgn, err );
	}

// Close file

	fclose( f );

// Print overall error

	int			istart,
				iend	= Epnt.size();
	double		rms		= sqrt( sum / iend ),
				big		= sqrt( biggest );
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

	printf( "Ten largest constraint errors---\n" );
	printf( "Error\n" );

	sort( Epnt.begin(), Epnt.end() );

	istart = max( 0, iend - 10 );

	for( int i = istart; i < iend; ++i )
		printf( "%f\n", sqrt( Epnt[i].amt ) );

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
		E.amt, I.z, I.id, I.rgn, I.GetName() );
	}

	printf( "\n" );

	IDBT2ICacheClear();
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
void EVL::Print_be_and_se_files( const vector<zsort> &zs )
{
	const int NPRNT = 10;
	const int NPLOT = 50;

	FILE	*fbe = FileOpenOrDie( "pf.be", "w" ),
			*fse = FileOpenOrDie( "pf.se", "w" );

	int		ne		= Epnt.size(),
			nc		= vAllC.size();
	double	bigpr	= (ne > NPRNT ? Epnt[ne - NPRNT].amt : 0.0),
			bigpl	= (ne > NPLOT ? Epnt[ne - NPLOT].amt : 0.0);

	printf( "Maximum layer number is %d\n\n", zs[zs.size()-1].z );

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

			const Til2Img	*m1, *m2;
			RGN::GetMeta( &m1, &m2, vRgn[C.r1], vRgn[C.r2] );

			printf( "%10.1f\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n",
			sqrt( err ),
			z1, m1->id, vRgn[C.r1].rgn,
			z2, m2->id, vRgn[C.r2].rgn );

			printf( "%s\n%s\n",
			m1->path.c_str(), m2->path.c_str() );
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

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* Print_errs_by_layer ------------------------------------------- */
/* --------------------------------------------------------------- */

void EVL::Print_errs_by_layer( const vector<zsort> &zs )
{
	FILE	*f = FileOpenOrDie( "errs_by_layer.txt", "w" );

	int	zmax = zs[zs.size()-1].z;

	for( int i = zs[0].z; i <= zmax; i++ ) {

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
/* ViseWriteXML -------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	visePix	128

static void ViseWriteXML(
	double					xmax,
	double					ymax,
	const vector<zsort>		&zs,
	const vector<double>	&X )
{
	FILE	*f = FileOpenOrDie( "visexml.xml", "w" );

	double	sclx = (double)gW / visePix,
			scly = (double)gH / visePix;
	int		oid  = 3;

	fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

	TrakEM2WriteDTD( f );

	fprintf( f, "<trakem2>\n" );

	fprintf( f,
	"\t<project\n"
	"\t\tid=\"0\"\n"
	"\t\ttitle=\"Project\"\n"
	"\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
	"\t\tn_mipmap_threads=\"8\"\n"
	"\t/>\n" );

	fprintf( f,
	"\t<t2_layer_set\n"
	"\t\toid=\"%d\"\n"
	"\t\ttransform=\"matrix(1,0,0,1,0,0)\"\n"
	"\t\ttitle=\"Top level\"\n"
	"\t\tlayer_width=\"%.2f\"\n"
	"\t\tlayer_height=\"%.2f\"\n"
	"\t>\n",
	oid++, xmax, ymax );

	int	prev	= -1;	// will be previously written layer
	int	nr		= vRgn.size();

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[zs[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		// changed layer
		if( zs[i].z != prev ) {

			if( prev != -1 )
				fprintf( f, "\t\t</t2_layer>\n" );

			fprintf( f,
			"\t\t<t2_layer\n"
			"\t\t\toid=\"%d\"\n"
			"\t\t\tthickness=\"0\"\n"
			"\t\t\tz=\"%d\"\n"
			"\t\t>\n",
			oid++, zs[i].z );

			prev = zs[i].z;
		}

		char			title[256];
		const Til2Img	*m;
		RGN::GetMeta( &m, NULL, I, I );

		if( m->col != -999 ) {

			sprintf( title, "ve_%d-%d__%d-%d-%d",
				I.z, I.id, m->col, m->row, m->cam );
		}
		else
			sprintf( title, "ve_%d-%d", I.z, I.id );

		TAffine	A = M->EqvAffine( X, I.itr );

		fprintf( f,
		"\t\t\t<t2_patch\n"
		"\t\t\t\toid=\"%d\"\n"
		"\t\t\t\twidth=\"%d\"\n"
		"\t\t\t\theight=\"%d\"\n"
		"\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
		"\t\t\t\ttitle=\"%s\"\n"
		"\t\t\t\ttype=\"4\"\n"
		"\t\t\t\tfile_path=\"viseimg/%d/%s.png\"\n"
		"\t\t\t\to_width=\"%d\"\n"
		"\t\t\t\to_height=\"%d\"\n"
		"\t\t\t\tmin=\"0\"\n"
		"\t\t\t\tmax=\"255\"\n"
		"\t\t\t/>\n",
		oid++, visePix, visePix,
		sclx*A.t[0], scly*A.t[3], sclx*A.t[1], scly*A.t[4], A.t[2], A.t[5],
		title, I.z, title, visePix, visePix );
	}

	if( nr > 0 )
		fprintf( f,"\t\t</t2_layer>\n");

	fprintf( f, "\t</t2_layer_set>\n");
	fprintf( f, "</trakem2>\n");
	fclose( f );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* ViseEval1 ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Each rgn gets a VisErr record that describes errors w.r.t.
// adjacent tiles in same layer {L,R,B,T} and with down {D}.
// In this first pass, each {L,R,B,T} bin gets a negative
// count of the correspondence points that map to that side.
// The issue is to make sure that corner points will map to
// a proper side. Sides having < 3 points are actually empty.
//
void EVL::ViseEval1(
	vector<VisErr>			&ve,
	const vector<double>	&X )
{
	Point	O( gW/2, gH/2 );	// local image center
	int		ne = Epnt.size();
	double	aTL, aTR;			// top-left, right angles

	aTR = atan2( O.y, O.x ) * 180/PI;
	aTL = 180 - aTR;

// zero all {L,R,B,T,D}
	memset( &ve[0], 0, vRgn.size() * sizeof(VisErr) );

	for( int i = 0; i < ne; ++i ) {

		const Constraint	&C = vAllC[Epnt[i].idx];

		int	z1 = vRgn[C.r1].z,
			z2 = vRgn[C.r2].z;

		if( z1 == z2 ) {

			// constraints 1 & 2 done symmetrically
			// so load into arrays and loop

			const Point*	p[2] = {&C.p1, &C.p2};
			int				r[2] = {C.r1, C.r2};

			for( int j = 0; j < 2; ++j ) {

				// the constraint points are in global coords
				// so to see which side-sector a point is in
				// we will get angle of vector from local center
				// to local point

				VisErr	&V = ve[r[j]];
				Point	L = *p[j];
				double	a, *s;

				M->G2LPoint( L, X, vRgn[r[j]].itr );

				a = atan2( L.y-O.y, L.x-O.x ) * 180/PI;

				// R (-aTR,aTR]
				// T (aTR,aTL]
				// B (-aTL,-aTR]
				// L else

				if( a > -aTL ) {

					if( a > -aTR ) {

						if( a > aTR ) {

							if( a > aTL )
								s = &V.L;
							else
								s = &V.T;
						}
						else
							s = &V.R;
					}
					else
						s = &V.B;
				}
				else
					s = &V.L;

				*s -= 1;	// counts are negative
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* ViseEval2 ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Each rgn gets a VisErr record that describes errors w.r.t.
// adjacent tiles in same layer {L,R,B,T} and with down {D}.
// In this second pass, we record the maximum error for each
// bin. If a point maps to side with low occupancy from pass
// one, then we remap it to a true side.
//
void EVL::ViseEval2(
	vector<VisErr>			&ve,
	const vector<double>	&X )
{
	Point	O( gW/2, gH/2 );	// local image center
	int		ne = Epnt.size();
	double	aTL, aTR;			// top-left, right angles

	aTR = atan2( O.y, O.x ) * 180/PI;
	aTL = 180 - aTR;

	for( int i = 0; i < ne; ++i ) {

		const Constraint	&C = vAllC[Epnt[i].idx];

		int	z1 = vRgn[C.r1].z,
			z2 = vRgn[C.r2].z;

		if( z1 == z2 ) {

			// constraints 1 & 2 done symmetrically
			// so load into arrays and loop

			const Point*	p[2] = {&C.p1, &C.p2};
			int				r[2] = {C.r1, C.r2};

			for( int j = 0; j < 2; ++j ) {

				// the constraint points are in global coords
				// so to see which side-sector a point is in
				// we will get angle of vector from local center
				// to local point

				VisErr	&V = ve[r[j]];
				Point	L = *p[j];
				double	a, *s;

				M->G2LPoint( L, X, vRgn[r[j]].itr );

				a = atan2( L.y-O.y, L.x-O.x ) * 180/PI;

				// R (-aTR,aTR]
				// T (aTR,aTL]
				// B (-aTL,-aTR]
				// L else

				if( a > -aTL ) {

					if( a > -aTR ) {

						if( a > aTR ) {

							if( a > aTL )
								s = &V.L;
							else
								s = &V.T;
						}
						else
							s = &V.R;
					}
					else
						s = &V.B;
				}
				else
					s = &V.L;

				// occupancy test and remap

				if( !*s )
					continue;
				else if( *s < -2 || *s > 0 )
					;
				else {
					// remap

					if( s == &V.L ) {

						if( a >= -180 ) {
							if( V.B < -3 || V.B > 0 )
								s = &V.B;
							else
								continue;
						}
						else {
							if( V.T < -3 || V.T > 0 )
								s = &V.T;
							else
								continue;
						}
					}
					else if( s == &V.R ) {

						if( a <= 0 ) {
							if( V.B < -3 || V.B > 0 )
								s = &V.B;
							else
								continue;
						}
						else {
							if( V.T < -3 || V.T > 0 )
								s = &V.T;
							else
								continue;
						}
					}
					else if( s == &V.B ) {

						if( a <= -90 ) {
							if( V.L < -3 || V.L > 0 )
								s = &V.L;
							else
								continue;
						}
						else {
							if( V.R < -3 || V.R > 0 )
								s = &V.R;
							else
								continue;
						}
					}
					else {

						if( a >= 90 ) {
							if( V.L < -3 || V.L > 0 )
								s = &V.L;
							else
								continue;
						}
						else {
							if( V.R < -3 || V.R > 0 )
								s = &V.R;
							else
								continue;
						}
					}
				}

				if( Epnt[i].amt > *s )
					*s = Epnt[i].amt;
			}
		}
		else {

			// for down, we will just lump all together into D

			double	&D = (z1 > z2 ? ve[C.r1].D : ve[C.r2].D);

			if( Epnt[i].amt > D )
				D = Epnt[i].amt;
		}
	}
}

/* --------------------------------------------------------------- */
/* ViseColor ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static uint32 ViseColor( double err )
{
	const uint32	ezero = 0xFF0000FF,
					eover = 0xFF00FFFF;

	uint32	ecolr;

	if( err <= 0 )
		ecolr = ezero;
	else if( err >= gArgs.viserr )
		ecolr = eover;
	else {
		uint32	c = (uint32)(255 * err / gArgs.viserr);
		ecolr = 0xFF000000 + (c << 8);
	}

	return ecolr;
}

/* --------------------------------------------------------------- */
/* VisePaintRect ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void VisePaintRect(
	vector<uint32>	&RGB,
	int				x0,
	int				xlim,
	int				y0,
	int				ylim,
	double			errsqr )
{
	uint32	ecolr = ViseColor( (errsqr > 0 ? sqrt( errsqr ) : 0) );

	++xlim;
	++ylim;

	for( int y = y0; y < ylim; ++y ) {

		for( int x = x0; x < xlim; ++x )
			RGB[x + visePix * y] = ecolr;
	}
}

/* --------------------------------------------------------------- */
/* BuildVise ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Build a set of error visualization images--
//
// visexml.xml is TrackEM2 file pointing at diagnostic images
// in viseimg/. Each RGB image is visePix^2 pixels. The central
// square in each image depicts down errors {red=none, green
// brightness proportional to max error}. The sides of the open
// box in each image are similarly scaled for in-plane errors.
//
void EVL::BuildVise(
	double					xmax,
	double					ymax,
	const vector<zsort>		&zs,
	const vector<double>	&X )
{
	ViseWriteXML( xmax, ymax, zs, X );

	int				nr = vRgn.size();
	vector<VisErr>	ve( nr );

	ViseEval1( ve, X );
	ViseEval2( ve, X );

	DskCreateDir( "viseimg", stdout );

	FILE	*f = FileOpenOrDie( "viseparse.txt", "w" );
	char	buf[256];
	int		prev = -1,	// will be previously written layer
			twv  = visePix / 12;

	fprintf( f, "Z\tID\tCol\tRow\tCam\tL\tR\tB\tT\tD\n" );

	for( int i = 0; i < nr; ++i ) {

		const RGN&	I = vRgn[zs[i].i];

		// skip unused tiles
		if( I.itr < 0 )
			continue;

		// changed layer
		if( zs[i].z != prev ) {
			sprintf( buf, "viseimg/%d", zs[i].z );
			DskCreateDir( buf, stdout );
			prev = zs[i].z;
		}
		// light gray image
		vector<uint32>	RGB( visePix * visePix, 0xFFD0D0D0 );

		const VisErr	&V = ve[zs[i].i];
		const Til2Img	*m;
		int				col = 0, row = 0, cam = 0;

		RGN::GetMeta( &m, NULL, I, I );

		if( m->col != -999 )
			col = m->col, row = m->row, cam = m->cam;

		fprintf( f, "%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n",
		I.z, I.id, col, row, cam,
		(V.L > 0 ? sqrt(V.L) : 0),
		(V.R > 0 ? sqrt(V.R) : 0),
		(V.B > 0 ? sqrt(V.B) : 0),
		(V.T > 0 ? sqrt(V.T) : 0),
		(V.D > 0 ? sqrt(V.D) : 0) );

		// down
		if( zs[i].z != zs[0].z ) {

			VisePaintRect( RGB,
				5 * twv, 7 * twv,
				5 * twv, 7 * twv,
				V.D );
		}

		//left
		VisePaintRect( RGB,
			3 * twv, 4 * twv,
			4 * twv, 8 * twv,
			V.L );

		//right
		VisePaintRect( RGB,
			8 * twv, 9 * twv,
			4 * twv, 8 * twv,
			V.R );

		//bot
		VisePaintRect( RGB,
			4 * twv, 8 * twv,
			3 * twv, 4 * twv,
			V.B );

		//top
		VisePaintRect( RGB,
			4 * twv, 8 * twv,
			8 * twv, 9 * twv,
			V.T );

		// border
		int	lim = visePix - 1;
		VisePaintRect( RGB, 0, twv/2, 0, lim, .01 );
		VisePaintRect( RGB, lim - twv/2, lim, 0, lim, .01 );
		VisePaintRect( RGB, 0, lim, 0, twv/2, .01 );
		VisePaintRect( RGB, 0, lim, lim - twv/2, lim, .01 );

		// store
		if( m->col != -999 ) {
			sprintf( buf, "viseimg/%d/ve_%d-%d__%d-%d-%d.png",
			I.z, I.z, I.id, col, row, cam );
		}
		else {
			sprintf( buf, "viseimg/%d/ve_%d-%d.png",
			I.z, I.z, I.id );
		}

		Raster32ToPngRGBA( buf, &RGB[0], visePix, visePix );
	}

	fclose( f );

	IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* Evaluate ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void EVL::Evaluate(
	double					xmax,
	double					ymax,
	const vector<zsort>		&zs,
	const vector<double>	&X )
{
	printf( "---- Evaluate errors ----\n" );

	Tabulate( zs, X );
//	Print_be_and_se_files( zs );
//	Print_errs_by_layer( zs );

	if( gArgs.viserr > 0 )
		BuildVise( xmax, ymax, zs, X );

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern void Test();

int main( int argc, char **argv )
{
	//{
	//	Test();
	//	exit( 0 );
	//}

	//{
	//	MAffine	M;
	//	M.Test();
	//	exit( 0 );
	//}

	clock_t	t0 = StartTiming();

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

// CNX: collect connection data
// SML: collect data to try similarity alignments

	CNX	*cnx = new CNX;
	SML	*sml = new SML;

	if( gArgs.strings ) {

		ReadPts_StrTags( FOUT, cnx, sml,
			IDFromName, gArgs.dir_file, gArgs.pts_file );
	}
	else
		ReadPts_NumTags( FOUT, cnx, sml, gArgs.pts_file,
			gArgs.davinocorn );

	t0 = StopTiming( stdout, "ReadPts", t0 );fflush( stdout );

/* ----------------- */
/* Sort regions by z */
/* ----------------- */

	int				nr = vRgn.size();
	vector<zsort>	zs( nr );

	for( int i = 0; i < nr; ++i )
		zs[i] = zsort( vRgn[i], i );

	sort( zs.begin(), zs.end() );

	t0 = StopTiming( stdout, "Sort", t0 );fflush( stdout );

	printf( "Z range [%d %d]\n\n", zs[0].z, zs[nr-1].z );

/* ------------------------- */
/* Try aligning region pairs */
/* ------------------------- */

// This writes logs about suspicious pairs
// but has no other impact on real solver.

//	sml->TestPairAlignments();

	delete sml;

	t0 = StopTiming( stdout, "TestPairs", t0 );fflush( stdout );

/* ------------ */
/* Select model */
/* ------------ */

	switch( gArgs.model ) {
		case 'T': M = new MTrans;  break;
		case 'S': M = new MSimlr;  break;
		case 'A': M = new MAffine; break;
		case 'H': M = new MHmgphy; break;
	}

/* ------------------------------------------- */
/* Decide which regions have valid connections */
/* ------------------------------------------- */

// Results mark the global RGN.itr fields

	int	nignored;

	gNTr = cnx->SelectIncludedImages(
			nignored, gArgs.minMtgLinks,
			(gArgs.use_all ? 0 : M->MinPairs()) );

	delete cnx;

	t0 = StopTiming( stdout, "CNX", t0 );fflush( stdout );

/* ----- */
/* Solve */
/* ----- */

// X are the packed global transform elements;
// E.g. six doubles per valid region for affines.

	vector<double>	X;

	M->SetModelParams(
		gW, gH,
		gArgs.same_strength,
		gArgs.square_strength,
		gArgs.scale_strength,
		gArgs.scaf_strength,
		gArgs.nproc,
		gArgs.unite_layer,
		gArgs.unt_file,
		gArgs.priorafftbl,
		&zs );

//static_cast<MAffine*>(M)->UpdateScaffold( X, gNTr );
//goto exit;

	IterateInliers( X, zs, nignored );
	ApplyLens( X, false );

/* ------------------ */
/* Calc global bounds */
/* ------------------ */

	double	xbnd, ybnd;

	M->Bounds( xbnd, ybnd, X,
		gArgs.lrbt, gArgs.degcw, FOUT );

/* ---------------- */
/* Write transforms */
/* ---------------- */

	M->WriteTransforms( X, gArgs.strings, FOUT );

	M->WriteTrakEM( xbnd, ybnd, X, gArgs.trim,
		gArgs.xml_type, gArgs.xml_min, gArgs.xml_max );

	M->WriteJython( X, gArgs.trim, gNTr );

/* ---------------------------------- */
/* Report any missing correspondences */
/* ---------------------------------- */

	ApplyLens( X, true );
//	NoCorrs( zs, X );

/* ------------------------ */
/* Assess and report errors */
/* ------------------------ */

	t0 = StartTiming();

	{
		EVL	evl;

		evl.Evaluate( xbnd, ybnd, zs, X );
	}

	StopTiming( stdout, "EvalErr", t0 );fflush( stdout );

/* ---- */
/* Done */
/* ---- */

exit:
	fclose( FOUT );
	VMStats( stdout );

	return 0;
}


