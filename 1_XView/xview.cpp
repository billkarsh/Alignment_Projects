//
// Convert X_A_BIN or X_H_BIN tform data to viewable text or xml.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"Timer.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Rgns {
// The rgns for given layer
// indexed by 0-based 'idx0'
public:
	vector<char>	used;	// which rgns used
	map<int,int>	m;		// map id -> idx0
	int				nr,		// num rgns
					z;		// common z
public:
	int Init( int iz );
};

class XArray {
public:
	int				NE;
	vector<double>	X;
public:
	void Load();
};

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*inpath,
				*idb;
	double		degcw;
	int			zilo,
				zihi,
				type;	// {'T','X'}

public:
	CArgs()
	{
		inpath	= NULL;
		idb		= NULL;
		degcw	= 0.0;
		zilo	= 0;
		zihi	= 32768;
		type	= 'T';
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static FILE*	flog = NULL;
static Rgns		R;
static XArray	X;
static int		gW, gH;
static bool		isAff;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "xview.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: xview BINpath -idb=idbpath -z=i,j\n" );
		exit( 42 );
	}

	vector<int>	vi;

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {
			inpath = argv[i];
			isAff = (NULL != strstr( FileNamePtr( inpath ), "X_A" ));
		}
		else if( GetArgStr( idb, "-idb=", argv[i] ) )
			;
		else if( GetArg( &degcw, "-degcw=%lf", argv[i] ) )
			;
		else if( GetArg( &type, "-type=%c", argv[i] ) )
			;
		else if( GetArgList( vi, "-z=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zilo = vi[0];
				zihi = vi[1];
				fprintf( flog, "z [%d %d]\n", zilo, zihi );
			}
			else {
				fprintf( flog, "Bad format in -z [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Rgns::Init ---------------------------------------------------- */
/* --------------------------------------------------------------- */

int Rgns::Init( int iz )
{
	z  = iz;
	nr = IDBGetIDRgnMap( m, gArgs.idb, z, flog );

	if( nr )
		used.resize( nr );

	return nr;
}

/* --------------------------------------------------------------- */
/* XArray::Load -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadXBin( vector<double> &x, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/X_%c_%d.bin",
		gArgs.inpath, (X.NE == 6 ? 'A' : 'H'), z );

	f = FileOpenOrDie( buf, "rb" );
	fread( &x[0], sizeof(double), x.size(), f );
	fclose( f );
}


static void ReadUBin( vector<char> &u, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/U_%d.bin", gArgs.inpath, z );
	f = FileOpenOrDie( buf, "rb" );
	fread( &u[0], sizeof(char), u.size(), f );
	fclose( f );
}


void XArray::Load()
{
	NE = (isAff ? 6 : 8);
	X.resize( R.nr * NE );
	ReadXBin( X, R.z );
	ReadUBin( R.used, R.z );
}

/* --------------------------------------------------------------- */
/* GetXY_Aff ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetXY_Aff( DBox &B, const TAffine &Trot )
{
	B.L = BIGD;
	B.R = -BIGD;
	B.B = BIGD;
	B.T = -BIGD;

	TAffine			T;
	vector<Point>	cnr;
	Set4Corners( cnr, gW, gH );

	for( int z = gArgs.zilo; z <= gArgs.zihi; ++z ) {

		if( !R.Init( z ) )
			continue;

		X.Load();

		for( int j = 0; j < R.nr; ++j ) {

			if( !R.used[j] )
				continue;

			vector<Point>	c( 4 );
			memcpy( &c[0], &cnr[0], 4*2*sizeof(double) );
			T = Trot * X_AS_AFF( X.X, j );
			T.Transform( c );

			for( int k = 0; k < 4; ++k ) {
				B.L = fmin( B.L, c[k].x );
				B.R = fmax( B.R, c[k].x );
				B.B = fmin( B.B, c[k].y );
				B.T = fmax( B.T, c[k].y );
			}
		}
	}

	B.R = ceil( B.R - B.L + 1 );
	B.T = ceil( B.T - B.B + 1 );
}

/* --------------------------------------------------------------- */
/* GetXY_Hmy ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetXY_Hmy( DBox &B, const THmgphy &Trot )
{
	B.L = BIGD;
	B.R = -BIGD;
	B.B = BIGD;
	B.T = -BIGD;

	THmgphy			T;
	vector<Point>	cnr;
	Set4Corners( cnr, gW, gH );

	for( int z = gArgs.zilo; z <= gArgs.zihi; ++z ) {

		if( !R.Init( z ) )
			continue;

		X.Load();

		for( int j = 0; j < R.nr; ++j ) {

			if( !R.used[j] )
				continue;

			vector<Point>	c( 4 );
			memcpy( &c[0], &cnr[0], 4*2*sizeof(double) );
			T = Trot * X_AS_HMY( X.X, j );
			T.Transform( c );

			for( int k = 0; k < 4; ++k ) {
				B.L = fmin( B.L, c[k].x );
				B.R = fmax( B.R, c[k].x );
				B.B = fmin( B.B, c[k].y );
				B.T = fmax( B.T, c[k].y );
			}
		}
	}

	B.R = ceil( B.R - B.L + 1 );
	B.T = ceil( B.T - B.B + 1 );
}

/* --------------------------------------------------------------- */
/* Update_Aff ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update_Aff( const TAffine &Trot, const DBox &B )
{
	for( int j = 0; j < R.nr; ++j ) {

		if( !R.used[j] )
			continue;

		TAffine&	T = X_AS_AFF( X.X, j );

		T = Trot * T;
		T.AddXY( -B.L, -B.B );
	}
}

/* --------------------------------------------------------------- */
/* Update_Hmy ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update_Hmy( const THmgphy &Trot, const DBox &B )
{
	THmgphy	M( 1,0,-B.L, 0,1,-B.B, 0,0 );

	for( int j = 0; j < R.nr; ++j ) {

		if( !R.used[j] )
			continue;

		THmgphy&	T = X_AS_HMY( X.X, j );

		T = M * (Trot * T);
	}
}

/* --------------------------------------------------------------- */
/* WriteT_Aff ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteT_Aff()
{
	char	buf[64];
	sprintf( buf, "X_A_TXT/X_A_%d.txt", R.z );
	FILE	*f = FileOpenOrDie( buf, "w" );

	map<int,int>::iterator	mi, en = R.m.end();

	for( mi = R.m.begin(); mi != en; ) {

		int	id		= mi->first,
			j0		= mi->second,
			jlim	= (++mi == en ? R.nr : mi->second);

		for( int j = j0; j < jlim; ++j ) {

			if( !R.used[j] )
				continue;

			TAffine&	T = X_AS_AFF( X.X, j );

			fprintf( f, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
			id, j - j0 + 1,
			T.t[0], T.t[1], T.t[2],
			T.t[3], T.t[4], T.t[5] );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* WriteT_Hmy ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteT_Hmy()
{
	char	buf[64];
	sprintf( buf, "X_H_TXT/X_H_%d.txt", R.z );
	FILE	*f = FileOpenOrDie( buf, "w" );

	map<int,int>::iterator	mi, en = R.m.end();

	for( mi = R.m.begin(); mi != en; ) {

		int	id		= mi->first,
			j0		= mi->second,
			jlim	= (++mi == en ? R.nr : mi->second);

		for( int j = j0; j < jlim; ++j ) {

			if( !R.used[j] )
				continue;

			THmgphy&	T = X_AS_HMY( X.X, j );

			fprintf( f,
			"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t.12g\t%.12g\n",
			id, j - j0 + 1,
			T.t[0], T.t[1], T.t[2],
			T.t[3], T.t[4], T.t[5],
			T.t[6], T.t[7] );
		}
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ConvertA ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ConvertA()
{
	TAffine	Trot;
	DBox	B;

	if( gArgs.degcw || gArgs.type == 'X' ) {
		Trot.NUSetRot( gArgs.degcw * PI/180.0 );
		GetXY_Aff( B, Trot );
	}

	if( gArgs.type == 'T' )
		DskCreateDir( "X_A_TXT", flog );

	for( int z = gArgs.zilo; z <= gArgs.zihi; ++z ) {

		if( !R.Init( z ) )
			continue;

		X.Load();

		if( gArgs.degcw || gArgs.type == 'X' )
			Update_Aff( Trot, B );

		if( gArgs.type == 'T' )
			WriteT_Aff();
	}
}

/* --------------------------------------------------------------- */
/* ConvertH ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ConvertH()
{
	THmgphy	Trot;
	DBox	B;

	if( gArgs.degcw || gArgs.type == 'X' ) {
		Trot.NUSetRot( gArgs.degcw * PI/180.0 );
		GetXY_Hmy( B, Trot );
	}

	if( gArgs.type == 'T' )
		DskCreateDir( "X_H_TXT", flog );

	for( int z = gArgs.zilo; z <= gArgs.zihi; ++z ) {

		if( !R.Init( z ) )
			continue;

		X.Load();

		if( gArgs.degcw || gArgs.type == 'X' )
			Update_Hmy( Trot, B );

		if( gArgs.type == 'T' )
			WriteT_Hmy();
	}
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ------- */
/* Convert */
/* ------- */

	clock_t	t0 = StartTiming();

	if( !IDBGetImageDims( gW, gH, gArgs.idb, flog ) )
		exit( 42 );

	if( isAff )
		ConvertA();
	else
		ConvertH();

	StopTiming( flog, "Convert", t0 );

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


