//
// Read in one or two 'Error' folders (A,B) from lsqw output.
//
// Create output file 'eview.txt' recording histogram of
// errors for:
//
//	Err	Aall	Asam	Adwn	Ball	Bsam	Bdwn
//
// All histograms are binwidth=1/div, nbins = div*lim + 1overflow.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CHst {
public:
	long	*all, *sam, *dwn;
public:
	void Read( const char *path );
};

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {
public:
	char	*inA, *inB;
	int		zilo, zihi,
			lim,  div;
public:
	CArgs() :
	inA(NULL), inB(NULL),
	zilo(0),   zihi(32768),
	lim(500),  div(10) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static FILE*	flog	= NULL;
static CHst		A, B;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "eview.txt", "w" );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: eview fileA [fileB] -z=i,j [options].\n" );
		exit( 42 );
	}

	vector<int>	vi;

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( !inA )
				inA = argv[i];
			else
				inB = argv[i];
		}
		else if( GetArgList( vi, "-z=", argv[i] ) ) {

			if( 2 == vi.size() ) {
				zilo = vi[0];
				zihi = vi[1];
			}
			else {
				fprintf( flog, "Bad format in -z [%s].\n", argv[i] );
				exit( 42 );
			}
		}
		else if( GetArg( &lim, "-lim=%d", argv[i] ) )
			;
		else if( GetArg( &div, "-div=%d", argv[i] ) )
			;
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}

// headers

	fprintf( flog, "Err\tAall\tAsam\tAdwn" );

	if( inB )
		fprintf( flog, "\tBall\tBsam\tBdwn\n" );
	else
		fprintf( flog, "\n" );
}

/* --------------------------------------------------------------- */
/* Read ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CHst::Read( const char *path )
{
	int	emax	= gArgs.lim * gArgs.div;
	int	bytes	= (emax + 1)*sizeof(long);

	all = (long*)malloc( bytes );
	sam = (long*)malloc( bytes );
	dwn = (long*)malloc( bytes );

	memset( all, 0, bytes );
	memset( sam, 0, bytes );
	memset( dwn, 0, bytes );

	vector<float>	ve;

	for( int z = gArgs.zilo; z <= gArgs.zihi; ++z ) {

		char	buf[2048];
		FILE	*f;
		long	n;

		/* ---- */
		/* Same */
		/* ---- */

		sprintf( buf, "%s/Err_S_%d.bin", path, z );
		n = (long)DskBytes( buf ) / sizeof(float);

		if( !n )
			continue;

		ve.resize( n );

		f = FileOpenOrDie( buf, "rb" );
		fread( &ve[0], sizeof(float), n, f );
		fclose( f );

		for( int i = 0; i < n; ++i ) {

			double	e    = gArgs.div * ve[i];
			int		ibin = (e < emax ? int(e) : emax);

			++all[ibin];
			++sam[ibin];
		}

		/* ---- */
		/* Down */
		/* ---- */

		sprintf( buf, "%s/Err_D_%d.bin", path, z );
		n = (long)DskBytes( buf ) / sizeof(float);

		if( !n )
			continue;

		ve.resize( n );

		f = FileOpenOrDie( buf, "rb" );
		fread( &ve[0], sizeof(float), n, f );
		fclose( f );

		for( int i = 0; i < n; ++i ) {

			double	e    = gArgs.div * ve[i];
			int		ibin = (e < emax ? int(e) : emax);

			++all[ibin];
			++dwn[ibin];
		}
	}
}

/* --------------------------------------------------------------- */
/* Record -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Record()
{
	int	n = gArgs.lim * gArgs.div + 1;

	for( int i = 0; i < n; ++i ) {

		if( gArgs.inB ) {

			fprintf( flog,
			"%.2f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
			double(i + 1)/gArgs.div,
			A.all[i], A.sam[i], A.dwn[i],
			B.all[i], B.sam[i], B.dwn[i] );
		}
		else {
			fprintf( flog,
			"%.2f\t%ld\t%ld\t%ld\n",
			double(i + 1)/gArgs.div,
			A.all[i], A.sam[i], A.dwn[i] );
		}
	}
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

/* ---- */
/* Read */
/* ---- */

	A.Read( gArgs.inA );

	if( gArgs.inB )
		B.Read( gArgs.inB );

/* ----- */
/* Print */
/* ----- */

	Record();

/* ---- */
/* Done */
/* ---- */

	fclose( flog );

	return 0;
}


