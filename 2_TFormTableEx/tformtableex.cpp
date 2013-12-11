//
// Merge {TAffineTable or THmgphyTable}
// and idb into new Bill file.
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	const char	*inpath,
				*idb;
	double		degcw;

public:
	CArgs()
	{
		inpath	= NULL;
		idb		= NULL;
		degcw	= 0.0;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;
static FILE*	flog = NULL;
static int		gW, gH;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "TFormTable.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: tformtable table idb\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !inpath )
				inpath = argv[i];
			else
				idb = argv[i];
		}
		else if( GetArg( &degcw, "-degcw=%lf", argv[i] ) )
			;
		else {
			printf( "Did not understand option [%s].\n", argv[i] );
			exit( 42 );
		}
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* MergeH -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetXYH(
	double			&x0,
	double			&y0,
	const THmgphy	&R,
	FILE			*fi )
{
	CLineScan	LS;

	x0 = BIGD;
	y0 = BIGD;

	while( LS.Get( fi ) > 0 ) {

		THmgphy	T;
		int		z, id, rgn;

		sscanf( LS.line, "%d\t%d\t%d"
		"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&z, &id, &rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5],
		&T.t[6], &T.t[7] );

		T = R * T;

		vector<Point>	cnr( 4 );

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			x0 = fmin( x0, cnr[k].x );
			y0 = fmin( y0, cnr[k].y );
		}
	}

	fseek( fi, 0, SEEK_SET );
}


static void MergeH()
{
	FILE		*fi = FileOpenOrDie( gArgs.inpath, "r" );
	FILE		*fo = FileOpenOrDie( "THmgphyTableEx.txt", "w" );
	CLineScan	LS;
	string		idb = gArgs.idb;
	THmgphy		R;
	double		x0 = 0, y0 = 0;

	if( gArgs.degcw ) {
		R.NUSetRot( gArgs.degcw * PI/180.0 );
		GetXYH( x0, y0, R, fi );
	}

	THmgphy	S( 1, 0, -x0, 0, 1, -y0, 0, 0 );

	while( LS.Get( fi ) > 0 ) {

		THmgphy	T;
		int		z, id, rgn;

		sscanf( LS.line, "%d\t%d\t%d"
		"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&z, &id, &rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5],
		&T.t[6], &T.t[7] );

		const Til2Img	*t2i;

		if( !IDBT2ICacheNGet1( t2i, idb, z, id, flog ) )
			continue;

		if( gArgs.degcw )
			T = S * (R * T);

		fprintf( fo, "%d\t%d"
		"\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
		"\t%d\t%d\t%d\t%s\n",
		z, id,
		T.t[0], T.t[1], T.t[2],
		T.t[3], T.t[4], T.t[5],
		T.t[6], T.t[7],
		t2i->col, t2i->row, t2i->cam, t2i->path.c_str() );
	}

	fclose( fo );
	fclose( fi );
}

/* --------------------------------------------------------------- */
/* MergeA -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetXYA(
	double			&x0,
	double			&y0,
	const TAffine	&R,
	FILE			*fi )
{
	CLineScan	LS;

	x0 = BIGD;
	y0 = BIGD;

	while( LS.Get( fi ) > 0 ) {

		TAffine	T;
		int		z, id, rgn;

		sscanf( LS.line, "%d\t%d\t%d"
		"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&z, &id, &rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5] );

		T = R * T;

		vector<Point>	cnr( 4 );

		cnr[0] = Point(  0.0, 0.0 );
		cnr[1] = Point( gW-1, 0.0 );
		cnr[2] = Point( gW-1, gH-1 );
		cnr[3] = Point(  0.0, gH-1 );

		T.Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			x0 = fmin( x0, cnr[k].x );
			y0 = fmin( y0, cnr[k].y );
		}
	}

	fseek( fi, 0, SEEK_SET );
}


static void MergeA()
{
	FILE		*fi = FileOpenOrDie( gArgs.inpath, "r" );
	FILE		*fo = FileOpenOrDie( "TAffineTableEx.txt", "w" );
	CLineScan	LS;
	string		idb = gArgs.idb;
	TAffine		R;
	double		x0, y0;

	if( gArgs.degcw ) {
		R.NUSetRot( gArgs.degcw * PI/180.0 );
		GetXYA( x0, y0, R, fi );
	}

	while( LS.Get( fi ) > 0 ) {

		TAffine	T;
		int		z, id, rgn;

		sscanf( LS.line, "%d\t%d\t%d"
		"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&z, &id, &rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5] );

		const Til2Img	*t2i;

		if( !IDBT2ICacheNGet1( t2i, idb, z, id, flog ) )
			continue;

		if( gArgs.degcw ) {
			T = R * T;
			T.t[2] -= x0;
			T.t[5] -= y0;
		}

		fprintf( fo, "%d\t%d"
		"\t%f\t%f\t%f\t%f\t%f\t%f"
		"\t%d\t%d\t%d\t%s\n",
		z, id,
		T.t[0], T.t[1], T.t[2], T.t[3], T.t[4], T.t[5],
		t2i->col, t2i->row, t2i->cam, t2i->path.c_str() );
	}

	fclose( fo );
	fclose( fi );
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

/* ----- */
/* Merge */
/* ----- */

	if( !IDBGetImageDims( gW, gH, gArgs.idb, flog ) )
		exit( 42 );

	if( FileNamePtr( gArgs.inpath )[1] == 'H' )
		MergeH();
	else
		MergeA();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


