//
// Get one binary histogram file from
// one 16-bit gray image.
//

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"

#include	<string.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_hist1 --------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_hist1 {
public:
	char	*img, *hst;
public:
	CArgs_hist1() : img(NULL), hst(NULL) {};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_hist1	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_hist1::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "Hist1.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: Hist1 <img-file> <hst_file>.\n" );
		exit( 42 );
	}

	img = argv[1];
	hst = argv[2];

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );
	}

	fprintf( flog, "\n\n" );
	fflush( flog );
}

/* --------------------------------------------------------------- */
/* Hist ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Hist()
{
// gather histogram

	const int		nbins = 65536;	// also max val
	vector<double>	bins( nbins, 0.0 );
	double			uflo = 0.0,
					oflo = 0.0;
	uint32			w, h;
	uint16*			ras;

	ras = Raster16FromTif16( gArgs.img, w, h, flog );

	Histogram( uflo, oflo, &bins[0], nbins,
		0.0, nbins, ras, w * h, false );

	RasterFree( ras );

// write binary hist file

	FILE	*f;

	if( f = fopen( gArgs.hst, "wb" ) ) {
		fwrite( &uflo, sizeof(double), 1, f );
		fwrite( &oflo, sizeof(double), 1, f );
		fwrite( &bins[0], sizeof(double), nbins, f );
		fclose( f );
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

/* --------- */
/* Histogram */
/* --------- */

	Hist();

/* ---- */
/* Done */
/* ---- */

	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


