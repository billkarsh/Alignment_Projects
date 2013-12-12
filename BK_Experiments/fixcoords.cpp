//
// Update a Rick file by replacing the XY coords with those
// read from an ImageMeta.txt file.
//

#include	"Cmdline.h"
#include	"File.h"

#include	<string.h>

#include	<map>
#include	<string>
#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CTile {
public:
	string	name;
	double	x,
			y;
	int		z;
};


class CXY {
public:
	double	x, y;
};

/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {

public:
	double	umperpix;
	char	*rick,
			*meta;

public:
	CArgs()
	{
		umperpix	= 0.645;
		rick		= NULL;
		meta		= NULL;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs			gArgs;
static FILE*			flog	= NULL;
static vector<CTile>	vT;
static map<string,CXY>	M;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "fixcoords.log", "w" );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: fixcoords rickfile metafile [-umperpix=].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' ) {

			if( !rick )
				rick = argv[i];
			else
				meta = argv[i];
		}
		else if( GetArg( &umperpix, "-umperpix=%lf", argv[i] ) )
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
/* ReadRick ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ReadRick()
{
	FILE	*f = FileOpenOrDie( gArgs.rick, "r", flog );

/* ---------- */
/* Scan lines */
/* ---------- */

	for( ;; ) {

		CTile	T;
		char	name[2048];

		/* ---------- */
		/* Get a line */
		/* ---------- */

		if( fscanf( f, "%s%lf%lf%d\n", name, &T.x, &T.y, &T.z ) != 4 )
			break;

		/* --------- */
		/* Set entry */
		/* --------- */

		T.name = name;

		vT.push_back( T );
	}

/* ----- */
/* Close */
/* ----- */

	fclose( f );
}

/* --------------------------------------------------------------- */
/* ReadMeta ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void ReadMeta()
{
	FILE	*f = FileOpenOrDie( gArgs.meta, "r", flog );

/* ----------- */
/* Skip header */
/* ----------- */

	char	*buf;
	size_t	dum;

	getline( &buf, &dum, f );

/* ---------- */
/* Scan lines */
/* ---------- */

	for( ;; ) {

		CXY		xy;
		double	x, y;
		char	name[2048];

		/* ---------- */
		/* Get a line */
		/* ---------- */

		if( fscanf( f,
			"%*d%*d%*d"
			"%lf%lf"
			"%*lf%*lf%*lf%*lf%*lf%*lf%*lf"
			"%*d%*d"
			"\t%s",
			&x, &y, name ) != 3 ) {

			break;
		}

		/* --------- */
		/* Set entry */
		/* --------- */

		const char	*fname = strrchr( name, '\\' );

		if( fname )
			++fname;
		else
			fname = name;

		string	s = fname;

		xy.x = -y / gArgs.umperpix;
		xy.y = -x / gArgs.umperpix;

		M[s] = xy;
	}

/* ----- */
/* Close */
/* ----- */

	fclose( f );
}

/* --------------------------------------------------------------- */
/* Update -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Update()
{
	int	nt = vT.size();

	for( int i = 0; i < nt; ++i ) {

		CTile	&T = vT[i];

		map<string,CXY>::iterator	it = M.find( T.name );

		if( it != M.end() ) {
			T.x = it->second.x;
			T.y = it->second.y;
		}
	}

	FILE	*f = FileOpenOrDie( "newrick.txt", "w", flog );

	for( int i = 0; i < nt; ++i ) {

		const CTile	&T = vT[i];

		fprintf( f, "%s\t%f\t%f\t%d\n",
			T.name.c_str(), T.x, T.y, T.z );
	}

	fclose( f );
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

	ReadRick();
	ReadMeta();

	Update();

/* ---- */
/* Done */
/* ---- */

exit:
	fclose( flog );

	return 0;
}


