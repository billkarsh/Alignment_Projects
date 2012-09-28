//
// Collect scapeops results into low res stack for viewing/editing.
//


#include	"Cmdline.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"../1_Cross_LowRes/ScapeMeta.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_alnmon -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_alnmon {

public:
	int		zmin,
			zmax;

public:
	CArgs_alnmon()
	{
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_alnmon	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_alnmon::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "cross_lowres.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Assemble stack: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf(
		"Usage: cross_lowres -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
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
/* MakeTForms ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeTForms(
	vector<TForm>		&vT,
	const vector<CLog>	&vL )
{
	int	nL = vT.size();

// T[0] is just indentity, defining global coords.
// Create T[1]...T[n] taking image[j] -> image[0].

	for( int ia = 1; ia < nL; ++ia ) {

		TForm	Taa;
		int		ib = ia - 1;

		// Begin by constructing Tba (image[a] -> image[b])

		// A image -> A content
		Taa.AddXY( vL[ia].M.x0, vL[ia].M.y0 );

		// Image size -> strip size
		Taa.MulXY( (double)vL[ia].M.scl / vL[ib].A.scl );

		// A content -> A strip
		Taa.AddXY( -vL[ib].A.x0, -vL[ib].A.y0 );

		// A strip -> B strip
		MultiplyTrans( vT[ia], vL[ib].T, Taa );

		// B strip -> B content
		vT[ia].AddXY( vL[ib].B.x0, vL[ib].B.y0 );

		// Strip size -> image size
		vT[ia].MulXY( (double)vL[ib].B.scl / vL[ia].M.scl );

		// B content -> B image
		vT[ia].AddXY( -vL[ib].M.x0, -vL[ib].M.y0 );

		// Now convert Tba to global transform T0a
		// T0a = T01.T12...Tba.
		MultiplyTrans( vT[ia], vT[ib], TForm( vT[ia] ) );
	}
}

/* --------------------------------------------------------------- */
/* MakeBounds ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeBounds(
	DBox				&B,
	const vector<CLog>	&vL,
	const vector<TForm>	&vT )
{
	int	nL = vT.size();

	B.L = BIGD, B.R = -BIGD,
	B.B = BIGD, B.T = -BIGD;

	for( int i = 0; i < nL; ++i ) {

		vector<Point>	cnr( 4 );

		cnr[0] = Point(            0.0,            0.0 );
		cnr[1] = Point( vL[i].M.ws - 1,            0.0 );
		cnr[2] = Point( vL[i].M.ws - 1, vL[i].M.hs - 1 );
		cnr[3] = Point(            0.0, vL[i].M.hs - 1 );

		vT[i].Transform( cnr );

		for( int k = 0; k < 4; ++k ) {

			B.L = fmin( B.L, cnr[k].x );
			B.R = fmax( B.R, cnr[k].x );
			B.B = fmin( B.B, cnr[k].y );
			B.T = fmax( B.T, cnr[k].y );
		}
	}
}

/* --------------------------------------------------------------- */
/* WriteTrakEM2Layer --------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTrakEM2Layer(
	FILE*		f,
	int			&oid,
	const CLog	&L,
	const TForm	&T )
{
// Layer prologue

	fprintf( f,
	"\t\t<t2_layer\n"
	"\t\t\toid=\"%d\"\n"
	"\t\t\tthickness=\"0\"\n"
	"\t\t\tz=\"%d\"\n"
	"\t\t>\n",
	oid++, L.M.z );

// Tile - just tile 0

	fprintf( f,
	"\t\t\t<t2_patch\n"
	"\t\t\t\toid=\"%d\"\n"
	"\t\t\t\twidth=\"%d\"\n"
	"\t\t\t\theight=\"%d\"\n"
	"\t\t\t\ttransform=\"matrix(%f,%f,%f,%f,%f,%f)\"\n"
	"\t\t\t\ttitle=\"M_%d_0\"\n"
	"\t\t\t\ttype=\"0\"\n"
	"\t\t\t\tfile_path=\"montages/M_%d_0.png\"\n"
	"\t\t\t\to_width=\"%d\"\n"
	"\t\t\t\to_height=\"%d\"\n"
	"\t\t\t\tmin=\"0\"\n"
	"\t\t\t\tmax=\"255\"\n"
	"\t\t\t/>\n",
	oid++, L.M.ws, L.M.hs,
	T.t[0], T.t[3], T.t[1], T.t[4], T.t[2], T.t[5],
	L.M.z, L.M.z, L.M.ws, L.M.hs );

// Layer epilogue

	fprintf( f, "\t\t</t2_layer>\n" );
}

/* --------------------------------------------------------------- */
/* WriteTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteTrakEM2(
	const char			*path,
	DBox				&B,
	const vector<CLog>	&vL,
	const vector<TForm>	&vT )
{
// Open file

	FILE	*f = FileOpenOrDie( path, "w" );

// Prologue + bounds

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
	oid++, B.R - B.L, B.T - B.B );

// Layers

	int	nL = vT.size();

	for( int i = 0; i < nL; ++i )
		WriteTrakEM2Layer( f, oid, vL[i], vT[i] );

// Epilogue

	fprintf( f, "\t</t2_layer_set>\n" );
	fprintf( f, "</trakem2>\n" );
	fclose( f );
}

/* --------------------------------------------------------------- */
/* BuildStack ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void BuildStack()
{
// Get log data

	vector<CLog>	vL;
	int				nL = ReadLogs( vL, gArgs.zmin, gArgs.zmax );

	if( !nL )
		return;

// Make layer TForms

	vector<TForm>	vT( nL );

	MakeTForms( vT, vL );

// Calculate bounds

	DBox	B;

	MakeBounds( B, vL, vT );

// Write

	WriteTrakEM2( "LowRes.xml", B, vL, vT );
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

/* ----------- */
/* Build stack */
/* ----------- */

	BuildStack();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


