//
// Collect scapeops results into rough stack.
//


#include	"GenDefs.h"
#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"CTForm.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CScapeMeta {
public:
	DBox	B;
	double	x0, y0;
	int		z,
			deg,
			scl,
			ws, hs;
};

class CLog {
public:
	CScapeMeta	M, A, B;
	TForm		T;
};

/* --------------------------------------------------------------- */
/* CArgs_alnmon -------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_alnmon {

public:
	char	outdir[2048];
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

	flog = FileOpenOrDie( "alignmontages2.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Assemble stack: %s ", atime );

// parse command line args

	if( argc < 4 ) {
		printf(
		"Usage: alignmontages2 -d. -zmin=i -zmax=j"
		" [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		char	*_outdir;

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( GetArgStr( _outdir, "-d", argv[i] ) )
			DskAbsPath( outdir, sizeof(outdir), _outdir, flog );
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
/* ReadLog1 ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool ReadLog1( CLog &L, int z )
{
	char		buf[2048];
	sprintf( buf, "%s/scplogs/scp_%d.log", gArgs.outdir, z );

	if( !DskExists( buf ) )
		return false;

	FILE		*f = FileOpenOrDie( buf, "r" );
	CLineScan	LS;

	while( LS.Get( f ) > 0 ) {

		if( LS.line[0] != '*' )
			continue;

		char	key = LS.line[1];

		LS.Get( f );

		if( key == 'M' ) {

			CScapeMeta	&E = L.M;

			sscanf( LS.line,
			"%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
			&E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
			&E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );
		}
		else if( key == 'A' ) {

			CScapeMeta	&E = L.A;

			sscanf( LS.line,
			"%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
			&E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
			&E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );
		}
		else if( key == 'B' ) {

			CScapeMeta	&E = L.B;

			sscanf( LS.line,
			"%d %d [%lf,%lf,%lf,%lf] %d [%d,%d] [%lf,%lf]",
			&E.z, &E.deg, &E.B.L, &E.B.R, &E.B.B, &E.B.T,
			&E.scl, &E.ws, &E.hs, &E.x0, &E.y0 );
		}
		else if( key == 'T' ) {

			double	*t = L.T.t;

			sscanf( LS.line + 1,
			"%lf,%lf,%lf,%lf,%lf,%lf",
			&t[0], &t[1], &t[2],
			&t[3], &t[4], &t[5] );

			break;
		}
	}

	fclose( f );
	return true;
}

/* --------------------------------------------------------------- */
/* ReadLogs ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static int ReadLogs( vector<CLog> &vL )
{
	for( int z = gArgs.zmin; z <= gArgs.zmax; ) {

		CLog	L;

		L.A.z = -1;

		if( ReadLog1( L, z ) ) {

			vL.push_back( L );

			if( L.A.z == -1 )
				break;

			z = L.A.z;
		}
		else
			++z;
	}

	return vL.size();
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
		MultiplyTrans( vT[ia], vT[ib], TForm(vT[ia]) );
	}
}

/* --------------------------------------------------------------- */
/* MakeBounds ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeBounds(
	DBox	&B,
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
	oid++, B.R, B.T );

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
	int				nL = ReadLogs( vL );

	if( !nL )
		return;

// Make layer TForms

	vector<TForm>	vT( nL );

	MakeTForms( vT, vL );

// Calculate bounds

	DBox	B;

	MakeBounds( B, vL, vT );

// Write

	WriteTrakEM2( "LowResMons.xml", B, vL, vT );
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


