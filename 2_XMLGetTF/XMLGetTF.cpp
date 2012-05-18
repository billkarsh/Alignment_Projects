//
// Get TFormTable in given z range.
//

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"CTForm.h"

#include	"tinyxml.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_xml ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_xml {

private:
	// re_id used to extract tile id from image name.
	// "/N" used for EM projects, "_N_" for APIG images,
	// "_Nex.mrc" typical for Leginon files.
	CRegexID	re_id;

public:
	char	*infile;
	int		zmin, zmax;

public:
	CArgs_xml()
	{
		infile	= NULL;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );

	int IDFromPatch( TiXmlElement* p );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_xml	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "XMLGetTF.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	char	*pat;

	re_id.Set( "_N_" );

	if( argc < 4 ) {
		printf( "Usage: XMLGetTF <xml-file1> -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
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

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int CArgs_xml::IDFromPatch( TiXmlElement* p )
{
	const char	*name = p->Attribute( "title" );
	int			id;

	if( !re_id.Decode( id, name ) ) {
		printf( "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
}

/* --------------------------------------------------------------- */
/* GetSortedTForms ----------------------------------------------- */
/* --------------------------------------------------------------- */

class TS {
// Use for sorting TForms
public:
	TForm	T;
	int		id;
public:
	TS( int _i, TForm &_T )	{id=_i; T=_T;};

	bool operator < (const TS &rhs) const
		{return id < rhs.id;};
};


static void GetSortedTForms(
	vector<TS>		&ts,
	TiXmlElement*	layer )
{
	TiXmlElement*	p = layer->FirstChildElement( "t2_patch" );

	for( ; p; p = p->NextSiblingElement() ) {

		TForm	T;
		int		id = gArgs.IDFromPatch( p );

		T.ScanTrackEM2( p->Attribute( "transform" ) );
		ts.push_back( TS( id, T ) );
	}

	sort( ts.begin(), ts.end() );
}

/* --------------------------------------------------------------- */
/* PrintTForms --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void PrintTForms( FILE *f, const vector<TS> &ts, int z )
{
	int	nt = ts.size();

	for( int i = 0; i < nt; ++i ) {

		const double *t = ts[i].T.t;

		fprintf( f, "%d\t%d\t1\t%f\t%f\t%f\t%f\t%f\t%f\n",
		z, ts[i].id, t[0], t[1], t[2], t[3], t[4], t[5] );
	}
}

/* --------------------------------------------------------------- */
/* GetTF --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void GetTF()
{
/* ------------- */
/* Load document */
/* ------------- */

	TiXmlDocument	doc( gArgs.infile );

	if( !doc.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hdoc( &doc );
	TiXmlElement*	layer;

	if( !doc.FirstChild() ) {
		fprintf( flog, "No trakEM2 node [%s].\n", gArgs.infile );
		exit( 42 );
	}

	layer = hdoc.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !layer ) {
		fprintf( flog, "No t2_layer [%s].\n", gArgs.infile );
		exit( 42 );
	}

/* --- */
/* Get */
/* --- */

	char		buf[256];
	const char	*st = FileNamePtr( gArgs.infile );

	sprintf( buf, "%.*s_TF.txt", strlen( st ) - 4, st );

	FILE	*f = FileOpenOrDie( buf, "w", flog );

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		vector<TS>	ts;

		GetSortedTForms( ts, layer );
		PrintTForms( f, ts, z );
	}

	fclose( f );
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
/* Process */
/* ------- */

	GetTF();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


