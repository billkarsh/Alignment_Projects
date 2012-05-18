//
// Concatenate TrakEM2 file2 to file1, like XMLAppend, and,
// assuming the last layer of file1 is the same as the first
// layer of file2, adjust all transforms of file2 to map to
// the global space of file1.
//

#include	"Cmdline.h"
#include	"CRegexID.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"
#include	"CTForm.h"


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
	char	*infile1,
			*infile2;
	int		zolap,
			reftile;

public:
	CArgs_xml()
	{
		infile1	= NULL;
		infile2	= NULL;
		zolap	= 1;
		reftile	= 0;
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

	flog = FileOpenOrDie( "XMLUnite.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	char	*pat;

	re_id.Set( "_N_" );

	if( argc < 3 ) {
		printf( "Usage: XMLUnite <xml-file1> <xml-file2> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' ) {

			if( !infile1 )
				infile1 = argv[i];
			else
				infile2 = argv[i];
		}
		else if( GetArgStr( pat, "-p", argv[i] ) )
			re_id.Set( pat );
		else if( GetArg( &zolap, "-olap=%d", argv[i] ) )
			;
		else if( GetArg( &reftile, "-ref=%d", argv[i] ) )
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
		fprintf( flog, "No tile-id found in '%s'.\n", name );
		exit( 42 );
	}

	return id;
}

/* -------------------------------------------------------------- */
/* GetThisTForm ------------------------------------------------- */
/* -------------------------------------------------------------- */

static bool GetThisTForm( TForm &T, TiXmlElement* layer, int id )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		if( gArgs.IDFromPatch( ptch ) == id ) {

			T.ScanTrackEM2( ptch->Attribute( "transform" ) );
			return true;
		}
	}

	return false;
}

/* -------------------------------------------------------------- */
/* UpdateTiles -------------------------------------------------- */
/* -------------------------------------------------------------- */

static void UpdateTiles( TiXmlElement* layer, const TForm& T )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		char	buf[256];
		TForm	told, tnew;

		told.ScanTrackEM2( ptch->Attribute( "transform" ) );
		MultiplyTrans( tnew, T, told );

		sprintf( buf, "matrix(%f,%f,%f,%f,%f,%f)",
		tnew.t[0], tnew.t[3], tnew.t[1],
		tnew.t[4], tnew.t[2], tnew.t[5] );

		ptch->SetAttribute( "transform", buf );
	}
}

/* -------------------------------------------------------------- */
/* Unite -------------------------------------------------------- */
/* -------------------------------------------------------------- */

static void Unite()
{
/* -------------- */
/* Load documents */
/* -------------- */

	TiXmlDocument	doc1( gArgs.infile1 ),
					doc2( gArgs.infile2 );

	if( !doc1.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile1 );
		exit( 42 );
	}

	if( !doc2.LoadFile() ) {
		fprintf( flog,
		"Could not open XML file [%s].\n", gArgs.infile2 );
		exit( 42 );
	}

/* ---------------- */
/* Verify <trakem2> */
/* ---------------- */

	TiXmlHandle		hDoc1( &doc1 ),
					hDoc2( &doc2 );
	TiXmlNode*		lyrset1;
	TiXmlElement*	layer1;
	TiXmlElement*	layer2;

	if( !doc1.FirstChild() || !doc2.FirstChild() ) {
		fprintf( flog, "No trakEM2 node.\n" );
		exit( 42 );
	}

	lyrset1 = hDoc1.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.ToNode();

	layer2 = hDoc2.FirstChild( "trakem2" )
				.FirstChild( "t2_layer_set" )
				.FirstChild( "t2_layer" )
				.ToElement();

	if( !lyrset1 || !layer2 ) {
		fprintf( flog, "No first trakEM2 child.\n" );
		exit( 42 );
	}

/* -------------------------------- */
/* Advance layer2 to top of overlap */
/* -------------------------------- */

	for( int i = 1; i < gArgs.zolap; ++i )
		layer2 = layer2->NextSiblingElement();

	if( !layer2 ) {
		fprintf( flog, "File2 fully overlaps file1.\n" );
		exit( 42 );
	}

/* ----- */
/* Unite */
/* ----- */

	layer1 = lyrset1->LastChild( "t2_layer" )->ToElement();

// get last z to propagate to added layers

	int	z		= atoi( layer1->Attribute( "z" ) );
	int	nextoid	= NextOID( hDoc1 );

// get transform mapping T

	TForm	T;

	{
		TForm	T1, T2i, T2;

		if( !GetThisTForm( T1, layer1, gArgs.reftile ) ) {
			fprintf( flog, "Reference tile not found in file1.\n" );
			exit( 42 );
		}

		if( !GetThisTForm( T2, layer2, gArgs.reftile ) ) {
			fprintf( flog, "Reference tile not found in file2.\n" );
			exit( 42 );
		}

		InvertTrans( T2i, T2 );
		MultiplyTrans( T, T1, T2i );
	}

// update and add

	layer2 = layer2->NextSiblingElement();

	do {

		layer2->SetAttribute( "z", ++z );
		nextoid = SetOID( layer2, nextoid );
		UpdateTiles( layer2, T );

		lyrset1->InsertEndChild( *layer2 );

	} while( layer2 = layer2->NextSiblingElement() );

/* ---- */
/* Save */
/* ---- */

	doc1.SaveFile( "xmltmp.txt" );

/* ----------------- */
/* Copy !DOCTYPE tag */
/* ----------------- */

	CopyDTD( gArgs.infile1, "xmltmp.txt" );
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

/* ------------- */
/* Write new xml */
/* ------------- */

	Unite();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


