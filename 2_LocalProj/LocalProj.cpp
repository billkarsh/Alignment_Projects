//
// Create portable xml project from given xml file...
//
// E.g.
// > localproj aaa.xml -dst=/groups/projects/myproj -zmin=10 -zmax=15
//
// Result:
// 1) Create new folder: '/groups/projects/myproj'
// 2) Copy aaa.xml to local 'myproj/aaa.xml'
// 3) Create new folder: 'myproj/Images' to hold all local images
// 4) Create new folder: 'Images/Znnn' for each layer
// 5) Copy images to the new Znnn folders
// 6) Edit local aaa.xml to remove all layers outside [zmin.zmax]
// 6) Edit local aaa.xml to point to local images (relative paths)
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"


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
	char	locxml[2048],
			images[2048],
			curlyr[2048];

public:
	char	*infile,
			*dst;
	int		zmin, zmax;

public:
	CArgs_xml()
	{
		locxml[0]	= 0;
		images[0]	= 0;
		curlyr[0]	= 0;
		infile		= NULL;
		dst			= NULL;
		zmin		= 0;
		zmax		= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
	void TopLevel();
	void NewLayer( int z );
	int  CopyImg( TiXmlElement* ptch );
	void SetRelPath( TiXmlElement* ptch );
	char *GetXML()	{return locxml;};
	void RenameXML();
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

	flog = FileOpenOrDie( "LocalProj.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 5 ) {
		printf( "Usage: localproj <xml-file1>"
		" -dst=path -zmin=i -zmax=j.\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArgStr( dst, "-dst=", argv[i] ) )
			;
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
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
/* TopLevel ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CArgs_xml::TopLevel()
{
	char		buf[2048];
	const char	*s;

// create top directory
	DskCreateDir( dst, flog );

// copy xml
	s = FileNamePtr( infile );
	sprintf( locxml, "%s/%s", dst, s );
	sprintf( buf, "cp %s %s", infile, locxml );
	system( buf );

// create images folder
	sprintf( images, "%s/Images", dst );
	DskCreateDir( images, flog );
}

/* --------------------------------------------------------------- */
/* NewLayer ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void CArgs_xml::NewLayer( int z )
{
	sprintf( curlyr, "%s/Z%d", images, z );
	DskCreateDir( curlyr, flog );
}

/* --------------------------------------------------------------- */
/* CopyImg ------------------------------------------------------- */
/* --------------------------------------------------------------- */

int CArgs_xml::CopyImg( TiXmlElement* ptch )
{
	char	buf[4096];
	const char	*p = ptch->Attribute( "file_path" ),
				*n = FileNamePtr( p );

	sprintf( buf, "cp %s %s/%s", p, curlyr, n );
	system( buf );

	sprintf( buf, "%s/%s", curlyr, n );

	return DskExists( buf );
}

/* --------------------------------------------------------------- */
/* SetRelPath ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::SetRelPath( TiXmlElement* ptch )
{
	char	buf[2048];

	sprintf( buf, "%s/%s",
		strstr( curlyr, "Images/Z" ),
		FileNamePtr( ptch->Attribute( "file_path" ) ) );

	ptch->SetAttribute( "file_path", buf );
}

/* --------------------------------------------------------------- */
/* RenameXML ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_xml::RenameXML()
{
	char	buf[4096];
	int		len = strlen( locxml );

	sprintf( buf, "mv %.*s_v2.xml %s", len - 4, locxml, locxml );
	system( buf );
}

/* --------------------------------------------------------------- */
/* UpdateTiles --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpdateTiles( TiXmlElement* layer )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		if( gArgs.CopyImg( ptch ) )
			gArgs.SetRelPath( ptch );
	}
}

/* --------------------------------------------------------------- */
/* DoLayers ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void DoLayers()
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.GetXML(), flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* ------------------------- */
/* Kill layers outside range */
/* ------------------------- */

	TiXmlNode*		lyrset = layer->Parent();
	TiXmlElement*	next;

	for( ; layer; layer = next ) {

		// next layer0 before deleting anything
		next = layer->NextSiblingElement();

		int	z = atoi( layer->Attribute( "z" ) );

		if( z < gArgs.zmin || z > gArgs.zmax )
			lyrset->RemoveChild( layer );
	}

/* --------------------------- */
/* Copies for remaining layers */
/* --------------------------- */

	layer = lyrset->FirstChild( "t2_layer" )->ToElement();

	for( ; layer; layer = layer->NextSiblingElement() ) {

		int	z = atoi( layer->Attribute( "z" ) );

		gArgs.NewLayer( z );
		UpdateTiles( layer );
	}

/* ---- */
/* Save */
/* ---- */

	xml.Save( "xmltmp.txt", true );

/* ------------------ */
/* Rename version two */
/* ------------------ */

	gArgs.RenameXML();
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
/* Top level */
/* --------- */

	gArgs.TopLevel();

/* ------------------ */
/* Process the layers */
/* ------------------ */

	DoLayers();

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


