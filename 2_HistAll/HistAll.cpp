//
// Write script to submit Hist1 jobs for given:
//
// z range
// channel list
//

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<sys/stat.h>


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class Picture {

public:
	string	fname;	// name excluding '_chan.tif'
	int		z, id;

public:
	Picture( const char* name, int _z );

	bool operator < (const Picture &rhs) const
		{return z < rhs.z || (z == rhs.z && id < rhs.id);};
};

/* --------------------------------------------------------------- */
/* CArgs_hsta ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_hsta {

public:
	char		*infile;
	int			zmin,
				zmax;
	vector<int>	chn;

public:
	CArgs_hsta()
	{
		infile	= NULL;
		zmin	= 0;
		zmax	= 32768;
	};

	void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_hsta	gArgs;
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_hsta::SetCmdLine( int argc, char* argv[] )
{
// start log

	flog = FileOpenOrDie( "HistAll.log", "w" );

// log start time

	time_t	t0 = time( NULL );
	char	atime[32];

	strcpy( atime, ctime( &t0 ) );
	atime[24] = '\0';	// remove the newline

	fprintf( flog, "Start: %s ", atime );

// parse command line args

	if( argc < 3 ) {
		printf( "Usage: HistAll <xml-file> -chan=,, [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		// echo to log
		fprintf( flog, "%s ", argv[i] );

		if( argv[i][0] != '-' )
			infile = argv[i];
		else if( GetArg( &zmin, "-zmin=%d", argv[i] ) )
			;
		else if( GetArg( &zmax, "-zmax=%d", argv[i] ) )
			;
		else if( GetArgList( chn, "-chan=", argv[i] ) )
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
/* Picture ------------------------------------------------------- */
/* --------------------------------------------------------------- */

Picture::Picture( const char* name, int _z )
{
	char	buf[2048];
	char	*s;
	int		len;

	len	= sprintf( buf, "%s", name );
	s	= buf + len - 6;
	*s	= 0;	// remove '_chn.tif'

	while( isdigit( s[-1] ) )
		--s;

	fname	= buf;
	z		= _z;
	id		= atoi( s );
}

/* --------------------------------------------------------------- */
/* GetTiles ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetTiles(
	vector<Picture>	&vp,
	TiXmlElement*	layer,
	int				z )
{
	TiXmlElement*	ptch = layer->FirstChildElement( "t2_patch" );

	for( ; ptch; ptch = ptch->NextSiblingElement() ) {

		vp.push_back( Picture( ptch->Attribute( "file_path" ), z ) );
	}
}

/* --------------------------------------------------------------- */
/* ParseTrakEM2 -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ParseTrakEM2( vector<Picture> &vp )
{
/* ---- */
/* Open */
/* ---- */

	XML_TKEM		xml( gArgs.infile, flog );
	TiXmlElement*	layer	= xml.GetFirstLayer();

/* -------------- */
/* For each layer */
/* -------------- */

	for( ; layer; layer = layer->NextSiblingElement() ) {

		/* ----------------- */
		/* Layer-level stuff */
		/* ----------------- */

		int	z = atoi( layer->Attribute( "z" ) );

		if( z > gArgs.zmax )
			break;

		if( z < gArgs.zmin )
			continue;

		GetTiles( vp, layer, z );
	}
}

/* --------------------------------------------------------------- */
/* WriteScript --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScript( const vector<Picture> &vp )
{
// HST directory

	DskCreateDir( "HST", flog );

// open file

	FILE	*f = FileOpenOrDie( "HST/make.hst.sh", "w", flog );

// write

	int		np = vp.size(),
			nc = gArgs.chn.size();

	fprintf( f, "#!/bin/sh\n\n" );

	for( int ip = 0; ip < np; ++ip ) {

		const Picture&	P = vp[ip];

		for( int ic = 0; ic < nc; ++ic ) {

			int	chn = gArgs.chn[ic];

			fprintf( f,
			"qsub -N hst-%d-%d -j y -o out.txt -b y -cwd -V"
			" Hist1 '%s_%d.tif' 'HST_%d_%d_%d.bin'\n",
			ip*nc+ic, np*nc,
			P.fname.c_str(), chn, P.z, P.id, chn );
		}
	}

	fprintf( f, "\n" );

	fclose( f );

// make executable

	chmod( "HST/make.hst.sh", S_IRWXU | S_IRWXG | S_IRWXO );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	vector<Picture>	vp;

/* ------------------ */
/* Parse command line */
/* ------------------ */

	gArgs.SetCmdLine( argc, argv );

/* ---------------- */
/* Read source file */
/* ---------------- */

	ParseTrakEM2( vp );

	fprintf( flog, "Got %d images.\n", vp.size() );

	if( !vp.size() )
		goto exit;

/* ------------------ */
/* Sort by Z and tile */
/* ------------------ */

	sort( vp.begin(), vp.end() );

/* ------------ */
/* Write script */
/* ------------ */

	WriteScript( vp );

/* ---- */
/* Done */
/* ---- */

exit:
	fprintf( flog, "\n" );
	fclose( flog );

	return 0;
}


