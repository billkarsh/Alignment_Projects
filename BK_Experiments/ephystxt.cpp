

#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"GenDefs.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* CArgs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs {
private:
	const char	*meta;
private:
	char	path[2048];
	char	*basename;
public:
	CArgs() : meta(NULL) {};
	void SetCmdLine( int argc, char* argv[] );
	void ParseName();
	FILE* OpenMeta()
	{
		char	buf[2048];
		sprintf( buf, "%s/%s.meta", path, basename );
		return FileOpenOrDie( buf, "r" );
	};
	void BinPath( char *buf )
	{
		sprintf( buf, "%s/%s.bin", path, basename );
	};
	double BinSize()
	{
		char	buf[2048];
		BinPath( buf );
		return DskBytes( buf );
	}
	FILE* OpenBin()
	{
		char	buf[2048];
		BinPath( buf );
		return FileOpenOrDie( buf, "rb" );
	};
	FILE* OpenText()
	{
		char	buf[2048];
		sprintf( buf, "%s.txt", basename );
		return FileOpenOrDie( buf, "w" );
	};
};

class Meta {
public:
	int	Hz,
		nchan;
public:
	Meta() : Hz(0), nchan(0) {};
	bool Scan();
	void Convert();
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs	gArgs;






/* --------------------------------------------------------------- */
/* CArgs::SetCmdLine --------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::SetCmdLine( int argc, char* argv[] )
{
// Parse command line args

	if( argc < 2 ) {
		printf( "Usage: ephystxt <meta-file> [options].\n" );
		exit( 42 );
	}

	for( int i = 1; i < argc; ++i ) {

		if( argv[i][0] != '-' )
			meta = argv[i];
		else {
			printf( "Did not understand option '%s'.\n", argv[i] );
			exit( 42 );
		}
	}
}

/* --------------------------------------------------------------- */
/* CArgs::ParseName ---------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs::ParseName()
{
	strcpy( path, meta );

	char *s = strrchr( path, '/' );

	if( s )
		*s = 0;
	else {
		path[0] = '.';
		path[1] = 0;
	}

	basename = FileCloneNamePart( meta );
}

/* --------------------------------------------------------------- */
/* Meta::Scan ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool Meta::Scan()
{
	FILE		*f = gArgs.OpenMeta();
	CLineScan	LS;

	while( LS.Get( f ) > 0 ) {

		if( LS.line[0] == 'n' ) {

			if( 1 == sscanf( LS.line, "nChans = %d", &nchan ) )
				printf( "num channels: %d\n", nchan );
		}

		if( LS.line[0] == 's' ) {

			if( 1 == sscanf( LS.line, "sRateHz = %d", &Hz ) )
				printf( "Hertz: %d\n", Hz );
		}
	}

	fclose( f );

	bool gotall = true;

	if( !Hz ) {
		printf( "Meta: Missing sRateHz field.\n" );
		gotall = false;
	}

	if( !nchan ) {
		printf( "Meta: Missing nChans field.\n" );
		gotall = false;
	}

	return gotall;
}

/* --------------------------------------------------------------- */
/* Meta::Convert ------------------------------------------------- */
/* --------------------------------------------------------------- */

void Meta::Convert()
{
	const int maxsamples = 20000;

	double	filebytes = gArgs.BinSize();
	int		sampbytes = sizeof(int16) * nchan,
			nsamp = int(filebytes / sampbytes),
			readbytes;

	if( nsamp > maxsamples )
		nsamp = maxsamples;

	vector<char>	D( readbytes = nsamp * sampbytes );
	FILE			*f;

	f = gArgs.OpenBin();
	fread( &D[0], 1, readbytes, f );
	fclose( f );

	f = gArgs.OpenText();

	for( int is = 0; is < nsamp; ++is ) {

		fprintf( f, "%.4e", double(is)/Hz );

		char*	s0 = &D[is * sampbytes];

		for( int ic = 0; ic < nchan; ++ic ) {
			fprintf( f, "\t%d", *(int16*)(s0 + sizeof(int16)*ic) );
		}

		fprintf( f, "\n" );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
	gArgs.SetCmdLine( argc, argv );
	gArgs.ParseName();

	Meta	M;

	if( M.Scan() )
		M.Convert();

	return 0;
}


