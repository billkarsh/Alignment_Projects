

#include	"File.h"

#include	<string.h>

#include	<vector>
using namespace std;


static FILE	*flog = NULL;


class CRec {
public:
	char	name[32];
	double	x, y;
	int		z;
};






static void Fix( char *path, int len, int chn )
{
	path[len - 5] = '0' + chn;

	FILE	*f = fopen( path, "r" );

// exists?

	if( !f ) {
		fprintf( flog, "%d: nofile\n", chn );
		return;
	}

// load records

	vector<CRec>	V;

	for(;;) {

		CRec	C;
		char	name[2048];

		if( 4 != fscanf( f, "%s%lf%lf%d\n", name, &C.x, &C.y, &C.z ) )
			break;

		strcpy( C.name, name );
		V.push_back( C );
	}

	fclose( f );

// any records?

	int	nc = V.size();

	if( !nc ) {
		fprintf( flog, "%d: empty\n", chn );
		return;
	}

// backup as file.old

	sprintf( path + len - 3, "old" );
	f = FileOpenOrDie( path, "w" );
	for( int i = 0; i < nc; ++i ) {
		CRec&	C = V[i];
		fprintf( f, "%s\t%g\t%g\t%d\n",
			C.name, C.x, C.y, C.z );
	}
	fclose( f );

// write new

	const double	Q = 1.02 * 1.02 * 4.0;

	sprintf( path + len - 3, "txt" );
	f = FileOpenOrDie( path, "w" );
	for( int i = 0; i < nc; ++i ) {
		CRec&	C = V[i];
		fprintf( f, "%s\t%g\t%g\t%d\n",
			C.name, C.x/Q, C.y/Q, C.z );
	}
	fclose( f );

	fprintf( flog, "%d: --ok--\n", chn );
}


int main( int argc, char **argv )
{
	flog = FileOpenOrDie( "~~~temfixlog.txt", "w" );

	FILE	*fin = FileOpenOrDie( "runlist.txt", "r" );

	CLineScan	LS;

	for(;;) {

		// get line
		if( LS.Get( fin ) <= 0 )
			break;

		// remove trailing slashes
		char	path[2048];
		sscanf( LS.line, "%[^\r\n]", path );
		int	len = strlen( path );
		if( path[len - 1] == '/' )
			path[--len] = 0;
		fprintf( flog, "path=[%s]\n", path );

		// get runname
		char	run[64];
		char	*s = strrchr( path, '/' );
		if( !s )
			s = path;
		else
			++s;
		strcpy( run, s );
		fprintf( flog, "run=[%s]\n", run );

		// build fullpath
		len += sprintf( path + len, "/%s_TrackEM2_ChnX.txt", run );
		fprintf( flog, "base=[%s]\n", path );

		for( int chn = 0; chn < 4; ++chn )
			Fix( path, len, chn );

		// new paragraph
		fprintf( flog, "\n\n" );
	}

	fclose( fin );
	fclose( flog );

	return 0;
}


