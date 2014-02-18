

#include	"Disk.h"
#include	"Timer.h"

#include	<dirent.h>
#include	<stdlib.h>

#include	<set>
using namespace std;


const int	MAXZRANGE = 10;


class CSubdirCat {
public:
	int			sx, sy,
				dx, dy;
	set<int>	zdown;
public:
	CSubdirCat() : sx(-1), sy(-1), dx(-1), dy(-1) {};
};

const char	*_d;
CSubdirCat	*_C;
int			_z;


static void FreeNamelist( struct dirent** &namelist, int n )
{
	if( namelist ) {

		while( n-- > 0 )
			free( namelist[n] );

		free( namelist );
		namelist = NULL;
	}
}


static void ScanThmPairs( const char *subdir )
{
	for( int i = _z  - MAXZRANGE; i < _z; ++i ) {

		char	path[2048];

		sprintf( path, "%s/%s/ThmPair_%d^%d.txt",
			_d, subdir, _z, i );

		if( DskExists( path ) )
			_C->zdown.insert( i );
	}
}


static int SorD( const struct dirent* E )
{
	if( E->d_name[0] == 'S' ) {

		int	x, y;
		if( 2 == sscanf( E->d_name + 1, "%d_%d", &x, &y ) ) {

			if( x > _C->sx )
				_C->sx = x;

			if( y > _C->sy )
				_C->sy = y;
		}
	}

	if( E->d_name[0] == 'D' ) {

		int	x, y;
		if( 2 == sscanf( E->d_name + 1, "%d_%d", &x, &y ) ) {

			if( x > _C->dx )
				_C->dx = x;

			if( y > _C->dy )
				_C->dy = y;

			ScanThmPairs( E->d_name );
		}
	}

	return 0;
}


static bool ScanThisZ( CSubdirCat &C, const char *top, int z )
{
	char	dir[2048];
	sprintf( dir, "%s/%d", top, z );
	_d = dir;
	_C = &C;
	_z = z;

	struct dirent **namelist = NULL;

	int	n = scandir( dir, &namelist, SorD, alphasort );

	FreeNamelist( namelist, n );

	return (n >= 0);
}


void Test()
{
	clock_t	t0 = StartTiming();

	for( int iz = 879; iz < 950; ++iz ) {

		CSubdirCat	C;

		ScanThisZ( C, "../..", iz );

		printf( "S %d %d D %d %d Z %d :",
			C.sx, C.sy, C.dx, C.dy, C.zdown.size() );

		set<int>::iterator	it;

		for( it = C.zdown.begin(); it != C.zdown.end(); ++it )
			printf( " %d", *it );

		printf( "\n" );
	}

	t0 = StopTiming( stdout, "Cat", t0 );
}

