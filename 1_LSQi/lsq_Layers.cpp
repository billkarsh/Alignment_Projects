

#include	"../1_LSQi/lsq_Layers.h"
#include	"File.h"
#include	"Disk.h"
#include	"Timer.h"

#include	<dirent.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const char	*_d;
static Layer		*_L;






static void FreeNamelist( struct dirent** &namelist, int n )
{
	if( namelist ) {

		while( n-- > 0 )
			free( namelist[n] );

		free( namelist );
		namelist = NULL;
	}
}


static int ThmPair( const struct dirent* E )
{
	if( E->d_name[0] == 'T' && E->d_name[1] == 'h' ) {

		int	za, zb;

		if( 2 == sscanf( E->d_name, "ThmPair_%d_@_%d", &za, &zb )
			&& za == _L->z ) {

			_L->zdown.insert( zb );
		}
	}

	return 0;
}


static void ScanThmPairs( const char *subdir )
{
	char	dir[2048];
	sprintf( dir, "%s/%s", _d, subdir );

	struct dirent **namelist = NULL;

	int	n = scandir( dir, &namelist, ThmPair, alphasort );

	FreeNamelist( namelist, n );
}


static int SorD( const struct dirent* E )
{
	if( E->d_name[0] == 'S' ) {

		int	x, y;
		if( 2 == sscanf( E->d_name + 1, "%d_%d", &x, &y ) ) {

			if( x > _L->sx )
				_L->sx = x;

			if( y > _L->sy )
				_L->sy = y;
		}
	}

	if( E->d_name[0] == 'D' ) {

		int	x, y;
		if( 2 == sscanf( E->d_name + 1, "%d_%d", &x, &y ) ) {

			if( x > _L->dx )
				_L->dx = x;

			if( y > _L->dy )
				_L->dy = y;

			ScanThmPairs( E->d_name );
		}
	}

	return 0;
}


static bool ScanThisZ( Layer &L, const char *top, int z )
{
	char	dir[2048];
	sprintf( dir, "%s/%d", top, z );

	_d	= dir;
	_L	= &L;
	L.z	= z;

	struct dirent **namelist = NULL;

	int	n = scandir( dir, &namelist, SorD, alphasort );

	FreeNamelist( namelist, n );

	return (n >= 0);	// ok (no error)
}


static void NewCat(
	vector<Layer>	&vL,
	const char		*tempdir,
	const char		*cachedir,
	int				zolo,
	int				zohi )
{
	DskCreateDir( cachedir, stdout );

	char	buf[2048];
	sprintf( buf, "%s/catalog.txt", cachedir );
	FILE	*f = FileOpenOrDie( buf, "w" );

// Query range header
	fprintf( f, "zo %d %d\n", zolo, zohi );

	for( int z = zolo; z <= zohi; ++z ) {

		Layer	L;

		ScanThisZ( L, tempdir, z );

		if( L.sx > -1 ) {

			fprintf( f, "%d S %d %d D %d %d Z %ld :",
				z, L.sx, L.sy, L.dx, L.dy, L.zdown.size() );

			for( set<int>::iterator	it = L.zdown.begin();
				it != L.zdown.end();
				++it ) {

				fprintf( f, " %d", *it );
			}

			fprintf( f, "\n" );

			vL.push_back( L );
		}
	}

	fclose( f );
}


static bool LoadCat(
	vector<Layer>	&vL,
	const char		*cachedir,
	int				zolo,
	int				zohi )
{
	char	buf[2048];
	sprintf( buf, "%s/catalog.txt", cachedir );

	if( !DskExists( buf ) )
		return false;

	FILE		*f = FileOpenOrDie( buf, "r" );
	CLineScan	LS;

// Is file appropriate range?

	if( LS.Get( f ) > 0 ) {

		int	lo, hi;

		if( 2 != sscanf( LS.line, "zo %d %d", &lo, &hi ) )
			goto fail;

		if( zolo < lo || zohi > hi )
			goto exit;
	}
	else
		goto fail;

// Read entries

	while( LS.Get( f ) > 0 ) {

		Layer	L;
		int		K, k, nz;

		if( 1 != sscanf( LS.line, "%d%n", &L.z, &K ) )
			goto fail;

		if( L.z > zohi )
			break;

		if( L.z < zolo )
			continue;

		if( 5 != sscanf( LS.line+K, " S %d %d D %d %d Z %d :%n",
					&L.sx, &L.sy, &L.dx, &L.dy, &nz, &k ) ) {

			goto fail;
		}

		while( nz-- > 0 ) {

			int	z;

			K += k;

			if( 1 != sscanf( LS.line+K, "%d%n", &z, &k ) )
				goto fail;

			L.zdown.insert( z );
		}

		vL.push_back( L );
	}

exit:
	if( f )
		fclose( f );

	return !vL.empty();

fail:
	printf( "Catalog: Remaking due to format error.\n" );

	vL.clear();

	if( f )
		fclose( f );

	return false;
}


void LayerCat(
	vector<Layer>	&vL,
	const char		*tempdir,
	const char		*cachedir,
	int				zolo,
	int				zohi,
	bool			catclr )
{
	printf( "\n---- Cataloging ----\n" );

	clock_t	t0 = StartTiming();

	if( catclr || !LoadCat( vL, cachedir, zolo, zohi ) )
		NewCat( vL, tempdir, cachedir, zolo, zohi );

	if( vL.empty() ) {
		printf( "Catalog: No catalog data in range.\n" );
		exit( 42 );
	}

	StopTiming( stdout, "Catalog", t0 );
}


