

#include	"lsq_Catalog.h"
#include	"File.h"
#include	"Disk.h"

#include	<dirent.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static const char	*_d;
static CSubdirCat	*_C;






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
			&& za == _C->z ) {

			_C->zdown.insert( zb );
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

	_d	= dir;
	_C	= &C;
	C.z	= z;

	struct dirent **namelist = NULL;

	int	n = scandir( dir, &namelist, SorD, alphasort );

	FreeNamelist( namelist, n );

	return (n >= 0);	// ok (no error)
}


static void NewCat(
	vector<CSubdirCat>	&vC,
	const char			*tempdir,
	int					zolo,
	int					zohi )
{
	FILE	*f = FileOpenOrDie( "catalog.txt", "w" );

// Query range header
	fprintf( f, "zo %d %d\n", zolo, zohi );

	for( int z = zolo; z <= zohi; ++z ) {

		CSubdirCat	C;

		ScanThisZ( C, tempdir, z );

		if( C.sx > -1 ) {

			fprintf( f, "%d S %d %d D %d %d Z %d :",
				z, C.sx, C.sy, C.dx, C.dy, C.zdown.size() );

			for( set<int>::iterator	it = C.zdown.begin();
				it != C.zdown.end();
				++it ) {

				fprintf( f, " %d", *it );
			}

			fprintf( f, "\n" );

			vC.push_back( C );
		}
	}

	fclose( f );
}


static bool LoadCat(
	vector<CSubdirCat>	&vC,
	int					zolo,
	int					zohi )
{
	if( !DskExists( "catalog.txt" ) )
		return false;

	FILE		*f = FileOpenOrDie( "catalog.txt", "r" );
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

		CSubdirCat	C;
		int			K, k, nz;

		if( 1 != sscanf( LS.line, "%d%n", &C.z, &K ) )
			goto fail;

		if( C.z > zohi )
			break;

		if( C.z < zolo )
			continue;

		if( 5 != sscanf( LS.line+K, " S %d %d D %d %d Z %d :%n",
					&C.sx, &C.sy, &C.dx, &C.dy, &nz, &k ) ) {

			goto fail;
		}

		while( nz-- > 0 ) {

			int	z;

			K += k;

			if( 1 != sscanf( LS.line+K, "%d%n", &z, &k ) )
				goto fail;

			C.zdown.insert( z );
		}

		vC.push_back( C );
	}

exit:
	if( f )
		fclose( f );

	return !vC.empty();

fail:
	printf( "Catalog: Remaking due to format error.\n" );

	vC.clear();

	if( f )
		fclose( f );

	return false;
}


void CatPoints(
	vector<CSubdirCat>	&vC,
	const char			*tempdir,
	int					zolo,
	int					zohi,
	bool				clrcat )
{
	if( clrcat || !LoadCat( vC, zolo, zohi ) )
		NewCat( vC, tempdir, zolo, zohi );

	if( vC.empty() ) {
		printf( "Catalog: No catalog data in range.\n" );
		exit( 42 );
	}
}


