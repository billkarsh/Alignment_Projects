

#include	"lsq_Globals.h"

#include	"File.h"
#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

string			idb;
vector<Layer>	vL;
map<int,int>	mZ;
vector<Rgns>	vR;
vector<CorrPnt>	vC;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void GetIDB( const char *tempdir )
{
	char	buf[2048];
	sprintf( buf, "%s/imageparams.txt", tempdir );

	FILE	*f = fopen( buf, "r" );

	if( f ) {

		CLineScan	LS;

		while( LS.Get( f ) > 0 ) {

			if( 1 == sscanf( LS.line, "IDBPATH %[^\n]", buf ) ) {

				idb = buf;
				goto close;
			}
		}

		printf( "Globals: imageparams.txt missing IDBPATH tag.\n" );

close:
		fclose( f );
	}
	else
		printf( "Globals: Can't open imageparams.txt.\n" );

	if( idb.empty() )
		exit( 42 );
}


void MapZPair( int &ia, int &ib, int za, int zb )
{
	static int	cached_za = -1, cached_zb = -1,
				cached_ia,      cached_ib;

	if( za != cached_za ) {
		cached_ia = mZ.find( za )->second;
		cached_za = za;
	}

	ia = cached_ia;

	if( za == zb )
		ib = ia;
	else {

		if( zb != cached_zb ) {
			cached_ib = mZ.find( zb )->second;
			cached_zb = zb;
		}

		ib = cached_ib;
	}
}


Rgns::Rgns( int z )
{
	this->z		= z;
	cached_id	= -1;
	nr			= IDBGetIDRgnMap( m, idb, z );
	used.resize( nr, 0 );
}


int Rgns::Map( int id, int r )
{
	if( id != cached_id ) {
		cached_i  = m.find( id )->second + r - 1;
		cached_id = id;
	}

	return cached_i;
}


