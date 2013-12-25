

#include	"lsq_Globals.h"
#include	"lsq_XArray.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"Timer.h"

#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray*		ME;
static const char	*gpath;
static vector<int>	giz;
static int			nthr;






/* --------------------------------------------------------------- */
/* _AFromIDB ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _AFromIDB( void* ithr )
{
	int	nz = vR.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		vector<Til2Img>	t2i;

		// Get rgn #1 tforms

		if( !IDBT2IGet_JustIDandT( t2i, idb, R.z ) )
			exit( 42 );

		x.resize( R.nr * 6 );

		int						nt = t2i.size();
		map<int,int>::iterator	en = R.m.end();

		// For each transform in IDB...

		for( int it = 0; it < nt; ++it ) {

			// Get its block start and limit {j0,jlim}

			const Til2Img&			T = t2i[it];
			map<int,int>::iterator	mi = R.m.find( T.id );
			int						j0, jlim;

			if( mi == en )
				continue;

			j0		= mi->second;
			jlim	= (++mi != en ? mi->second : R.nr);

			// Propagate rgn #1 tform to all block members

			for( int j = j0; j < jlim; ++j ) {

				if( R.pts[j].size() >= 3 ) {

					T.T.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* _AFromTxt ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// There are two differences in reading 'TXT' vs 'MET' formats:
//
// (1) The folder path name; but that's absorbed in gpath.
//
// (2) We scan each line for {id,r,T} but there is more beyond
//		in the 'MET' case. That's handled in CIDRA::FromFile()
//		by format '%*[^\n]\n' which reads and tosses everything
//		up to and including the terminator.
//

class CIDRA {
public:
	TAffine	A;
	int		id, r;
public:
	inline bool FromFile( FILE *f )
	{
		return 8 == fscanf( f,
			" %d %d %lf %lf %lf %lf %lf %lf%*[^\n]\n",
			&id, &r,
			&A.t[0], &A.t[1], &A.t[2],
			&A.t[3], &A.t[4], &A.t[5] );
	};
};


static bool Read_vA( vector<CIDRA> &vA, int z )
{
	char	buf[2048];
	FILE	*f;
	CIDRA	A;
	bool	nf = true;	// default = no folds

	sprintf( buf, "%s/X_A_%d.txt", gpath, z );
	f = FileOpenOrDie( buf, "r" );

	while( A.FromFile( f ) ) {

		if( A.r > 1 )
			nf = false;

		vA.push_back( A );
	}

	fclose( f );

	return nf;
}


static void* _AFromTxt( void* ithr )
{
	int	nz = vR.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		vector<CIDRA>	vA;
		int				nf = Read_vA( vA, R.z );

		x.resize( R.nr * 6 );

		int						na = vA.size();
		map<int,int>::iterator	en = R.m.end();

		if( nf ) {	// Propagate rgn #1 to all block members

			// For each transform in vA...

			for( int ia = 0; ia < na; ++ia ) {

				// Get its block start and limit {j0,jlim}

				const CIDRA&			A = vA[ia];
				map<int,int>::iterator	mi = R.m.find( A.id );
				int						j0, jlim;

				if( mi == en )
					continue;

				j0		= mi->second;
				jlim	= (++mi != en ? mi->second : R.nr);

				// Propagate rgn #1 tform to all block members

				for( int j = j0; j < jlim; ++j ) {

					if( R.pts[j].size() >= 3 ) {

						A.A.CopyOut( &x[j * 6] );
						R.used[j] = true;
					}
				}
			}
		}
		else {	// Move each A to its specified position

			// For each transform in vA...

			for( int ia = 0; ia < na; ++ia ) {

				// Get its location j

				const CIDRA&			A = vA[ia];
				map<int,int>::iterator	mi = R.m.find( A.id );
				int						j;

				if( mi == en )
					continue;

				j = mi->second + A.r - 1;

				if( R.pts[j].size() >= 3 ) {

					A.A.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* _HFromTxt ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// There are two differences in reading 'TXT' vs 'MET' formats:
//
// (1) The folder path name; but that's absorbed in gpath.
//
// (2) We scan each line for {id,r,T} but there is more beyond
//		in the 'MET' case. That's handled in CIDRH::FromFile()
//		by format '%*[^\n]\n' which reads and tosses everything
//		up to and including the terminator.
//

class CIDRH {
public:
	THmgphy	H;
	int		id, r;
public:
	inline bool FromFile( FILE *f )
	{
		return 10 == fscanf( f,
			" %d %d %lf %lf %lf %lf %lf %lf %lf %lf%*[^\n]\n",
			&id, &r,
			&H.t[0], &H.t[1], &H.t[2],
			&H.t[3], &H.t[4], &H.t[5],
			&H.t[6], &H.t[7] );
	};
};


static bool Read_vH( vector<CIDRH> &vH, int z )
{
	char	buf[2048];
	FILE	*f;
	CIDRH	H;
	bool	nf = true;	// default = no folds

	sprintf( buf, "%s/X_H_%d.txt", gpath, z );
	f = FileOpenOrDie( buf, "r" );

	while( H.FromFile( f ) ) {

		if( H.r > 1 )
			nf = false;

		vH.push_back( H );
	}

	fclose( f );

	return nf;
}


static void* _HFromTxt( void* ithr )
{
	int	nz = vR.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		vector<CIDRH>	vH;
		int				nf = Read_vH( vH, R.z );

		x.resize( R.nr * 6 );

		int						nh = vH.size();
		map<int,int>::iterator	en = R.m.end();

		if( nf ) {	// Propagate rgn #1 to all block members

			// For each transform in vH...

			for( int ih = 0; ih < nh; ++ih ) {

				// Get its block start and limit {j0,jlim}

				const CIDRH&			H = vH[ih];
				map<int,int>::iterator	mi = R.m.find( H.id );
				int						j0, jlim;

				if( mi == en )
					continue;

				j0		= mi->second;
				jlim	= (++mi != en ? mi->second : R.nr);

				// Propagate rgn #1 tform to all block members

				for( int j = j0; j < jlim; ++j ) {

					if( R.pts[j].size() >= 4 ) {

						H.H.CopyOut( &x[j * 6] );
						R.used[j] = true;
					}
				}
			}
		}
		else {	// Move each H to its specified position

			// For each transform in vH...

			for( int ih = 0; ih < nh; ++ih ) {

				// Get its location j

				const CIDRH&			H = vH[ih];
				map<int,int>::iterator	mi = R.m.find( H.id );
				int						j;

				if( mi == en )
					continue;

				j = mi->second + H.r - 1;

				if( R.pts[j].size() >= 4 ) {

					H.H.CopyOut( &x[j * 6] );
					R.used[j] = true;
				}
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* _XFromBin ----------------------------------------------------- */
/* --------------------------------------------------------------- */

static void ReadXBin( vector<double> &x, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/X_%c_%d.bin",
		gpath, (ME->NE == 6 ? 'A' : 'H'), z );

	f = FileOpenOrDie( buf, "rb" );
	fread( &x[0], sizeof(double), x.size(), f );
	fclose( f );
}


static void ReadUBin( vector<char> &u, int z )
{
	char	buf[2048];
	FILE	*f;

	sprintf( buf, "%s/U_%d.bin", gpath, z );
	f = FileOpenOrDie( buf, "rb" );
	fread( &u[0], sizeof(char), u.size(), f );
	fclose( f );
}


static void* _XFromBin( void* ithr )
{
	int	nz = vR.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[iz];
		vector<double>&	x = ME->X[iz];
		int				minpts = ME->NE / 2;

		x.resize( R.nr * ME->NE );
		ReadXBin( x, R.z );
		ReadUBin( R.used, R.z );

		for( int j = 0; j < R.nr; ++j )
			R.used[j] &= R.pts[j].size() >= minpts;
	}
}

/* --------------------------------------------------------------- */
/* _Updt --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void* _Updt( void* ithr )
{
	int	nz = giz.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[giz[iz]];
		vector<double>&	x = ME->X[giz[iz]];
		int				minpts = ME->NE / 2;

		ReadXBin( x, R.z );

		vector<char>	u( R.nr );
		ReadUBin( u, R.z );

		for( int j = 0; j < R.nr; ++j )
			R.used[j] &= u[j] && R.pts[j].size() >= minpts;
	}
}

/* --------------------------------------------------------------- */
/* _Save --------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void SaveXBin( vector<double> &x, int z )
{
	char	buf[64];
	FILE	*f;

	sprintf( buf, "%s/X_%c_%d.bin",
		gpath, (ME->NE == 6 ? 'A' : 'H'), z );

	f = FileOpenOrDie( buf, "wb" );
	fwrite( &x[0], sizeof(double), x.size(), f );
	fclose( f );
}


static void SaveUBin( vector<char> &u, int z )
{
	char	buf[64];
	FILE	*f;

	sprintf( buf, "%s/U_%d.bin", gpath, z );
	f = FileOpenOrDie( buf, "wb" );
	fwrite( &u[0], sizeof(char), u.size(), f );
	fclose( f );
}


static void* _Save( void* ithr )
{
	int	nz = giz.size();

	for( int iz = (long)ithr; iz < nz; iz += nthr ) {

		Rgns&			R = vR[giz[iz]];
		vector<double>&	x = ME->X[giz[iz]];

		SaveXBin( x, R.z );
		SaveUBin( R.used, R.z );
	}
}

/* --------------------------------------------------------------- */
/* Size ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Size( int ne )
{
	NE = ne;

	int	nz = vR.size();

	X.resize( nz );

	for( int iz = 0; iz < nz; ++iz )
		X[iz].resize( vR[iz].nr * ne );
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Load( const char *path )
{
	clock_t	t0 = StartTiming();

	int	nz = vR.size();

	X.resize( nz );

	ME		= this;
	gpath	= path;
	nthr	= (zolo != zohi ? 16 : 1);

	if( nthr > nz )
		nthr = nz;

	EZThreadproc	proc;
	const char		*sproc;

	if( !path ) {
		NE		= 6;
		proc	= _AFromIDB;
		sproc	= "AFromIDB";
	}
	else {

		const char *name = FileNamePtr( path );

		if( strstr( name, "X_A" ) ) {

			NE = 6;

			if( strstr( name, "X_A_TXT" ) ||
				strstr( name, "X_A_MET" ) ) {

				proc	= _AFromTxt;
				sproc	= "AFromTxt";
			}
			else if( strstr( name, "X_A_BIN" ) ) {
				proc	= _XFromBin;
				sproc	= "XFromBin";
			}
			else
				goto error;
		}
		else if( strstr( name, "X_H" ) ) {

			NE = 8;

			if( strstr( name, "X_H_TXT" ) ||
				strstr( name, "X_H_MET" ) ) {

				proc	= _HFromTxt;
				sproc	= "HFromTxt";
			}
			else if( strstr( name, "X_H_BIN" ) ) {
				proc	= _XFromBin;
				sproc	= "XFromBin";
			}
			else
				goto error;
		}
		else {
error:
			printf( "Unknown prior type [%s].\n", name );
			exit( 42 );
		}
	}

	if( !EZThreads( proc, nthr, 1, sproc ) )
		exit( 42 );

	StopTiming( stdout, sproc, t0 );
}

/* --------------------------------------------------------------- */
/* Updt ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Updt()
{
	giz.clear();

	for( int i = zolo; i < zilo; ++i )
		giz.push_back( i );

	for( int i = zihi + 1; i <= zohi; ++i )
		giz.push_back( i );

	int	nz = giz.size();

	if( !nz )
		return;

	clock_t	t0 = StartTiming();

	char	buf[32];
	sprintf( buf, "X_%c_BIN", (ME->NE == 6 ? 'A' : 'H') );

	ME		= this;
	gpath	= buf;
	nthr	= 16;

	if( nthr > nz )
		nthr = nz;

	if( !EZThreads( _Updt, nthr, 1, "Updt" ) )
		exit( 42 );

	StopTiming( stdout, "Updt", t0 );
}

/* --------------------------------------------------------------- */
/* Save ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void XArray::Save()
{
	giz.clear();

	for( int i = zilo; i <= zihi; ++i )
		giz.push_back( i );

	int	nz = giz.size();

	clock_t	t0 = StartTiming();

	char	buf[32];
	sprintf( buf, "X_%c_BIN", (ME->NE == 6 ? 'A' : 'H') );
	DskCreateDir( buf, stdout );

	ME		= this;
	gpath	= buf;
	nthr	= 16;

	if( nthr > nz )
		nthr = nz;

	if( !EZThreads( _Save, nthr, 1, "Save" ) )
		exit( 42 );

	StopTiming( stdout, "Save", t0 );
}

/* --------------------------------------------------------------- */
/* PriorIsAffine ------------------------------------------------- */
/* --------------------------------------------------------------- */

bool XArray::PriorIsAffine( const char *path )
{
	return (NULL != strstr( FileNamePtr( path ), "X_A" ));
}


