

#include	"lsq_ReadPts.h"

#include	"File.h"
#include	"PipeFiles.h"


/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

map<MZID,int>	nConRgn;	// # connected-rgn this image
int				gW = 4056,	// image coords
				gH = 4056;






/* --------------------------------------------------------------- */
/* FindOrAdd ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return zero-based index for given RGN.
//
// If already stored, that index is returned. Else, new
// entries are created with index (nr); nr is incremented.
//
static int FindOrAdd( map<MZIDR,int> &m, int &nr, const RGN &R )
{
// Already mapped?

	MZIDR						key( R.z, R.id, R.rgn );
	map<MZIDR,int>::iterator	it = m.find( key );

	if( it != m.end() )
		return it->second;

// No - add to map

	m[key] = nr;

// Add to vRgn

	vRgn.push_back( R );

	return nr++;
}

/* --------------------------------------------------------------- */
/* ReadPts_StrTags ----------------------------------------------- */
/* --------------------------------------------------------------- */

void ReadPts_StrTags(
	FILE		*FOUT,
	CNX			*cnx,
	SML			*sml,
	int			(*IDFromName)( const char *name ),
	const char	*dirfile,
	const char	*ptsfile )
{
	printf( "---- Read pts ----\n" );

	FILE		*f = FileOpenOrDie( ptsfile, "r" );
	CLineScan	LS;
	DIR			dir;	// map name strings to z layers

	dir.ReadDIRFile( dirfile, FOUT );

	map<MZIDR,int>	mapRGN;
	int				nr = 0, nlines = 0;

	for(;;) {

		char	name1[2048], name2[2048];

		if( LS.Get( f ) <= 0 )
			break;

		++nlines;

		if( !strncmp( LS.line, "CPOINT", 6 ) ) {

			char	key1[32], key2[32];
			Point	p1, p2;

			if( 8 != sscanf( LS.line + 7,
						"'%[^']' %s %lf %lf '%[^']' %s %lf %lf",
						name1, key1, &p1.x, &p1.y,
						name2, key2, &p2.x, &p2.y ) ) {

				printf(
				"WARNING: 'CPOINT' format error; line %d\n",
				nlines );

				continue;
			}

			RGN	R1( name1, key1 );
			RGN	R2( name2, key2 );
			int r1 = FindOrAdd( mapRGN, nr, R1 );
			int r2 = FindOrAdd( mapRGN, nr, R2 );

			cnx->AddCorrespondence( r1, r2 );
			sml->AddPOINTPair( r1, p1, r2, p2 );

			vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
		}
		else if( !strncmp( LS.line, "POINT", 5 ) ) {

			Point	p1, p2;

			if( 6 != sscanf( LS.line + 6,
						"%s %lf %lf %s %lf %lf",
						name1, &p1.x, &p1.y,
						name2, &p2.x, &p2.y ) ) {

				printf(
				"WARNING: 'POINT' format error; line %d\n",
				nlines );

				continue;
			}

			RGN	R1( name1, dir, IDFromName( name1 ) );
			RGN	R2( name2, dir, IDFromName( name2 ) );
			int r1 = FindOrAdd( mapRGN, nr, R1 );
			int r2 = FindOrAdd( mapRGN, nr, R2 );

			cnx->AddCorrespondence( r1, r2 );
			sml->AddPOINTPair( r1, p1, r2, p2 );

			vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
		}
		else if( !strncmp( LS.line, "FOLDMAP", 7 ) ) {

			int	z, id, nrgn = -1;

			sscanf( LS.line + 8, "'%*[^']' %s %d", name1, &nrgn );
			ZIDFromFMPath( z, id, name1 );
			nConRgn[MZID( z, id )] = nrgn;

			fprintf( FOUT, LS.line );
		}
		else if( !strncmp( LS.line, "IMAGESIZE", 9 ) ) {

			if( 2 != sscanf( LS.line + 10, "%d %d", &gW, &gH ) ) {
				printf( "Bad IMAGESIZE line '%s'.\n", LS.line );
				exit( 42 );
			}

			fprintf( FOUT, LS.line );
			printf( LS.line );
		}
		else {

			char	*s = strchr( LS.line, ' ' );

			if( s )
				*s = 0;

			printf(
			"WARNING: Unknown entry type; '%s' line %d\n",
			LS.line, nlines );
		}
	}

	fclose( f );

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* IsDaviCorner -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true to reject corner-corner matches.
//
static bool IsDaviCorner( const RGN &R1, const RGN &R2 )
{
	if( R1.z != R2.z )
		return false;

	const char	*c, *n;
	int			row1, col1, row2, col2;

	n = FileNamePtr( R1.GetName() );
	c = strstr( n, "col" );
	sscanf( c, "col%d_row%d", &col1, &row1 );

	n = FileNamePtr( R2.GetName() );
	c = strstr( n, "col" );
	sscanf( c, "col%d_row%d", &col2, &row2 );

	return ((row2 - row1) * (col2 - col1) != 0);
}

/* --------------------------------------------------------------- */
/* RejectPair ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Chance to optionally apply rejection criteria against
// point pairs.
//
// Return true to reject.
//
static bool RejectPair( const RGN &R1, const RGN &R2 )
{
// ------------------------------------
// reject layer 30 and 38
#if 0
	if( R1.z == 30 || R1.z == 38 || R2.z == 30 || R2.z == 38 )
		return true;
#endif
// ------------------------------------

// ------------------------------------
// accept only col/row subsection
#if 0
	const char	*c, *n;
	int			row, col;

	n = FileNamePtr( R1.GetName() );
	if( c = strstr( n, "col" ) ) {
		sscanf( c, "col%d_row%d", &col, &row );
//		if( col < 33 || col > 39 || row < 85 || row > 91 )
		if( col < 25-4 || col > 25+4 || row < 19-4 || row > 19+4 )
//		if( row > 2 || col > 2 )
			return true;

//		if( col == 44 && row == 19 )
//		if( col == 48 && row == 16 )
//			return true;
	}

	n = FileNamePtr( R2.GetName() );
	if( c = strstr( n, "col" ) ) {
		sscanf( c, "col%d_row%d", &col, &row );
//		if( col < 33 || col > 39 || row < 85 || row > 91 )
		if( col < 25-4 || col > 25+4 || row < 19-4 || row > 19+4 )
//		if( row > 2 || col > 2 )
			return true;

//		if( col == 44 && row == 19 )
//		if( col == 48 && row == 16 )
//			return true;
	}
#endif
// ------------------------------------

	return false;
}

/* --------------------------------------------------------------- */
/* ReadPts_NumTags ----------------------------------------------- */
/* --------------------------------------------------------------- */

void ReadPts_NumTags(
	FILE		*FOUT,
	CNX			*cnx,
	SML			*sml,
	const char	*ptsfile,
	int			davinocorn )
{
	printf( "---- Read pts ----\n" );

	FILE		*f = FileOpenOrDie( ptsfile, "r" );
	CLineScan	LS;

	map<MZIDR,int>	mapRGN;
	int				nr = 0, nlines = 0;

	for(;;) {

		char	name1[2048], name2[2048];

		if( LS.Get( f ) <= 0 )
			break;

		++nlines;

		if( !strncmp( LS.line, "CPOINT2", 7 ) ) {

			char	key1[32], key2[32];
			Point	p1, p2;

			if( 6 != sscanf( LS.line + 8,
						"%s %lf %lf %s %lf %lf",
						key1, &p1.x, &p1.y,
						key2, &p2.x, &p2.y ) ) {

				printf(
				"WARNING: 'CPOINT2' format error; line %d\n",
				nlines );

				continue;
			}

			RGN	R1( key1 );
			RGN	R2( key2 );

			if( davinocorn && IsDaviCorner( R1, R2 ) )
				continue;

			if( RejectPair( R1, R2 ) )
				continue;

			int r1 = FindOrAdd( mapRGN, nr, R1 );
			int r2 = FindOrAdd( mapRGN, nr, R2 );

			cnx->AddCorrespondence( r1, r2 );
			sml->AddPOINTPair( r1, p1, r2, p2 );

			vAllC.push_back( Constraint( r1, p1, r2, p2 ) );
		}
		else if( !strncmp( LS.line, "FOLDMAP2", 8 ) ) {

			int	z, id, nrgn;

			sscanf( LS.line + 9, "%d.%d %d", &z, &id, &nrgn );
			nConRgn[MZID( z, id )] = nrgn;
		}
		else if( !strncmp( LS.line, "IDBPATH", 7 ) ) {

			char	buf[2048];

			if( !sscanf( LS.line + 8, "%s", buf ) ) {

				printf( "Bad IDBPATH line '%s'.\n", LS.line );
				exit( 42 );
			}

			idb = buf;
			fprintf( FOUT, LS.line );
			printf( LS.line );
		}
		else if( !strncmp( LS.line, "IMAGESIZE", 9 ) ) {

			if( 2 != sscanf( LS.line + 10, "%d %d", &gW, &gH ) ) {
				printf( "Bad IMAGESIZE line '%s'.\n", LS.line );
				exit( 42 );
			}

			fprintf( FOUT, LS.line );
			printf( LS.line );
		}
		else {

			char	*s = strchr( LS.line, ' ' );

			if( s )
				*s = 0;

			printf(
			"WARNING: Unknown entry type; '%s' line %d\n",
			LS.line, nlines );
		}
	}

	fclose( f );

	printf( "\n" );

	IDBTil2ImgClear();
}


