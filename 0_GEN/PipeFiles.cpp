

#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"






/* --------------------------------------------------------------- */
/* OpenPairLog --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Rename stdout using image labels
//
void OpenPairLog( int alr, int atl, int blr, int btl )
{
	char	slog[256];

	sprintf( slog, "pair_%d_%d_@_%d_%d.log",
		alr, atl, blr, btl );

	freopen( slog, "a", stdout );
}

/* --------------------------------------------------------------- */
/* GetPrm -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool GetPrm(
	void		*v,
	const char	*fmt,
	const char	*caller,
	FILE		*fin,
	FILE		*flog )
{
	CLineScan	LS;
	bool		ok = false;

	while( LS.Get( fin ) > 0 ) {

		if( isspace( LS.line[0] ) || LS.line[0] == '#' )
			continue;
		else if( 1 == sscanf( LS.line, fmt, v ) ) {
			fprintf( flog, LS.line );
			ok = true;
			goto exit;
		}
		else {
			fprintf( flog, "%s: Expected [%s], found [%.23s...].\n",
			caller, fmt, LS.line );
			goto exit;
		}
	}

	fprintf( flog, "%s: Can't find [%s].\n", caller, fmt );

exit:
	return ok;
}

/* --------------------------------------------------------------- */
/* ReadMatchParams ----------------------------------------------- */
/* --------------------------------------------------------------- */

#define	GETPRM_MCH( addr, pat )									\
	if( !GetPrm( addr, pat, "ReadMatchParams", f, flog ) )		\
		goto exit

// Read parameter file governing thumbnail and mesh matching.
//
// The layer params specify an override file. E.g., if alr=5 and
// blr=4 the search order (in alignment 'temp' directory) is first
// 'matchparams_5_4.txt' then 'matchparams.txt' (no indices).
//
bool ReadMatchParams(
	MatchParams		&M,
	int				alr,
	int				blr,
	FILE			*flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../../matchparams_%d_%d.txt", alr, blr );
	f = fopen( name, "r" );

	if( !f ) {
		sprintf( name, "../../matchparams.txt" );
		f = fopen( name, "r" );
	}

	if( f ) {

		fprintf( flog, "\n---- Match parameters ----\n" );

		GETPRM_MCH( &M.PXBRO, "PXBRO=%d" );
		GETPRM_MCH( &M.PXLENS, "PXLENS=%c" );
		GETPRM_MCH( &M.PXDOG, "PXDOG=%c" );
		GETPRM_MCH( &M.PXDOG_R1, "PXDOG_R1=%d" );
		GETPRM_MCH( &M.PXDOG_R2, "PXDOG_R2=%d" );
		GETPRM_MCH( &M.FLD, "FLD=%c" );
		GETPRM_MCH( &M.PRETWEAK, "PRETWEAK=%c" );
		GETPRM_MCH( &M.SCALE, "SCALE=%lf" );
		GETPRM_MCH( &M.XSCALE, "XSCALE=%lf" );
		GETPRM_MCH( &M.YSCALE, "YSCALE=%lf" );
		GETPRM_MCH( &M.SKEW, "SKEW=%lf" );
		GETPRM_MCH( &M.MODE_SL, "MODE_SL=%c" );
		GETPRM_MCH( &M.MODE_XL, "MODE_XL=%c" );
		GETPRM_MCH( &M.TAB2DFM_SL, "TAB2DFM_SL=%c" );
		GETPRM_MCH( &M.TAB2DFM_XL, "TAB2DFM_XL=%c" );
		GETPRM_MCH( &M.XYCONF_SL, "XYCONF_SL=%lf" );
		GETPRM_MCH( &M.XYCONF_XL, "XYCONF_XL=%lf" );
		GETPRM_MCH( &M.THMDEC_SL, "THMDEC_SL=%d" );
		GETPRM_MCH( &M.THMDEC_XL, "THMDEC_XL=%d" );
		GETPRM_MCH( &M.OLAP1D_SL, "OLAP1D_SL=%d" );
		GETPRM_MCH( &M.OLAP1D_XL, "OLAP1D_XL=%d" );
		GETPRM_MCH( &M.OLAP2D_SL, "OLAP2D_SL=%d" );
		GETPRM_MCH( &M.OLAP2D_XL, "OLAP2D_XL=%d" );
		GETPRM_MCH( &M.NBMXHT_SL, "NBMXHT_SL=%lf" );
		GETPRM_MCH( &M.NBMXHT_XL, "NBMXHT_XL=%lf" );
		GETPRM_MCH( &M.HFANGDN_SL, "HFANGDN_SL=%lf" );
		GETPRM_MCH( &M.HFANGDN_XL, "HFANGDN_XL=%lf" );
		GETPRM_MCH( &M.HFANGPR_SL, "HFANGPR_SL=%lf" );
		GETPRM_MCH( &M.HFANGPR_XL, "HFANGPR_XL=%lf" );
		GETPRM_MCH( &M.RTRSH_SL, "RTRSH_SL=%lf" );
		GETPRM_MCH( &M.RTRSH_XL, "RTRSH_XL=%lf" );
		GETPRM_MCH( &M.TWEAKS, "TWEAKS=%c" );
		GETPRM_MCH( &M.LIMXY_SL, "LIMXY_SL=%d" );
		GETPRM_MCH( &M.LIMXY_XL, "LIMXY_XL=%d" );
		GETPRM_MCH( &M.OPT_SL, "OPT_SL=%c" );
		GETPRM_MCH( &M.RIT_SL, "RIT_SL=%lf" );
		GETPRM_MCH( &M.RIT_XL, "RIT_XL=%lf" );
		GETPRM_MCH( &M.RFA_SL, "RFA_SL=%lf" );
		GETPRM_MCH( &M.RFA_XL, "RFA_XL=%lf" );
		GETPRM_MCH( &M.RFT_SL, "RFT_SL=%lf" );
		GETPRM_MCH( &M.RFT_XL, "RFT_XL=%lf" );
		GETPRM_MCH( &M.TMC, "TMC=%lf" );
		GETPRM_MCH( &M.TSC, "TSC=%lf" );
		GETPRM_MCH( &M.MNL, "MNL=%d" );
		GETPRM_MCH( &M.MTA, "MTA=%d" );
		GETPRM_MCH( &M.MMA, "MMA=%d" );
		GETPRM_MCH( &M.ONE, "ONE=%c" );
		GETPRM_MCH( &M.IFM, "IFM=%lf" );
		GETPRM_MCH( &M.FFM, "FFM=%lf" );
		GETPRM_MCH( &M.FYL, "FYL=%lf" );
		GETPRM_MCH( &M.CPD, "CPD=%lf" );
		GETPRM_MCH( &M.EMM, "EMM=%c" );
		GETPRM_MCH( &M.EMT, "EMT=%lf" );
		GETPRM_MCH( &M.WDI, "WDI=%c" );
		GETPRM_MCH( &M.LDA, "LDA=%lf" );
		GETPRM_MCH( &M.LDR, "LDR=%lf" );
		GETPRM_MCH( &M.LDC, "LDC=%lf" );
		GETPRM_MCH( &M.DXY, "DXY=%lf" );
		GETPRM_MCH( &M.WMT, "WMT=%c" );
		GETPRM_MCH( &M.WTT, "WTT=%c" );

		// ensure upper case
		M.FLD		= toupper( M.FLD );
		M.MODE_SL	= toupper( M.MODE_SL );
		M.MODE_XL	= toupper( M.MODE_XL );

		// finish Y/N booleans
		M.PXLENS		= (toupper( M.PXLENS ) == 'Y');
		M.PXDOG			= (toupper( M.PXDOG ) == 'Y');
		M.PRETWEAK		= (toupper( M.PRETWEAK ) == 'Y');
		M.TAB2DFM_SL	= (toupper( M.TAB2DFM_SL ) == 'Y');
		M.TAB2DFM_XL	= (toupper( M.TAB2DFM_XL ) == 'Y');
		M.TWEAKS		= (toupper( M.TWEAKS ) == 'Y');
		M.OPT_SL		= (toupper( M.OPT_SL ) == 'Y');
		M.ONE			= (toupper( M.ONE ) == 'Y');
		M.EMM			= (toupper( M.EMM ) == 'Y');
		M.WDI			= (toupper( M.WDI ) == 'Y');
		M.WMT			= (toupper( M.WMT ) == 'Y');
		M.WTT			= (toupper( M.WTT ) == 'Y');

		fprintf( flog, "\n" );

		ok = true;
	}
	else
		fprintf( flog, "ReadMatchParams: Can't open [%s].\n", name );

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* IDBReadImgParams ---------------------------------------------- */
/* --------------------------------------------------------------- */

void IDBReadImgParams( string &idbpath, FILE *flog )
{
	idbpath.clear();

	FILE	*f = fopen( "../../imageparams.txt", "r" );

	if( f ) {

		CLineScan	LS;

		while( LS.Get( f ) > 0 ) {

			char	buf[2048];

			if( 1 == sscanf( LS.line, "IDBPATH %[^\n]", buf ) ) {

				idbpath = buf;
				goto close;
			}
		}

		fprintf( flog,
		"IDB: WARNING: imageparams.txt missing IDBPATH tag.\n" );

close:
		fclose( f );
	}
	else {
		fprintf( flog,
		"IDB: WARNING: Can't open imageparams.txt.\n" );
	}
}

/* --------------------------------------------------------------- */
/* IDBAllTil2Img ------------------------------------------------- */
/* --------------------------------------------------------------- */

bool IDBAllTil2Img(
	vector<Til2Img>	&t2i,
	const string	&idb,
	int				layer,
	FILE			*flog )
{
	char	name[2048];
	FILE	*f;
	int		ok = false;

	t2i.clear();

	sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), layer );

	if( f = fopen( name, "r" ) ) {

		CLineScan	LS;

		if( LS.Get( f ) <= 0 ) {
			fprintf( flog, "IDBAllTil2Img: Empty file [%s].\n", name );
			goto exit;
		}

		while( LS.Get( f ) > 0 ) {

			Til2Img	E;
			char	buf[2048];

			sscanf( LS.line,
			"%d\t%lf\t%lf\t%lf"
			"\t%lf\t%lf\t%lf\t%[^\t\n]",
			&E.tile,
			&E.T.t[0], &E.T.t[1], &E.T.t[2],
			&E.T.t[3], &E.T.t[4], &E.T.t[5],
			buf );

			E.path = buf;

			t2i.push_back( E );
		}

		ok = true;
	}
	else
		fprintf( flog, "IDBAllTil2Img: Can't open [%s].\n", name );

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* IDBTil2Img ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given IDBPATH/TileToImage file for this tile's data.
//
bool IDBTil2Img(
	Til2Img			&t2i,
	const string	&idb,
	int				layer,
	int				tile,
	const char		*forcepath,
	FILE			*flog )
{
// override provided

	if( forcepath ) {
		t2i.tile	= tile;
		t2i.T		= TAffine( 1,0,0,0,1,0 );
		t2i.path	= forcepath;
		return true;
	}

// standard way using idb

	char	name[2048];
	FILE	*f;
	int		ok = false;

	if( idb.empty() )
		sprintf( name, "../%d/TileToImage.txt", layer );
	else
		sprintf( name, "%s/%d/TileToImage.txt", idb.c_str(), layer );

	if( f = fopen( name, "r" ) ) {

		CLineScan	LS;

		if( LS.Get( f ) <= 0 ) {
			fprintf( flog, "IDBTil2Img: Empty file [%s].\n", name );
			goto exit;
		}

		while( LS.Get( f ) > 0 ) {

			char	buf[2048];

			t2i.tile = -1;

			sscanf( LS.line, "%d", &t2i.tile );

			if( t2i.tile != tile )
				continue;

			sscanf( LS.line,
			"%d\t%lf\t%lf\t%lf"
			"\t%lf\t%lf\t%lf\t%[^\t\n]",
			&t2i.tile,
			&t2i.T.t[0], &t2i.T.t[1], &t2i.T.t[2],
			&t2i.T.t[3], &t2i.T.t[4], &t2i.T.t[5],
			buf );

			t2i.path	= buf;
			ok			= true;
			goto exit;
		}

		fprintf( flog, "IDBTil2Img: No entry for [%d %d].\n",
		layer, tile );
	}
	else
		fprintf( flog, "IDBTil2Img: Can't open [%s].\n", name );

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* IDBTil2FM ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given IDBPATH/TileToFM file for this tile's data.
//
bool IDBTil2FM(
	Til2FM			&t2f,
	const string	&idb,
	int				layer,
	int				tile,
	FILE			*flog )
{
	char	name[2048];

// old-style layer/tile dir hierarchy

	if( idb.empty() ) {

		int	len;

		// try name as tif
		len = sprintf( name, "../%d/%d/fm.tif", layer, tile );

		if( !DskExists( name ) ) {
			// assume png
			name[len-3] = 'p';
			name[len-2] = 'n';
			name[len-1] = 'g';
		}

		t2f.tile	= tile;
		t2f.path	= name;
		return true;
	}

// new way using idb

	FILE	*f;
	int		ok = false;

	sprintf( name, "%s/%d/TileToFM.txt", idb.c_str(), layer );

	if( f = fopen( name, "r" ) ) {

		CLineScan	LS;

		if( LS.Get( f ) <= 0 ) {
			fprintf( flog, "IDBTil2FM: Empty file [%s].\n", name );
			goto exit;
		}

		while( LS.Get( f ) > 0 ) {

			char	buf[2048];

			t2f.tile = -1;

			sscanf( LS.line, "%d", &t2f.tile );

			if( t2f.tile != tile )
				continue;

			sscanf( LS.line, "%d\t%[^\t\n]", &t2f.tile, buf );
			t2f.path	= buf;
			ok			= true;
			goto exit;
		}

		fprintf( flog, "IDBTil2FM: No entry for [%d %d].\n",
		layer, tile );
	}
	else
		fprintf( flog, "IDBTil2FM: Can't open [%s].\n", name );

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* IDBTil2FMD ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given IDBPATH/TileToFMd file for this tile's data.
//
// Unlike IDBTil2FM, this function does not print error messages
// upon failure; it just returns false. In that case, the caller
// should attempt to get the standard fm.
//
bool IDBTil2FMD(
	Til2FM			&t2f,
	const string	&idb,
	int				layer,
	int				tile )
{
	char	name[2048];

// old-style layer/tile dir hierarchy

	if( idb.empty() ) {

		int	len;

		// try name as tif
		len = sprintf( name, "../%d/%d/fmd.tif", layer, tile );

		if( !DskExists( name ) ) {
			// assume png
			name[len-3] = 'p';
			name[len-2] = 'n';
			name[len-1] = 'g';
		}

		t2f.tile	= tile;
		t2f.path	= name;
		return true;
	}

// new way using idb

	FILE	*f;
	int		ok = false;

	sprintf( name, "%s/%d/TileToFMD.txt", idb.c_str(), layer );

	if( f = fopen( name, "r" ) ) {

		CLineScan	LS;

		if( LS.Get( f ) <= 0 )
			goto exit;

		while( LS.Get( f ) > 0 ) {

			char	buf[2048];

			t2f.tile = -1;

			sscanf( LS.line, "%d", &t2f.tile );

			if( t2f.tile != tile )
				continue;

			sscanf( LS.line, "%d\t%[^\t\n]", &t2f.tile, buf );
			t2f.path	= buf;
			ok			= true;
			goto exit;
		}
	}

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* PrintTil2Img -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print Til2Img entry; cAB = {'A', 'B'}.
//
void PrintTil2Img( FILE *flog, int cAB, const Til2Img &t2i )
{
	fprintf( flog, "Til2Img entry: %c"
	" T=[%7.4f %7.4f %8.2f %7.4f %7.4f %8.2f]"
	" path=[%s].\n",
	cAB,
	t2i.T.t[0], t2i.T.t[1], t2i.T.t[2],
	t2i.T.t[3], t2i.T.t[4], t2i.T.t[5],
	t2i.path.c_str() );
}

/* --------------------------------------------------------------- */
/* PrintTil2FM --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Print Til2FM entry; cAB = {'A', 'B'}.
//
void PrintTil2FM( FILE *flog, int cAB, const Til2FM &t2f )
{
	fprintf( flog, "Til2FM entry: %c path=[%s].\n",
	cAB, t2f.path.c_str() );
}

/* --------------------------------------------------------------- */
/* ReadThmPair --------------------------------------------------- */
/* --------------------------------------------------------------- */

bool ReadThmPair(
	ThmPair	&tpr,
	int		alr,
	int		atl,
	int		acr,
	int		blr,
	int		btl,
	int		bcr,
	FILE	*flog )
{
	CMutex	M;
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "tpr_%d_%d", alr, blr );

	if( M.Get( name ) ) {

		sprintf( name, "ThmPair_%d_@_%d.txt", alr, blr );
		f = fopen( name, "r" );

		if( f ) {

			CLineScan	LS;

			if( LS.Get( f ) <= 0 ) {

				fprintf( flog,
				"ReadThmPair: Empty file [%s].\n", name );
				goto exit;
			}

			while( LS.Get( f ) > 0 ) {

				sscanf( LS.line,
				"%d\t%d\t%d\t%d\t%d"
				"\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&tpr.atl, &tpr.btl, &tpr.acr, &tpr.bcr, &tpr.err,
				&tpr.A, &tpr.R,
				&tpr.T.t[0], &tpr.T.t[1], &tpr.T.t[2],
				&tpr.T.t[3], &tpr.T.t[4], &tpr.T.t[5] );

				if( tpr.atl != atl || tpr.btl != btl ||
					tpr.acr != acr || tpr.bcr != bcr ) {

					continue;
				}

				fprintf( flog, "ReadThmPair: Got entry: "
				"A=%f, R=%f, T=[%f %f %f %f %f %f].\n",
				tpr.A, tpr.R,
				tpr.T.t[0], tpr.T.t[1], tpr.T.t[2],
				tpr.T.t[3], tpr.T.t[4], tpr.T.t[5] );

				ok = true;
				goto exit;
			}

			fprintf( flog,
			"ReadThmPair: No entry for %d_%d@%d_%d; cr:%d-%d\n",
			alr, atl, blr, btl, acr, bcr );
		}
		else
			fprintf( flog, "ReadThmPair: Can't open [%s].\n", name );

exit:
		if( f )
			fclose( f );
	}

	M.Release();

	return ok;
}

/* --------------------------------------------------------------- */
/* ReadAllThmPair ------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return true if successfully read file. Even so, entry count
// may be zero, so always check tpr.size().
//
bool ReadAllThmPair(
	vector<ThmPair>	&tpr,
	int				alr,
	int				blr,
	FILE			*flog )
{
	CMutex	M;
	char	name[256];
	FILE	*f;
	int		ok = false;

	tpr.clear();

	sprintf( name, "tpr_%d_%d", alr, blr );

	if( M.Get( name ) ) {

		sprintf( name, "ThmPair_%d_@_%d.txt", alr, blr );
		f = fopen( name, "r" );

		if( f ) {

			CLineScan	LS;

			if( LS.Get( f ) <= 0 ) {

				fprintf( flog,
				"ReadThmPair: Empty file [%s].\n", name );
				goto exit;
			}

			while( LS.Get( f ) > 0 ) {

				ThmPair		P;

				sscanf( LS.line,
				"%d\t%d\t%d\t%d\t%d"
				"\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&P.atl, &P.btl, &P.acr, &P.bcr, &P.err,
				&P.A, &P.R,
				&P.T.t[0], &P.T.t[1], &P.T.t[2],
				&P.T.t[3], &P.T.t[4], &P.T.t[5] );

				tpr.push_back( P );
			}

			ok = true;
		}
		else
			fprintf( flog, "ReadThmPair: Can't open [%s].\n", name );

exit:
		if( f )
			fclose( f );
	}

	M.Release();

	return ok;
}

/* --------------------------------------------------------------- */
/* WriteThmPairHdr ----------------------------------------------- */
/* --------------------------------------------------------------- */

void WriteThmPairHdr( FILE *f )
{
	fprintf( f,
	"Atl\tBtl\tAcr\tBcr\tErr\tDeg\tR"
	"\tT0\tT1\tX\tT3\tT4\tY\n" );
}

/* --------------------------------------------------------------- */
/* CreateJobsDir ------------------------------------------------- */
/* --------------------------------------------------------------- */

void CreateJobsDir(
	const char	*lyrdir,
	int			ix,
	int			iy,
	int			za,
	int			zb,
	FILE		*flog )
{
    char	name[2048];
    int		len;

// Create dir
	len = sprintf( name, "%s/%c%d_%d",
			lyrdir, (za == zb ? 'S' : 'D'), ix, iy );
	DskCreateDir( name, flog );

// Create ThmPair file
	if( zb >= 0 ) {
		sprintf( name + len, "/ThmPair_%d_@_%d.txt", za, zb );
		FILE	*f = FileOpenOrDie( name, "w", flog );
		WriteThmPairHdr( f );
		fclose( f );
	}
}

/* --------------------------------------------------------------- */
/* WriteThmPair -------------------------------------------------- */
/* --------------------------------------------------------------- */

void WriteThmPair(
	const ThmPair	&tpr,
	int				alr,
	int				atl,
	int				acr,
	int				blr,
	int				btl,
	int				bcr )
{
	CMutex	M;
	char	name[256];

	sprintf( name, "tpr_%d_%d", alr, blr );

	if( M.Get( name ) ) {

		sprintf( name, "ThmPair_%d_@_%d.txt", alr, blr );
		FILE *f = fopen( name, "a" );

		if( f ) {
			fprintf( f,
				"%d\t%d\t%d\t%d\t%d"
				"\t%f\t%f"
				"\t%f\t%f\t%f\t%f\t%f\t%f\n",
				atl, btl, acr, bcr, tpr.err,
				tpr.A, tpr.R,
				tpr.T.t[0], tpr.T.t[1], tpr.T.t[2],
				tpr.T.t[3], tpr.T.t[4], tpr.T.t[5] );
			fflush( f );
			fclose( f );
		}
	}

	M.Release();
}

/* --------------------------------------------------------------- */
/* ZIDFromFMPath ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Given path to a foldmask (or other) file located in the
// standard working directory hierarchy:
//
//	e.g. '/groups/.../temp/z/id/fm.tif'
//	e.g. '../z/id/fm.png'
//
// attempt to extract the z and tile-id values.
//
// Return true if success.
//
bool ZIDFromFMPath( int &z, int &id, const char *path )
{
	const char	*s;
	int			nprev;

// back over filename

	s = path + (nprev = strlen( path ));

	while( nprev-- > 0 && *--s != '/' )
		;

	if( *s != '/' )
		return false;

// back over id and translate

	while( nprev-- > 0 && *--s != '/' )
		;

	if( nprev <= 0 )
		return false;

	if( 1 != sscanf( s + 1, "%d", &id ) )
		return false;

// back over z and translate

	while( nprev-- > 0 && *--s != '/' )
		;

	if( *s != '/' )
		return false;

	if( 1 != sscanf( s + 1, "%d", &z ) )
		return false;

	return true;
}

/* --------------------------------------------------------------- */
/* LoadTAffineTbl_AllZ ------------------------------------------- */
/* --------------------------------------------------------------- */

// From an LSQ-style TAffineTable file, fill out a table of
// TAffine[mapped by (z,id)] and a table of all unique z-values.
//
void LoadTAffineTbl_AllZ(
	map<MZID,TAffine>	&Tmap,
	set<int>			&Zset,
	const char			*path,
	FILE				*flog )
{
	FILE		*f	= FileOpenOrDie( path, "r", flog );
	CLineScan	LS;

	for(;;) {

		if( LS.Get( f ) <= 0 )
			break;

		MZID	zid;
		TAffine	T;
		int		rgn;

		sscanf( LS.line, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&zid.z, &zid.id, &rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5] );

		if( rgn != 1 )
			continue;

		Zset.insert( zid.z );

		Tmap[zid] = T;
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* LoadTAffineTbl_ThisZ ------------------------------------------ */
/* --------------------------------------------------------------- */

void LoadTAffineTbl_ThisZ(
	map<MZIDR,TAffine>	&Tmap,
	int					z,
	const char			*path,
	FILE				*flog )
{
	FILE		*f	= FileOpenOrDie( path, "r", flog );
	CLineScan	LS;

	for(;;) {

		if( LS.Get( f ) <= 0 )
			break;

		MZIDR	zir;
		TAffine	T;

		zir.z = atoi( LS.line );

		if( zir.z < z )
			continue;

		if( zir.z > z )
			break;

		sscanf( LS.line, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		&zir.z, &zir.id, &zir.rgn,
		&T.t[0], &T.t[1], &T.t[2],
		&T.t[3], &T.t[4], &T.t[5] );

		Tmap[zir] = T;
	}

	fclose( f );
}


