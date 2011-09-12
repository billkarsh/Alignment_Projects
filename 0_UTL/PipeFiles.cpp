

#include	"PipeFiles.h"

#include	"Disk.h"
#include	"File.h"






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
/* ReadThmParams ------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	GETPRM_THM( addr, pat )									\
	if( !GetPrm( addr, pat, "ReadThmParams", f, flog ) )		\
		goto exit

// Read parameter file governing thumbnail matching.
//
// The layer index specifies an override file for a given layer.
// For example, if layer=5, the search order (in temp directory)
// is first '../thmparams_5.txt' then '../thmparams.txt' with
// no index.
//
bool ReadThmParams( ThmParams &T, int layer, FILE *flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../thmparams_%d.txt", layer );
	f = fopen( name, "r" );

	if( !f ) {
		sprintf( name, "../thmparams.txt" );
		f = fopen( name, "r" );
	}

	if( f ) {

		fprintf( flog, "\n---- Thumb parameters ----\n" );

		GETPRM_THM( &T.pxp.PXBRO, "PXBRO=%d" );
		GETPRM_THM( &T.pxp.PXDOG, "PXDOG=%c" );
		GETPRM_THM( &T.pxp.PXDOG_R1, "PXDOG_R1=%d" );
		GETPRM_THM( &T.pxp.PXDOG_R2, "PXDOG_R2=%d" );

		// finish Y/N booleans
		T.pxp.PXDOG = (T.pxp.PXDOG == 'Y');

		GETPRM_THM( &T.FLD, "FLD=%c" );
		GETPRM_THM( &T.SCALE, "SCALE=%lf" );
		GETPRM_THM( &T.XSCALE, "XSCALE=%lf" );
		GETPRM_THM( &T.YSCALE, "YSCALE=%lf" );
		GETPRM_THM( &T.SKEW, "SKEW=%lf" );
		GETPRM_THM( &T.SLOPPY_SL, "SLOPPY_SL=%c" );
		GETPRM_THM( &T.OLAP1D, "OLAP1D=%d" );
		GETPRM_THM( &T.OLAP2D_SL, "OLAP2D_SL=%d" );
		GETPRM_THM( &T.OLAP2D_XL, "OLAP2D_XL=%d" );
		GETPRM_THM( &T.HFANGDN_SL, "HFANGDN_SL=%lf" );
		GETPRM_THM( &T.HFANGDN_XL, "HFANGDN_XL=%lf" );
		GETPRM_THM( &T.HFANGPR_SL, "HFANGPR_SL=%lf" );
		GETPRM_THM( &T.HFANGPR_XL, "HFANGPR_XL=%lf" );
		GETPRM_THM( &T.QTRSH_SL, "QTRSH_SL=%lf" );
		GETPRM_THM( &T.QTRSH_XL, "QTRSH_XL=%lf" );
		GETPRM_THM( &T.RTRSH_SL, "RTRSH_SL=%lf" );
		GETPRM_THM( &T.RTRSH_XL, "RTRSH_XL=%lf" );
		GETPRM_THM( &T.TWEAKS, "TWEAKS=%c" );
		GETPRM_THM( &T.INPALN, "INPALN=%c" );
		GETPRM_THM( &T.DINPUT, "DINPUT=%d" );

		// finish Y/N booleans
		T.FLD		= (T.FLD == 'Y');
		T.SLOPPY_SL	= (T.SLOPPY_SL == 'Y');
		T.TWEAKS	= (T.TWEAKS == 'Y');
		T.INPALN	= (T.INPALN == 'Y');

		fprintf( flog, "\n" );

		ok = true;
	}
	else
		fprintf( flog, "ReadThmParams: Can't open [%s].\n", name );

exit:
	if( f )
		fclose( f );

	return ok;
}

/* --------------------------------------------------------------- */
/* ReadMeshParams ------------------------------------------------ */
/* --------------------------------------------------------------- */

#define	GETPRM_MSH( addr, pat )									\
	if( !GetPrm( addr, pat, "ReadMeshParams", f, flog ) )		\
		goto exit

// Read parameter file governing mesh deformation.
//
// The layer index specifies an override file for a given layer.
// For example, if layer=5, the search order (in temp directory)
// is first '../dmeshparams_5.txt' then '../dmeshparams.txt' with
// no index.
//
bool ReadMeshParams( MeshParams &M, int layer, FILE *flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../dmeshparams_%d.txt", layer );
	f = fopen( name, "r" );

	if( !f ) {
		sprintf( name, "../dmeshparams.txt" );
		f = fopen( name, "r" );
	}

	if( f ) {

		fprintf( flog, "\n---- dmesh parameters ----\n" );

		GETPRM_THM( &M.pxp.PXBRO, "PXBRO=%d" );
		GETPRM_THM( &M.pxp.PXDOG, "PXDOG=%c" );
		GETPRM_THM( &M.pxp.PXDOG_R1, "PXDOG_R1=%d" );
		GETPRM_THM( &M.pxp.PXDOG_R2, "PXDOG_R2=%d" );

		// finish Y/N booleans
		M.pxp.PXDOG = (M.pxp.PXDOG == 'Y');

		GETPRM_MSH( &M.FLD, "FLD=%c" );
		GETPRM_MSH( &M.DSL, "DSL=%c" );
		GETPRM_MSH( &M.DIT, "DIT=%lf" );
		GETPRM_MSH( &M.DAF, "DAF=%lf" );
		GETPRM_MSH( &M.DFT, "DFT=%lf" );
		GETPRM_MSH( &M.TMC, "TMC=%lf" );
		GETPRM_MSH( &M.TSC, "TSC=%lf" );
		GETPRM_MSH( &M.MNL, "MNL=%d" );
		GETPRM_MSH( &M.MTA, "MTA=%d" );
		GETPRM_MSH( &M.MMA, "MMA=%d" );
		GETPRM_MSH( &M.ONE, "ONE=%c" );
		GETPRM_MSH( &M.IFM, "IFM=%lf" );
		GETPRM_MSH( &M.FFM, "FFM=%lf" );
		GETPRM_MSH( &M.FYL, "FYL=%lf" );
		GETPRM_MSH( &M.CPD, "CPD=%lf" );
		GETPRM_MSH( &M.EMM, "EMM=%c" );
		GETPRM_MSH( &M.EMT, "EMT=%lf" );
		GETPRM_MSH( &M.WDI, "WDI=%c" );
		GETPRM_MSH( &M.LDA, "LDA=%lf" );
		GETPRM_MSH( &M.LDR, "LDR=%lf" );
		GETPRM_MSH( &M.LDC, "LDC=%lf" );
		GETPRM_MSH( &M.DXY, "DXY=%lf" );
		GETPRM_MSH( &M.WMT, "WMT=%c" );
		GETPRM_MSH( &M.WTT, "WTT=%c" );

		// finish Y/N booleans
		M.FLD = (M.FLD == 'Y');
		M.DSL = (M.DSL == 'Y');
		M.ONE = (M.ONE == 'Y');
		M.EMM = (M.EMM == 'Y');
		M.WDI = (M.WDI == 'Y');
		M.WMT = (M.WMT == 'Y');
		M.WTT = (M.WTT == 'Y');

		fprintf( flog, "\n" );

		ok = true;
	}
	else
		fprintf( flog, "ReadMeshParams: Can't open [%s].\n", name );

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
	Til2Img		&t2i,
	const char	*idb,
	int			layer,
	int			tile,
	FILE		*flog )
{
	char	name[2048];
	FILE	*f;
	int		ok = false;

	if( idb && idb[0] )
		sprintf( name, "%s/%d/TileToImage.txt", idb, layer );
	else
		sprintf( name, "../%d/TileToImage.txt", layer );

	f = fopen( name, "r" );

	if( f ) {

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

			sscanf( LS.line, "%d\t%lf\t%lf\t%lf"
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
	Til2FM		&t2f,
	const char	*idb,
	int			layer,
	int			tile,
	FILE		*flog )
{
	char	name[2048];
	FILE	*f;
	int		ok = false;

	if( idb && idb[0] )
		sprintf( name, "%s/%d/TileToFM.txt", idb, layer );
	else
		sprintf( name, "../%d/TileToFM.txt", layer );

	f = fopen( name, "r" );

	if( f ) {

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
	Til2FM		&t2f,
	const char	*idb,
	int			layer,
	int			tile )
{
	char	name[2048];
	FILE	*f;
	int		ok = false;

	if( idb && idb[0] )
		sprintf( name, "%s/%d/TileToFMD.txt", idb, layer );
	else
		sprintf( name, "../%d/TileToFMD.txt", layer );

	f = fopen( name, "r" );

	if( f ) {

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

				sscanf( LS.line, "%d\t%d\t%d\t%d\t%d"
				"\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&tpr.atl, &tpr.btl, &tpr.acr, &tpr.bcr, &tpr.err,
				&tpr.A, &tpr.Q, &tpr.R,
				&tpr.T.t[0], &tpr.T.t[1], &tpr.T.t[2],
				&tpr.T.t[3], &tpr.T.t[4], &tpr.T.t[5] );

				if( tpr.atl != atl || tpr.btl != btl ||
					tpr.acr != acr || tpr.bcr != bcr ) {

					continue;
				}

				fprintf( flog, "ReadThmPair: Got entry: "
				"A=%f, Q=%f, R=%f, T=[%f %f %f %f %f %f].\n",
				tpr.A, tpr.Q, tpr.R,
				tpr.T.t[0], tpr.T.t[1], tpr.T.t[2],
				tpr.T.t[3], tpr.T.t[4], tpr.T.t[5] );

				ok = true;
				goto exit;
			}

			fprintf( flog,
			"ReadThmPair: No entry for %d_%d@%d_%d; cr:%d-%d.\n",
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

				sscanf( LS.line, "%d\t%d\t%d\t%d\t%d"
				"\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&P.atl, &P.btl, &P.acr, &P.bcr, &P.err,
				&P.A, &P.Q, &P.R,
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
			fprintf( f, "%d\t%d\t%d\t%d\t%d"
				"\t%f\t%f\t%f"
				"\t%f\t%f\t%f\t%f\t%f\t%f\n",
				atl, btl, acr, bcr, tpr.err,
				tpr.A, tpr.Q, tpr.R,
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


