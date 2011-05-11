

#include	"PipeFiles.h"

#include	"Disk.h"






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
	char	line[1024];
	bool	ok = false;

	while( fgets( line, sizeof(line), fin ) ) {

		if( 1 == sscanf( line, fmt, v ) ) {
			fprintf( flog, line );
			ok = true;
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
	if( !GetPrm(addr, pat, "ReadThmParams", f, flog) )			\
		goto exit

bool ReadThmParams( ThmParams &T, FILE *flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../thmparams.txt" );
	f = fopen( name, "r" );

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
		GETPRM_THM( &T.OLAP1D, "OLAP1D=%d" );
		GETPRM_THM( &T.OLAP2D_SL, "OLAP2D_SL=%d" );
		GETPRM_THM( &T.OLAP2D_XL, "OLAP2D_XL=%d" );
		GETPRM_THM( &T.HFANG_SL, "HFANG_SL=%lf" );
		GETPRM_THM( &T.HFANG_XL, "HFANG_XL=%lf" );
		GETPRM_THM( &T.QTRSH_SL, "QTRSH_SL=%lf" );
		GETPRM_THM( &T.QTRSH_XL, "QTRSH_XL=%lf" );
		GETPRM_THM( &T.RTRSH_SL, "RTRSH_SL=%lf" );
		GETPRM_THM( &T.RTRSH_XL, "RTRSH_XL=%lf" );
		GETPRM_THM( &T.TWEAKS, "TWEAKS=%c" );
		GETPRM_THM( &T.INPALN, "INPALN=%c" );
		GETPRM_THM( &T.DINPUT, "DINPUT=%d" );

		// finish Y/N booleans
		T.FLD		= (T.FLD == 'Y');
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
	if( !GetPrm(addr, pat, "ReadMeshParams", f, flog) )			\
		goto exit

bool ReadMeshParams( MeshParams &M, FILE *flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../dmeshparams.txt" );
	f = fopen( name, "r" );

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
/* ReadTil2Img --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Scan given TileToImage file for this tile's data.
//
bool ReadTil2Img(
	Til2Img	&t2i,
	int		layer,
	int		tile,
	FILE	*flog )
{
	char	name[256];
	FILE	*f;
	int		ok = false;

	sprintf( name, "../%d/TileToImage.txt", layer );
	f = fopen( name, "r" );

	if( f ) {

		char	line[2048];

		if( !fgets( line, sizeof(line), f ) ) {
			fprintf( flog, "ReadTil2Img: Empty file [%s].\n", name );
			goto exit;
		}

		while( fgets( line, sizeof(line), f ) ) {

			t2i.tile = -1;

			sscanf( line, "%d", &t2i.tile );

			if( t2i.tile != tile )
				continue;

			sscanf( line, "%d\t%lf\t%lf\t%lf"
			"\t%lf\t%lf\t%lf\t%s",
			&t2i.tile,
			&t2i.T.t[0], &t2i.T.t[1], &t2i.T.t[2],
			&t2i.T.t[3], &t2i.T.t[4], &t2i.T.t[5],
			t2i.path );

			ok = true;
			goto exit;
		}

		fprintf( flog, "ReadTil2Img: No entry for [%d %d].\n",
		layer, tile );
	}
	else
		fprintf( flog, "ReadTil2Img: Can't open [%s].\n", name );

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
	t2i.path );
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

			char	line[256];

			if( !fgets( line, sizeof(line), f ) ) {

				fprintf( flog,
				"ReadThmPair: Empty file [%s].\n", name );
				goto exit;
			}

			while( fgets( line, sizeof(line), f ) ) {

				sscanf( line, "%d\t%d\t%d\t%d\t%d"
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

			char	line[256];

			if( !fgets( line, sizeof(line), f ) ) {

				fprintf( flog,
				"ReadThmPair: Empty file [%s].\n", name );
				goto exit;
			}

			while( fgets( line, sizeof(line), f ) ) {

				ThmPair		P;

				sscanf( line, "%d\t%d\t%d\t%d\t%d"
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


