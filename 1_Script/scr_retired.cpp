/* --------------------------------------------------------------- */
/* Make_MakePts -------------------------------------------------- */
/* --------------------------------------------------------------- */

// This writes a script telling program align to scan tiny and
// ptest output: {foldmasks, n/m.tf.txt, n/m.map.tif} and create
// text lists {pts.same, pts.down, ...} of FOLDMAP entries and
// POINT entries, which are triangle centroid->centroid mappings.
//
// However, align and this script are retired. Rather, tiny makes
// file fm.same (with FOLDMAPs) in each layer. Program ptest makes
// outputs {pts.same, pts.up, pts.down} according to the jobs
// submitted to it.
//
// Shell script 'combine' is still used to select among these
// and combine them into pts.all for the lsq program.
//
static void Make_MakePts( const char *lyrdir, int L )
{
	char	name[1024];
	FILE	*f;

	fprintf( flog, "--Make_MakePts: layer %d\n", L );

	sprintf( name, "%s/make.pts", lyrdir );

	f = FileOpenOrDie( name, "w", flog );

// master target depends on all others

	//fprintf( f, "all: pts.same pts.same.nf pts.up pts.down\n\n" );
	//fprintf( f, "all: pts.same pts.up pts.down\n\n" );
	fprintf( f, "all: pts.same pts.down\n\n" );

// subtargets and rules

    fprintf( f, "pts.same:\n" );
    fprintf( f, "\talign -fpts.same make.same\n\n");

    //fprintf( f, "pts.same.nf:\n" );
    //fprintf( f, "\talign -fpts.same.nf make.same -nf\n\n");

    //fprintf( f, "pts.up:\n" );
    //fprintf( f, "\talign -fpts.up make.up\n\n");

    fprintf( f, "pts.down:\n" );
    fprintf( f, "\talign -fpts.down make.down\n");

	fclose( f );
}


