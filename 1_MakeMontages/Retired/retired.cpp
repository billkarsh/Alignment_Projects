

/* --------------------------------------------------------------- */
/* WriteSFinishFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSFinishFile()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/sfinish.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Run finish.sht script on its own cluster node.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > ./sfinish.sht\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "QSUB_1NODE.sht 9 \"finish\" \"\" 1 8 \"./finish.sht\"\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteFinishFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteFinishFile()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/finish.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Gather all points to stack folder, cd there, run lsq solver.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# ./finish.sht\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "./combine.sht %d %d\n",
    gArgs.zmin, gArgs.zmax );
    fprintf( f, "cd stack\n" );
    fprintf( f, "./runlsq.sht \"\"\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteCombineFile ---------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCombineFile()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/combine.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Gather all point pair files into stack/pts.all\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > ./combine.sht <zmin> <zmax>\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "rm -f pts.all\n" );
    fprintf( f, "\n" );
    fprintf( f, "#get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
    fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < imageparams.txt)\n" );
    fprintf( f, "\n" );
    fprintf( f, "cp imageparams.txt pts.all\n" );
    fprintf( f, "\n" );
    fprintf( f, "for lyr in $(seq $1 $2)\n" );
    fprintf( f, "do\n" );
    fprintf( f, "\tcat $idb/$lyr/fm.same >> pts.all\n" );
    fprintf( f, "done\n" );
    fprintf( f, "\n" );
    fprintf( f, "for lyr in $(seq $1 $2)\n" );
    fprintf( f, "do\n" );
    fprintf( f, "\techo $lyr\n" );
    fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
    fprintf( f, "\tthen\n" );
    fprintf( f, "\t\tcd $lyr\n" );
    fprintf( f, "\n" );
    fprintf( f, "\t\tfor jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
    fprintf( f, "\t\tdo\n" );
    fprintf( f, "\t\t\tcat $jb/pts.same >> ../pts.all\n" );
    fprintf( f, "\t\tdone\n" );
    fprintf( f, "\n" );
    fprintf( f, "\t\tif (($lyr != $1))\n" );
    fprintf( f, "\t\tthen\n" );
    fprintf( f, "\t\t\tfor jb in $(ls -d * | grep -E 'D[0-9]{1,}_[0-9]{1,}')\n" );
    fprintf( f, "\t\t\tdo\n" );
    fprintf( f, "\t\t\t\tcat $jb/pts.down >> ../pts.all\n" );
    fprintf( f, "\t\t\tdone\n" );
    fprintf( f, "\t\tfi\n" );
    fprintf( f, "\n" );
    fprintf( f, "\t\tcd ..\n" );
    fprintf( f, "\tfi\n" );
    fprintf( f, "done\n" );
    fprintf( f, "\n" );
    fprintf( f, "mv pts.all stack\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteMontage1File --------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteMontage1File()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/montage1.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# For given layer, gather all point pair files into montage/pts.all,\n" );
    fprintf( f, "# cd there, run lsq solver.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > ./montage1.sht <z>\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "cd $1\n" );
    fprintf( f, "\n" );
    fprintf( f, "rm -f pts.all\n" );
    fprintf( f, "\n" );
    fprintf( f, "# get line 1, subst 'IDBPATH=xxx' with 'xxx'\n" );
    fprintf( f, "idb=$(sed -n -e 's|IDBPATH \\(.*\\)|\\1|' -e '1p' < ../imageparams.txt)\n" );
    fprintf( f, "\n" );
    fprintf( f, "cp ../imageparams.txt pts.all\n" );
    fprintf( f, "cat $idb/$1/fm.same >> pts.all\n" );
    fprintf( f, "\n" );
    fprintf( f, "for jb in $(ls -d * | grep -E 'S[0-9]{1,}_[0-9]{1,}')\n" );
    fprintf( f, "do\n" );
    fprintf( f, "\tcat $jb/pts.same >> pts.all\n" );
    fprintf( f, "done\n" );
    fprintf( f, "\n" );
    fprintf( f, "mv pts.all montage\n" );
    fprintf( f, "\n" );
    fprintf( f, "cd montage\n" );
    fprintf( f, "./runlsq.sht \"\"\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmonFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmonFile()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/submon.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# For each layer in range, gather points, cd to layer's montage dir,\n" );
    fprintf( f, "# run lsq there.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > ./submon.sht <zmin> [zmax]\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "if (($# == 1))\n" );
    fprintf( f, "then\n" );
    fprintf( f, "\tlast=$1\n" );
    fprintf( f, "else\n" );
    fprintf( f, "\tlast=$2\n" );
    fprintf( f, "fi\n" );
    fprintf( f, "\n" );
    fprintf( f, "for lyr in $(seq $1 $last)\n" );
    fprintf( f, "do\n" );
    fprintf( f, "\techo $lyr\n" );
    fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
    fprintf( f, "\tthen\n" );
    fprintf( f, "\t\tQSUB_1NODE.sht 4 \"mon-$lyr\" \"\" 1 8 \"./montage1.sht $lyr\"\n" );
    fprintf( f, "\tfi\n" );
    fprintf( f, "done\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

/* --------------------------------------------------------------- */
/* WriteSubmosFile ----------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteSubmosFile()
{
    char	buf[2048];
    FILE	*f;

    sprintf( buf, "%s/mosaic/submos.sht", gArgs.outdir );
    f = FileOpenOrDie( buf, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Heal montage seams and update superpixel data.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > mos <simple-file> <xmin,ymin,xsize,ysize> <zmin,zmax> [options]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Default sizing is 0,0,-1,-1 meaning natural size.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Options:\n" );
    fprintf( f, "# -d\t\t\t\t;debug\n" );
    fprintf( f, "# -strings\t\t\t;simple-file data labeled by path strings\n" );
    fprintf( f, "# -warp\t\t\t\t;heal seams\n" );
    fprintf( f, "# -nf\t\t\t\t;no folds\n" );
    fprintf( f, "# -a\t\t\t\t;annotate\n" );
    fprintf( f, "# -tiles\t\t\t;make raveler tiles\n" );
    fprintf( f, "# -noflat\t\t\t;no flat image ('before')\n" );
    fprintf( f, "# -nomap\t\t\t;no map image (where data from)\n" );
    fprintf( f, "# -matlab\t\t\t;matlab/closeness order\n" );
    fprintf( f, "# -drn\t\t\t\t;don't renumber superpixels\n" );
    fprintf( f, "# -dms=0.01\t\t\t;don't move strength\n" );
    fprintf( f, "# -fold_dir=path\t;prepended fm location, default=CWD\n" );
    fprintf( f, "# -region_dir=path\t;results go here, default=CWD\n" );
    fprintf( f, "# -gray_dir=path\t;gray images go here, default=CWD\n" );
    fprintf( f, "# -grey_dir=path\t;gray images go here, default=CWD\n" );
    fprintf( f, "# -sp_dir=path\t\t;superpixel maps go here, default=CWD\n" );
    fprintf( f, "# -inv_dir=path\t\t;inverse maps go here, default=CWD\n" );
    fprintf( f, "# -rav_dir=path\t\t;raveler tiles go here, default=CWD\n" );
    fprintf( f, "# -bmap_dir=path\t;boundary maps go here, default=CWD\n" );
    fprintf( f, "# -s=1\t\t\t\t;scale down by this integer\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "export MRC_TRIM=12\n" );
    fprintf( f, "\n" );
    fprintf( f, "if (($# == 1))\n" );
    fprintf( f, "then\n" );
    fprintf( f, "\tlast=$1\n" );
    fprintf( f, "else\n" );
    fprintf( f, "\tlast=$2\n" );
    fprintf( f, "fi\n" );
    fprintf( f, "\n" );
    fprintf( f, "for lyr in $(seq $1 $last)\n" );
    fprintf( f, "do\n" );
    fprintf( f, "\techo $lyr\n" );
    fprintf( f, "\tif [ -d \"$lyr\" ]\n" );
    fprintf( f, "\tthen\n" );
    fprintf( f, "\t\tQSUB_1NODE.sht 10 \"mos-$lyr\" \"\" 1 8 \"mos ../stack/simple 0,0,-1,-1 $lyr,$lyr -warp%s > mos_$lyr.txt\"\n",
    (gArgs.NoFolds ? " -nf" : "") );
    fprintf( f, "\tfi\n" );
    fprintf( f, "done\n" );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( buf );
}

