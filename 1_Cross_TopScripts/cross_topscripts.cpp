//
// Write scripts governing cross layer alignment.
//


#include	"Cmdline.h"
#include	"Disk.h"
#include	"File.h"
#include	"PipeFiles.h"
#include	"TrakEM2_UTL.h"


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* CArgs_cross --------------------------------------------------- */
/* --------------------------------------------------------------- */

class CArgs_cross {

public:
    char		srcmons[2048];
    const char	*script;
    int			zmin,
                zmax;

public:
    CArgs_cross()
    {
        script	= NULL;
        zmin	= 0;
        zmax	= 32768;

        strcpy( srcmons, "X_A_BIN_mons" );
    };

    void SetCmdLine( int argc, char* argv[] );
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static CArgs_cross	gArgs;
static ScriptParams	scr;
static string		idb;
static const char	*gtopdir = "cross_wkspc";
static FILE*		flog = NULL;






/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CArgs_cross::SetCmdLine( int argc, char* argv[] )
{
// start log

    flog = FileOpenOrDie( "cross_topscripts.log", "w" );

// log start time

    time_t	t0 = time( NULL );
    char	atime[32];

    strcpy( atime, ctime( &t0 ) );
    atime[24] = '\0';	// remove the newline

    fprintf( flog, "Make scapeops scripts: %s ", atime );

// parse command line args

    if( argc < 4 ) {
        printf(
        "Usage: cross_topscripts srcmons"
        " -script=scriptpath -z=i,j.\n" );
        exit( 42 );
    }

    vector<int>	vi;

    for( int i = 1; i < argc; ++i ) {

        const char	*_outdir;

        // echo to log
        fprintf( flog, "%s ", argv[i] );

        if( argv[i][0] != '-' )
            DskAbsPath( srcmons, sizeof(srcmons), argv[i], flog );
        else if( GetArgStr( script, "-script=", argv[i] ) )
            ;
        else if( GetArgList( vi, "-z=", argv[i] ) ) {

            if( 2 == vi.size() ) {
                zmin = vi[0];
                zmax = vi[1];
            }
            else {
                fprintf( flog,
                "Bad format in -z [%s].\n", argv[i] );
                exit( 42 );
            }
        }
        else {
            printf( "Did not understand option '%s'.\n", argv[i] );
            exit( 42 );
        }
    }

    fprintf( flog, "\n\n" );

    if( !DskExists( srcmons ) ) {

        fprintf( flog,
        "Can't find [%s].\n"
        "(Did you run gathermons.sht yet?)\n", srcmons );

        exit( 42 );
    }

    fflush( flog );
}

/* --------------------------------------------------------------- */
/* GetZList ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static void GetZList( vector<int> &zlist )
{
    const char	*name = FileNamePtr( gArgs.srcmons );
    char		path[2048];

    if( strstr( name, "X_A_BIN" ) ) {

        for( int z = gArgs.zmin; z <= gArgs.zmax; ++z ) {

            sprintf( path, "%s/X_A_%d.bin", gArgs.srcmons, z );

            if( DskExists( path ) )
                zlist.push_back( z );
        }
    }
    else if( strstr( name, "X_A_TXT" ) ) {

        for( int z = gArgs.zmin; z <= gArgs.zmax; ++z ) {

            sprintf( path, "%s/X_A_%d.txt", gArgs.srcmons, z );

            if( DskExists( path ) )
                zlist.push_back( z );
        }
    }
    else {
        fprintf( flog,
        "Unsupported montage folder type [%s].\n", gArgs.srcmons );
        exit( 42 );
    }
}

/* --------------------------------------------------------------- */
/* CreateTopDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void CreateTopDir()
{
// create top subdir
    DskCreateDir( gtopdir, flog );
}

/* --------------------------------------------------------------- */
/* WriteSubscapes ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void WriteSubscapes( vector<int> &zlist )
{
// compose common argument string

    char	sopt[2048];

    sprintf( sopt,
    "'%s' -script=%s -idb=%s -mb",
    gArgs.srcmons, gArgs.script, idb.c_str() );

// open file

    char	path[2048];
    FILE	*f;

    sprintf( path, "%s/subscapes.sht", gtopdir );
    f = FileOpenOrDie( path, "w", flog );

// header

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# First step in cross-layer alignment.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > scapeops srcmons -script=scriptpath -idb=idbpath [options]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Required:\n" );
    fprintf( f, "# srcmons\t\t\t\t;collected independent montages\n" );
    fprintf( f, "# -script=scriptpath\t;alignment pipeline params file\n" );
    fprintf( f, "# -idb=idbpath\t\t\t;path to idb directory\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Scapeops does montage drawing and/or strip aligning as follows:\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# If drawing a montage...\n" );
    fprintf( f, "#\n" );
    fprintf( f, "#\t-mb -zb=%%d\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# If aligning strips...\n" );
    fprintf( f, "#\n" );
    fprintf( f, "#\t-ab -za=%%d -zb=%%d\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Options:\n" );
    fprintf( f, "# -mb\t\t\t;make montage for layer zb\n" );
    fprintf( f, "# -ab\t\t\t;align layer za to zb\n" );
    fprintf( f, "# -za\t\t\t;layer za used only with -ab option\n" );
    fprintf( f, "# -zb\t\t\t;required layer zb\n" );
    fprintf( f, "# -abdbg\t\t;make diagnostic strip images and exit\n" );
    fprintf( f, "# -abdbgfull\t;make diagnostic full images and exit\n" );
    fprintf( f, "# -abctr=0\t\t;debug at this a-to-b angle\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Create output subdirs\n" );
    fprintf( f, "mkdir -p strips\n" );
    fprintf( f, "mkdir -p montages\n" );
    fprintf( f, "mkdir -p scplogs\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Submit layer pairs\n" );

// write all but last layer

    int	nz = zlist.size();

    for( int iz = 1; iz < nz; ++iz ) {

        fprintf( f,
        "QSUB_1NODE.sht 5 \"sc-%d\" \"-j y -o out.txt\" %d"
        " \"scapeops %s -ab -za=%d -zb=%d\"\n",
        zlist[iz - 1], scr.stripslots,
        sopt, zlist[iz], zlist[iz - 1] );
    }

// last layer

    fprintf( f, "\n" );
    fprintf( f, "# Just montage last layer\n" );

    fprintf( f,
    "QSUB_1NODE.sht 6 \"sc-%d\" \"-j y -o out.txt\" %d"
    " \"scapeops %s -zb=%d\"\n",
    zlist[nz - 1], scr.stripslots,
    sopt, zlist[nz - 1] );

    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteLowresgo ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteLowresgo()
{
    char	path[2048];
    FILE	*f;

    sprintf( path, "%s/lowresgo.sht", gtopdir );
    f = FileOpenOrDie( path, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Second step in cross-layer alignment.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Collects output from scapeops into 'LowRes.xml' stack\n" );
    fprintf( f, "# having one reduced scale montage per layer.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# You must subsequently view and edit LowRes.xml using TrakEM2.\n" );
    fprintf( f, "# Correct mistakes using 'Align with Manual Landmarks' feature\n" );
    fprintf( f, "# and then Save result with the same name.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > cross_lowres -z=i,j [options]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Required:\n" );
    fprintf( f, "# -z=i,j\t\t;assemble layers in range z=[i..j]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Options:\n" );
    fprintf( f, "# -table\t\t;write debugging striptable.txt and exit\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "cross_lowres -z=%d,%d\n",
    gArgs.zmin, gArgs.zmax );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* WriteScafgo --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteScafgo()
{
    char	path[2048];
    FILE	*f;

    sprintf( path, "%s/scafgo.sht", gtopdir );
    f = FileOpenOrDie( path, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Third step in cross-layer alignment.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Apply the coarse layer-layer tforms from 'LowRes.xml' to the\n" );
    fprintf( f, "# individual tiles in srcmons. The new coarsely aligned stack\n" );
    fprintf( f, "# will be created in cross_wkspc and will be named X_A_BIN_scaf\n" );
    fprintf( f, "# or X_A_TXT_scaf following the type of srcmons. The scaffold\n" );
    fprintf( f, "# serves both as the input for the block-block alignment, and\n" );
    fprintf( f, "# as the starting guess for the final LSQ alignment.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > cross_scaffold srcmons -z=i,j [options]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Required:\n" );
    fprintf( f, "# srcmons\t\t\t\t;collected independent montages\n" );
    fprintf( f, "# -z=i,j\t\t\t\t;align layers in range z=[i..j]\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Options:\n" );
    fprintf( f, "# -lowres=LowRes.xml\t;alternate coarse alignment reference\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "cross_scaffold %s -z=%d,%d\n",
    gArgs.srcmons, gArgs.zmin, gArgs.zmax );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* BuildScafPath ------------------------------------------------- */
/* --------------------------------------------------------------- */

static char* BuildScafPath( char srcscaf[2048] )
{
    DskAbsPath( srcscaf, 2048, gtopdir, flog );

    if( strstr( FileNamePtr( gArgs.srcmons ), "X_A_BIN" ) )
        strcat( srcscaf, "/X_A_BIN_scaf" );
    else
        strcat( srcscaf, "/X_A_TXT_scaf" );

    return srcscaf;
}

/* --------------------------------------------------------------- */
/* WriteCarvego -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void WriteCarvego()
{
    char	path[2048], srcscaf[2048];
    FILE	*f;

    sprintf( path, "%s/carvego.sht", gtopdir );
    f = FileOpenOrDie( path, "w", flog );

    fprintf( f, "#!/bin/sh\n" );
    fprintf( f, "\n" );
    fprintf( f, "# Purpose:\n" );
    fprintf( f, "# Fourth step in cross-layer alignment.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Carve each scaffold layer into blocks of size crossblocksize\n" );
    fprintf( f, "# and create new script 'bsub.sht' to distribute block-block\n" );
    fprintf( f, "# alignment jobs to cluster.\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# > cross_carveblocks srcscaf -script=scriptpath -z=i,j\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Required:\n" );
    fprintf( f, "# srcscaf\t\t\t\t;scaffold created by scafgo.sht\n" );
    fprintf( f, "# -script=scriptpath\t;alignment pipeline params file\n" );
    fprintf( f, "# -z=i,j\t\t\t\t;align layers in range z=[i..j]\n" );
    fprintf( f, "\n" );
    fprintf( f, "\n" );
    fprintf( f, "cross_carveblocks %s -script=%s -z=%d,%d\n",
    BuildScafPath( srcscaf ), gArgs.script,
    gArgs.zmin, gArgs.zmax );
    fprintf( f, "\n" );

    fclose( f );
    FileScriptPerms( path );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main( int argc, char* argv[] )
{
    vector<int>	zlist;

/* ------------------ */
/* Parse command line */
/* ------------------ */

    gArgs.SetCmdLine( argc, argv );

    IDBFromTemp( idb, ".", flog );

    if( idb.empty() )
        exit( 42 );

    if( !ReadScriptParams( scr, gArgs.script, flog ) )
        exit( 42 );

/* ---------------- */
/* Read source data */
/* ---------------- */

    GetZList( zlist );

    if( zlist.size() < 2 ) {
        fprintf( flog, "Fewer than 2 layers -- do nothing.\n" );
        goto exit;
    }

/* -------------- */
/* Create content */
/* -------------- */

    CreateTopDir();

    WriteSubscapes( zlist );
    WriteLowresgo();
    WriteScafgo();
    WriteCarvego();

/* ---- */
/* Done */
/* ---- */

exit:
    fprintf( flog, "\n" );
    fclose( flog );

    return 0;
}


