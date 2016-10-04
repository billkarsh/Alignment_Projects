

#include	"lsq_MTrans.h"

#include	"TrakEM2_UTL.h"
#include	"File.h"


/* --------------------------------------------------------------- */
/* RotateAll ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::RotateAll(
    vector<double>	&X,
    double			degcw )
{
    return;
}

/* --------------------------------------------------------------- */
/* NewOriginAll -------------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::NewOriginAll(
    vector<double>	&X,
    double			xorg,
    double			yorg )
{
    int	nr	= vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        int	itr = vRgn[i].itr;

        if( itr < 0 )
            continue;

        itr *= NX;

        X[itr  ] -= xorg;
        X[itr+1] -= yorg;
    }
}

/* --------------------------------------------------------------- */
/* SolveSystem --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::SolveSystem( vector<double> &X, int nTr )
{
    int	nvars	= nTr * NX,
        nc		= vAllC.size();

    printf( "Trans: %d unknowns; %d constraints.\n", nvars, nc );

    vector<double> RHS( nvars, 0.0 );
    vector<LHSCol> LHS( nvars );

    X.resize( nvars );

// SetPointPairs

    for( int i = 0; i < nc; ++i ) {

        const Constraint &C = vAllC[i];

        if( !C.used || !C.inlier )
            continue;

        int	j  = vRgn[C.r1].itr * NX,
            k  = vRgn[C.r2].itr * NX;

        // T1(p1) - T2(p2) = 0

        double	v[2]  = {  1,  -1};
        int		i1[2] = {  j,   k},
                i2[2] = {j+1, k+1};

        AddConstraint( LHS, RHS, 2, i1, v, C.p2.x - C.p1.x );
        AddConstraint( LHS, RHS, 2, i2, v, C.p2.y - C.p1.y );
    }

// SetIdentityTForm

    double	one	= 1;
    int		j	= nTr;	// = nTr/2 * 2

    AddConstraint( LHS, RHS, 1, &j, &one, 0 );	j++;
    AddConstraint( LHS, RHS, 1, &j, &one, 0 );

// SolveFromPoints

    WriteSolveRead( X, LHS, RHS, "T-FrmPts", nproc, false );
}

/* --------------------------------------------------------------- */
/* WriteTransforms ----------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::WriteTransforms(
    const vector<double>	&X,
    int						bstrings,
    FILE					*FOUT )
{
    printf( "---- Write transforms ----\n" );

    FILE	*f   = FileOpenOrDie( "TAffineTable.txt", "w" );
    int		nr   = vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        if( I.itr < 0 )
            continue;

        int	j = I.itr * NX;

        fprintf( f, "%d\t%d\t%d\t1\t0\t%f\t0\t1\t%f\n",
        I.z, I.id, I.rgn,
        X[j], X[j+1] );

        if( !bstrings ) {

            fprintf( FOUT, "TAFFINE %d.%d-%d 1 0 %f 0 1 %f\n",
            I.z, I.id, I.rgn,
            X[j], X[j+1] );
        }
        else {
            fprintf( FOUT, "TRANSFORM '%s::%d' 1 0 %f 0 1 %f\n",
            I.GetName(), I.rgn,
            X[j], X[j+1] );
        }
    }

    fclose( f );

    printf(
    "Average magnitude=1, min=1, max=1, max/min=1\n\n" );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteTrakEM --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::WriteTrakEM(
    double					xmax,
    double					ymax,
    const vector<double>	&X,
    double					trim,
    int						xml_type,
    int						xml_min,
    int						xml_max )
{
    FILE	*f = FileOpenOrDie( "MultLayAff.xml", "w" );

    int	oid = 3;

    fprintf( f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" );

    TrakEM2WriteDTD( f );

    fprintf( f, "<trakem2>\n" );

    fprintf( f,
    "\t<project\n"
    "\t\tid=\"0\"\n"
    "\t\ttitle=\"Project\"\n"
    "\t\tmipmaps_folder=\"trakem2.mipmaps/\"\n"
    "\t\tn_mipmap_threads=\"8\"\n"
    "\t/>\n" );

    fprintf( f,
    "\t<t2_layer_set\n"
    "\t\toid=\"%d\"\n"
    "\t\ttransform=\"matrix(1,0,0,1,0,0)\"\n"
    "\t\ttitle=\"Top level\"\n"
    "\t\tlayer_width=\"%.2f\"\n"
    "\t\tlayer_height=\"%.2f\"\n"
    "\t>\n",
    oid++, xmax, ymax );

    int	prev	= -1;	// will be previously written layer
    int	offset	= int(2 * trim + 0.5);
    int	nr		= vRgn.size();

    for( int i = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        // skip unused tiles
        if( I.itr < 0 )
            continue;

        // changed layer
        if( (*zs)[i].z != prev ) {

            if( prev != -1 )
                fprintf( f, "\t\t</t2_layer>\n" );

            fprintf( f,
            "\t\t<t2_layer\n"
            "\t\t\toid=\"%d\"\n"
            "\t\t\tthickness=\"0\"\n"
            "\t\t\tz=\"%d\"\n"
            "\t\t>\n",
            oid++, (*zs)[i].z );

            prev = (*zs)[i].z;
        }

        const char	*path;
        char		title[128];
        DisplayStrings( title, path, I );

        // fix origin : undo trimming
        int		j = I.itr * NX;
        double	x_orig = trim + X[j];
        double	y_orig = trim + X[j+1];

        fprintf( f,
        "\t\t\t<t2_patch\n"
        "\t\t\t\toid=\"%d\"\n"
        "\t\t\t\twidth=\"%d\"\n"
        "\t\t\t\theight=\"%d\"\n"
        "\t\t\t\ttransform=\"matrix(1,0,0,1,%f,%f)\"\n"
        "\t\t\t\ttitle=\"%s\"\n"
        "\t\t\t\ttype=\"%d\"\n"
        "\t\t\t\tfile_path=\"%s\"\n"
        "\t\t\t\to_width=\"%d\"\n"
        "\t\t\t\to_height=\"%d\"\n",
        oid++, gW - offset, gH - offset,
        x_orig, y_orig,
        title, xml_type, path, gW - offset, gH - offset );

        if( xml_min < xml_max ) {

            fprintf( f,
            "\t\t\t\tmin=\"%d\"\n"
            "\t\t\t\tmax=\"%d\"\n"
            "\t\t\t/>\n",
            xml_min, xml_max );
        }
        else
            fprintf( f, "\t\t\t/>\n" );
    }

    if( nr > 0 )
        fprintf( f, "\t\t</t2_layer>\n" );

    fprintf( f, "\t</t2_layer_set>\n" );
    fprintf( f, "</trakem2>\n" );
    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* WriteJython --------------------------------------------------- */
/* --------------------------------------------------------------- */

void MTrans::WriteJython(
    const vector<double>	&X,
    double					trim,
    int						Ntr )
{
    FILE	*f = FileOpenOrDie( "JythonTransforms.txt", "w" );

    fprintf( f, "transforms = {\n" );

    int	nr = vRgn.size();

    for( int i = 0, itrf = 0; i < nr; ++i ) {

        const RGN&	I = vRgn[(*zs)[i].i];

        // skip unused tiles
        if( I.itr < 0 )
            continue;

        ++itrf;

        const char	*path;
        DisplayStrings( NULL, path, I );

        // fix origin : undo trimming
        int		j = I.itr * NX;
        double	x_orig = trim + X[j];
        double	y_orig = trim + X[j+1];

        fprintf( f, "\"%s\" : [1, 0, 0, 1, %f, %f]%s\n",
            path, x_orig, y_orig,
            (itrf == Ntr ? "" : ",") );
    }

    fprintf( f, "}\n" );
    fclose( f );

    IDBT2ICacheClear();
}

/* --------------------------------------------------------------- */
/* G2LPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MTrans::G2LPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    int	j = itr * NX;
    p.x -= X[j];
    p.y -= X[j+1];
}

/* --------------------------------------------------------------- */
/* L2GPoint ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void MTrans::L2GPoint(
    Point					&p,
    const vector<double>	&X,
    int						itr )
{
    int	j = itr * NX;
    p.x += X[j];
    p.y += X[j+1];
}


void MTrans::L2GPoint(
    vector<Point>			&p,
    const vector<double>	&X,
    int						itr )
{
    int	np = p.size(),
        j  = itr * NX;

    for( int i = 0; i < np; ++i ) {
        p[i].x += X[j];
        p[i].y += X[j+1];
    }
}


