

#include	"CGBL_dmesh.h"
#include	"InSectionOverlap.h"

#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"CPicBase.h"

#include	<stdlib.h>


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Inside -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Tells if a point is inside a given rectangle
static bool Inside( int x, int y, int x1, int y1, int x2, int y2 )
{
    return x > x1 && x < x2 && y > y1 && y < y2;
}

/* --------------------------------------------------------------- */
/* InSectionLegal ------------------------------------------------ */
/* --------------------------------------------------------------- */

// This is the function that defines a legal overlap
// for this application.  nx and ny are the sizes of
// the overlap region.
//
static bool InSectionLegal( int nx, int ny, void *a )
{
    double	n = nx * ny;

// 30,000 total pixels, and at least one side is 1500 pixels long
    return n > 30000 && (nx >= 1500 || ny > 1500);
}

/* --------------------------------------------------------------- */
/* MakeLambda ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void MakeLambda(
    vector<double>	&l,
    Point			p,
    vector<Point>	&cpts,
    int				i0,
    int				i1,
    int				i2,
    int				i3 )
{
// 0 2
// 1 3
    double alpha = (p.x-cpts[i0].x)/(cpts[i2].x-cpts[i0].x);
    double  beta = (p.y-cpts[i0].y)/(cpts[i1].y-cpts[i0].y);

    l[i0] = (1-alpha)*(1-beta);
    l[i1] = (1-alpha)*   beta;
    l[i2] =    alpha *(1-beta);
    l[i3] =    alpha *   beta;
}

/* --------------------------------------------------------------- */
/* PointAvg ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static Point PointAvg(
    vector<Point>	&pts,
    int				i1,
    int				i2,
    int				i3,
    int				i4 )
{
    Point P( (pts[i1].x + pts[i2].x + pts[i3].x + pts[i4].x)/4,
             (pts[i1].y + pts[i2].y + pts[i3].y + pts[i4].y)/4 );

    return P;
}

/* --------------------------------------------------------------- */
/* InSectionOverlap ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Here we seek a mapping from image a to image b in the same layer.
// The basic idea is that they are are overlapping at the edge by
// some unknown amount.  Find a few exact correspondence points.
//
// Npts		- number of correspondence points found
// apts		- their coordinates in image 1
// bpts		- their coordinates in image 2
// px		- the image data (calculates with scaled versions)
// flog		- log text file
//
void InSectionOverlap(
    int				&Npts,
    double*			&apts,
    double*			&bpts,
    const PixPair	&px,
    FILE*			flog )
{
    const vector<double>&	av_aln	= *px.avs_aln;
    const vector<double>&	bv_aln	= *px.bvs_aln;
    int						w		= px.ws,
                            h		= px.hs;

// Find how each of these regions overlap at the edges,
// if indeed they do at all.

    const double	frame	= 0.10 ;		// size of the frame
    int				x1		= int(frame*w);	// coord of center hole
    int				x2		= w - 1 - x1;
    int				y1		= int(frame*h);
    int				y2		= h - 1 - y1;

    fprintf( flog,
    "Olap: w %d, h %d, x1 %d, y1 %d, x2, %d, y2 %d\n",
    w, h, x1, y1, x2, y2 );

// Collect normalized points from the frame (border) area
// of the images.

    vector<double>	image2( 4096 * 4096, 0.0 );
    vector<Point>	ap, bp;
    vector<double>	av, bv;
    MeanStd			ma, mb;
    double			meana, stda, meanb, stdb;

    for( int ix = 0; ix < w; ++ix ) {

        for( int iy = 0; iy < h; ++iy ) {

            if( !Inside( ix, iy, x1, y1, x2, y2 ) ) {

                ma.Element( av_aln[ix + w*iy] );
                mb.Element( bv_aln[ix + w*iy] );
            }
        }
    }

    ma.Stats( meana, stda );
    mb.Stats( meanb, stdb );

    fprintf( flog,
    "Olap: Frame size %d pixels, %f percent.\n",
    ma.HowMany(), ma.HowMany()*100.0 / (w*h) );

    fprintf( flog,
    "Olap: A frame mean = %f and std dev = %f\n",
    meana, stda );

    fprintf( flog,
    "Olap: B frame mean = %f and std dev = %f\n",
    meanb, stdb );

    for( int ix = 0; ix < w; ++ix ) {

        for( int iy = 0; iy < h; ++iy ) {

            if( !Inside( ix, iy, x1, y1, x2, y2 ) ) {

                Point	p( ix, iy );
                double	apx = (av_aln[ix + w*iy] - meana) / stda,
                        bpx = (bv_aln[ix + w*iy] - meanb) / stdb;

                ap.push_back( p );
                bp.push_back( p );

                av.push_back( apx );
                bv.push_back( bpx );

                image2[ix + 4096*iy] = bpx;
            }
        }
    }

// Now call the FFT routine to do the actual work....

    double		dx, dy;
    vector<CD>	ftc;
    double		c = CorrPatches(
                        flog, false, dx, dy,
                        ap, av, bp, bv, 0, 0, 4000,
                        InSectionLegal, NULL,
                        NULL, NULL, ftc );

    if( c < GBL.ctx.RTRSH ) {	// empirical

        fprintf( flog,
        "Olap: Can't find a good starting point"
        " - maximum c = %f\n", c );

        Npts = 0;
        return;
    }

    fprintf( flog, "Olap: dx, dy = %f %f\n", dx, dy );

    int	idx = int(floor(dx+0.5));	// round
    int	idy = int(floor(dy+0.5));

    IBox	OL1, OL2;
    BoxesFromShifts( OL1, OL2, w, w, w, w, idx, idy );

// here we assume the box is horizontal.

    vector<Point>			cpts( 8, Point(0.0, 0.0) );
    vector<vector<double> >	lambdas;
    vector<double>			spv;	// source pixel values
    bool					horizontal = OL1.R-OL1.L > OL1.T-OL1.B;
    double					xm1 = OL1.L + 0.25*(OL1.R-OL1.L);
    double					xm2 = OL1.L + 0.75*(OL1.R-OL1.L);
    double					ym1 = OL1.B + 0.25*(OL1.T-OL1.B);
    double					ym2 = OL1.B + 0.75*(OL1.T-OL1.B);

    if( horizontal ) {

        // Arrange the points like this:
        // 0    2         4    6
        // 1    3         5    7
        cpts[1] = Point(OL1.L, OL1.T);
        cpts[3] = Point(xm1,   OL1.T);
        cpts[5] = Point(xm2,   OL1.T);
        cpts[7] = Point(OL1.R, OL1.T);
        cpts[0] = Point(OL1.L, OL1.B);
        cpts[2] = Point(xm1,   OL1.B);
        cpts[4] = Point(xm2,   OL1.B);
        cpts[6] = Point(OL1.R, OL1.B);
    }
    else {  // box is more vertical

        // 0 1
        // 2 3
        // 4 5
        // 6 7
        cpts[6] = Point(OL1.L, OL1.T);
        cpts[7] = Point(OL1.R, OL1.T);
        cpts[4] = Point(OL1.L, ym2);
        cpts[5] = Point(OL1.R, ym2);
        cpts[2] = Point(OL1.L, ym1);
        cpts[3] = Point(OL1.R, ym1);
        cpts[0] = Point(OL1.L, OL1.B);
        cpts[1] = Point(OL1.R, OL1.B);
    }

    for( int iy = OL1.B; iy <= OL1.T; ++iy ) {

        for( int ix = OL1.L; ix <= OL1.R; ++ix ) {

            spv.push_back( av_aln[ix + w*iy] );
            vector<double>	l( 8, 0.0 );
            Point			p( ix, iy );

            if( horizontal ) {

                if( ix < xm1 )
                    MakeLambda( l, p, cpts, 0, 1, 2, 3 );
                else if( ix < xm2 )
                    MakeLambda( l, p, cpts, 2, 3, 4, 5 );
                else
                    MakeLambda( l, p, cpts, 4, 5, 6, 7 );
            }
            else {  // box is vertical

                if( iy < ym1 )
                    MakeLambda( l, p, cpts, 0, 2, 1, 3 );
                else if( iy < ym2 )
                    MakeLambda( l, p, cpts, 2, 4, 3, 5 );
                else
                    MakeLambda( l, p, cpts, 4, 6, 5, 7 );
            }

            lambdas.push_back( l );
        }
    }

    vector<Point>	before = cpts;	// start with 8 control points
    before.push_back( PointAvg(cpts, 0, 1, 2, 3) );	// add rect centers
    before.push_back( PointAvg(cpts, 2, 3, 4, 5) );
    before.push_back( PointAvg(cpts, 4, 5, 6, 7) );

// Now transfer the control points to image b's coordinate system.

    for( int i = 0; i < 8; ++i ) {

        cpts[i].x += dx;
        cpts[i].y += dy;
    }

// Now improve the correlation by gradient descent

    double	threshold = 0.5;

    c = ImproveControlPts(
            cpts, lambdas, spv, image2, 4096, 4096, flog,
            "in-section", GBL.ctx.RIT, threshold );

    if( c < threshold ) {

        fprintf( flog,
        "Olap: Correlation not good enough - %f\n", c );

        Npts = 0;
        return;
    }

    vector<Point>	after = cpts;
    after.push_back( PointAvg(cpts, 0, 1, 2, 3) );
    after.push_back( PointAvg(cpts, 2, 3, 4, 5) );
    after.push_back( PointAvg(cpts, 4, 5, 6, 7) );

// Keep only the points that map into both images.
// Source points will always be inside, by construction, but
// some of the target points may fall outside after detailed
// matching (rotation, for example).

    for( int i = 0; i < after.size(); ) {

        if( after[i].x >= 0 &&
            after[i].x <  w &&
            after[i].y >= 0 &&
            after[i].y <  w ) { // it's OK

            ++i;
        }
        else {	// image
            before.erase( before.begin() + i );
            after.erase(  after.begin()  + i );
        }
    }

// fill in points

    Npts = after.size();

    fprintf( flog, "Olap: %d matching points.\n", Npts );

    apts = (double*)malloc( 2*Npts * sizeof(double) );
    bpts = (double*)malloc( 2*Npts * sizeof(double) );

// scale back to original image size, if bigger than 2K.

    int		sc = px.scl;

    for( int i = 0; i < Npts; ++i ) {

        fprintf( flog,
        "Olap: Point (%7.2f %7.2f) in image 1"
        " maps to point (%7.2f %7.2f) in image 2,"
        " delta (%7.2f %7.2f).\n",
        sc*before[i].x, sc*before[i].y,
        sc*after[i].x,  sc*after[i].y,
        sc*before[i].x - sc*after[i].x,
        sc*before[i].y - sc*after[i].y );

        apts[2*i]		= sc*before[i].x;
        apts[2*i + 1]	= sc*before[i].y;

        bpts[2*i]		= sc*after[i].x;
        bpts[2*i + 1]	= sc*after[i].y;
    }
}


