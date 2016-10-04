

#include	"Draw.h"
#include	"Maths.h"

#include	<math.h>






/* --------------------------------------------------------------- */
/* DrawLine ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Draw a line in a raster.
//
void DrawLine(
    uint8*	r,
    int		w,
    int		h,
    double	x1,
    double	y1,
    double	x2,
    double	y2 )
{
    if( fabs( x1 - x2 ) > fabs( y1 - y2 ) ) {

        // mostly horizontal

        if( x1 > x2 ) {

            double	t;
            t = x1; x1 = x2; x2 = t;
            t = y1; y1 = y2; y2 = t;
        }

        for( int ix = ROUND( x1 ); ix <= ROUND( x2 ); ++ix ) {

            double	y		= y1 + (ix-x1)/(x2-x1)*(y2-y1);
            int		iy		= ROUND( y );

            if( 0 <= ix && ix < w && 0 <= iy && iy < h )
                r[ix + (uint32)w*iy] = 255;
        }
    }
    else {

        // mostly vertical

        if( y1 > y2 ) {

            double	t;
            t = x1; x1 = x2; x2 = t;
            t = y1; y1 = y2; y2 = t;
        }

        for( int iy = ROUND( y1 ); iy <= ROUND( y2 ); ++iy ) {

            double	x		= x1 + (iy-y1)/(y2-y1)*(x2-x1);
            int		ix		= ROUND( x );

            if( 0 <= ix && ix < w && 0 <= iy && iy < h )
                r[ix + (uint32)w*iy] = 255;
        }
    }
}

/* --------------------------------------------------------------- */
/* DrawCircle ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Draw a circle in a raster.
//
void DrawCircle(
    uint8*	r,
    int		w,
    int		h,
    double	xc,
    double	yc,
    double	radius )
{
    const double	dt = 2.0*PI/100.0;

    for( int i = 0; i < 100; ++i ) {

        double	theta = i * dt;
        int		ix = RND( xc + radius*cos( theta ) );
        int		iy = RND( yc + radius*sin( theta ) );

        if( 0 <= ix && ix < w && 0 <= iy && iy < h )
            r[ix + (uint32)w*iy] = 255;
    }
}


