

#include	"lsq_Globals.h"
#include	"lsq_LoadPoints.h"

#include	"File.h"
#include	"Timer.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* InitTablesToMaximum ------------------------------------------- */
/* --------------------------------------------------------------- */

static void InitTablesToMaximum()
{
	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		int	z = vL[iL].z;

		mZ[z] = iL;
		vR.push_back( Rgns( z ) );
	}
}

/* --------------------------------------------------------------- */
/* LoadPoints ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static void LoadPoints(
	const char	*tempdir,
	int			z,
	int			SorD,
	int			xhi,
	int			yhi )
{
	if( xhi < 0 )
		return;

	for( int y = 0; y <= yhi; ++y ) {

		for( int x = 0; x <= xhi; ++x ) {

			char	buf[2048];
			sprintf( buf, "%s/%d/%c%d_%d/pts.%s",
			tempdir, z, SorD, x, y,
			(SorD == 'S' ? "same" : "down") );

			FILE	*f = fopen( buf, "r" );

			if( f ) {

				double		x1, y1,
							x2, y2;
				int			z1, d1, r1,
							z2, d2, r2;

				while( 10 == fscanf( f, "CPOINT2"
					" %d.%d:%d %lf %lf"
					" %d.%d:%d %lf %lf\n",
					&z1, &d1, &r1, &x1, &y1,
					&z2, &d2, &r2, &x2, &y2 ) ) {

				}

				fclose( f );
			}
		}
	}
}

/* --------------------------------------------------------------- */
/* Destruct ------------------------------------------------------ */
/* --------------------------------------------------------------- */

CLoadPoints::~CLoadPoints()
{
}

/* --------------------------------------------------------------- */
/* Load ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CLoadPoints::Load( const char *tempdir, bool downtoo )
{
	printf( "\n---- Loading points ----\n" );

	clock_t	t0 = StartTiming();

	InitTablesToMaximum();

// Load the points

	int	nL = vL.size();

	for( int iL = 0; iL < nL; ++iL ) {

		const Layer&	L = vL[iL];

		LoadPoints( tempdir, L.z, 'S', L.sx, L.sy );

		if( downtoo )
			LoadPoints( tempdir, L.z, 'D', L.dx, L.dy );
printf( "loaded %d\n", L.z );fflush(stdout);
if(iL == 9 )break;
	}

	StopTiming( stdout, "Load", t0 );
}


