

#include	"ImageIO.h"
#include	"Maths.h"
#include	"Geometry.h"
#include	"Timer.h"





static double CRS(
    const vector<double>	&I,
    int						w,
    const IBox				&box,
    double					delx,
    double					dely )
{
    double	AB = 0.0;
    int		n = (box.T-box.B+1) * (box.R-box.L+1);

    for( int y = box.B; y <= box.T; ++y ) {

        for( int x = box.L; x <= box.R; ++x ) {

            double	a = I[x + w*y];
            double	b = InterpolatePixel( x + delx, y + dely, I, w );

            AB += a*b;
        }
    }

    return AB / n;
}


static double R(
    const vector<double>	&I,
    int						w,
    const IBox				&box,
    double					delx,
    double					dely )
{
    double	A = 0.0, B = 0.0, AA = 0.0, BB = 0.0, AB = 0.0;
    int		n = (box.T-box.B+1) * (box.R-box.L+1);

    for( int y = box.B; y <= box.T; ++y ) {

        for( int x = box.L; x <= box.R; ++x ) {

            double	a = I[x + w*y];
            double	b = InterpolatePixel( x + delx, y + dely, I, w );

            A	+= a;
            B	+= b;
            AA	+= a * a;
            BB	+= b * b;
            AB	+= a * b;
        }
    }

    return (n*AB - A*B) / sqrt( (n*AA - A*A) * (n*BB - B*B) );
}


static double DSQ(
    const vector<double>	&I,
    int						w,
    const IBox				&box,
    double					delx,
    double					dely )
{
    double	D = 0.0;
    int		n = (box.T-box.B+1) * (box.R-box.L+1);

    for( int y = box.B; y <= box.T; ++y ) {

        for( int x = box.L; x <= box.R; ++x ) {

            double	a = I[x + w*y];
            double	b = InterpolatePixel( x + delx, y + dely, I, w );

            b -= a;

            D += b*b;
        }
    }

    return D / n;
}


static void RunCorrTest()
{
    uint32	w, h;
    uint8*	I8 = Raster8FromAny(
            "/groups/apig/tomo/EX2/TIF1/A1_13_2.tif", w, h );
    int		nI = w * h;

    vector<double>	I( nI );

    for( int i = 0; i < nI; ++i )
        I[i] = I8[i];

    Normalize( I );

    IBox	B;
    double	delx;
    FILE	*f = fopen( "/groups/apig/tomo/EX2/corr.txt", "w" );

    fprintf( f, "delx\tcross\tr\tdifsq\n" );

    B.L	= B.B = 500 - 50;
    B.R	= B.T = 500 + 50;

    for( delx = -10.0; delx <= 10.0; delx += 0.1 ) {

        double	crs, r, dsq;

        crs	= CRS( I, w, B, delx, 10 );
        r	= R( I, w, B, delx, 10 );
        dsq	= DSQ( I, w, B, delx, 10 );

        fprintf( f, "%f\t%f\t%f\t%f\n", delx, crs, r, dsq );
    }

    fclose( f );
    RasterFree( I8 );
}


#if 0
// Convert SuperLU to Matrix Market Coordinate format
//
#include "File.h"
static void Splitfile()
{
const char *inpath = "/nobackup/bock/karsh/ZZ200/temp0/stackTest/triples-200";

FILE *fi = FileOpenOrDie( inpath, "r" );
FILE *fa = FileOpenOrDie( "mmtriples", "w" );
FILE *fb = FileOpenOrDie( "mmtriples.rhs", "w" );
CLineScan	LS;

long	n, nnz, nnzL = 0;
LS.Get( fi );
sscanf( LS.line, "%ld%ld", &n, &nnz );
for( long i = 0; i < nnz; ++i ) {
long	r, c;
LS.Get( fi );
sscanf( LS.line, "%ld%ld", &r, &c );
if( r >= c ) ++nnzL;
}
fclose( fi );


fi = FileOpenOrDie( inpath, "r" );
LS.Get( fi );
fprintf( fa, "%%%%MatrixMarket matrix coordinate real symmetric\n" );
fprintf( fa, "%ld %ld %ld\n", n, n, nnzL );

for( long i = 0; i < nnz; ++i ) {
double	v;
long	r, c;
LS.Get( fi );
sscanf( LS.line, "%ld%ld%lf", &r, &c, &v );
if( r >= c )
    fprintf( fa, "%ld %ld %.16f\n", r, c, v );
}

for( long i = 0; i < n; ++i ) {
double	v;
LS.Get( fi );
sscanf( LS.line, "%lf", &v );
fprintf( fb, "%.16f\n", v );
}

fclose( fi );
fclose( fa );
fclose( fb );
}
#endif


int main( int argc, char **argv )
{
    clock_t	t0 = StartTiming();

    RunCorrTest();

    StopTiming( stdout, "exper", t0 );

    return 0;
}


