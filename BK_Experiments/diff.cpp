

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"CTForm.h"
#include	"Timer.h"
#include	"CTemplate.h"
#include	"CPicBase.h"


typedef struct {
	char			path[256], name[64];
	uint32			w, h, n;
	uint8*			ras;
	vector<double>	V;
	vector<Point>	P;

	void LoadPic( const char *path );
} Pic;


typedef struct {
	Pic				a, b;
	vector<double>	Ib;

	void PrepImages( const char *patha, const char *pathb );
	void Release();
	void TestPatchToImage( double& dx, double &dy );
	void RefineDxyByXCor( double& dx, double &dy );
	void RefineDxyBySumSqDif( double& dx, double &dy );
	void WriteDiffImage( double dx, double dy );
} Raw;


typedef struct {
	vector<double>	av, bv;
	vector<Point>	apts, bpts;
	vector<CD>		ftc;	// fourier transform cache
	long			reqArea;
	int				scl;
} Thumbs;


typedef struct {
	TForm	T;
	double	X, Y,
			R, A;
} CorRec;






void Pic::LoadPic( const char *path )
{
	strcpy( this->path, path );

	const char *s = FileNamePtr( path );

	strncpy( name, s, strlen( s ) - 4 );

	ras	= Raster8FromTif( path, w, h, stdout );
//----------------
//{
//	if( strstr( name, "A1_3_2" ) ) {
//		h = h/2;
//	}
//}
//----------------
	n	= w * h;

	V.resize( n );
	P.resize( n );

	for( int i = 0; i < n; ++i ) {

		int		y = i / w;
		int		x = i - w*y;

		P[i]	= Point( x, y );
		V[i]	= ras[i];
	}
}


void Raw::PrepImages( const char *patha, const char *pathb )
{
	a.LoadPic( patha );
	b.LoadPic( pathb );

	Normalize( a.V );
	Normalize( b.V );

	uint32	szfft = FFTSizeSQR( a.w, a.h, b.w, b.h );
	Ib.resize( szfft*szfft, 0.0 );

	CopyRaster( &Ib[0], szfft, &b.V[0], b.w, b.w, b.h );
}


void Raw::Release()
{
	RasterFree( a.ras );
	RasterFree( b.ras );
}


void Raw::TestPatchToImage( double& dx, double &dy )
{
	CorrPatchToImage( dx, dy, a.P, a.V, Ib,
		0, 0, (int)sqrt( a.w*a.w + a.h*a.h ), false );
}


void Raw::RefineDxyByXCor( double& dx, double &dy )
{
	double	x0	= floor(dx) - 1, xlim = x0 + 2.0,
			y0	= floor(dy) - 1, ylim = y0 + 2.0;
	double	B	= 0.0;

	for( double	tx = x0; tx < xlim; tx += 0.05 ) {

		for( double	ty = y0; ty < ylim; ty += 0.05 ) {

			double sum = 0.0;

			for( int i = 0; i < a.n; ++i ) {

				int		y	= i / a.w;
				int		x	= i - a.w*y;
				double	xb	= x + tx;
				double	yb	= y + ty;

				if( xb >= 0 && xb < b.w && yb >= 0 && yb < b.h ) {

					sum += a.V[i] *
					InterpolatePixel( xb, yb, b.V, b.w );
				}
			}

			if( sum > B ) {
				B = sum;
				dx = tx;
				dy = ty;
			}
		}
	}

	printf( "\nBest by sumxcor= %f at (%f, %f).\n\n", B, dx, dy );
}


void Raw::RefineDxyBySumSqDif( double& dx, double &dy )
{
	double	x0	= floor(dx) - 1, xlim = x0 + 2.0,
			y0	= floor(dy) - 1, ylim = y0 + 2.0;
	double	B	= 1E30;

	for( double	tx = x0; tx < xlim; tx += 0.05 ) {

		for( double	ty = y0; ty < ylim; ty += 0.05 ) {

			double sum = 0.0;

			for( int i = 0; i < a.n; ++i ) {

				int		y	= i / a.w;
				int		x	= i - a.w*y;
				double	xb	= x + tx;
				double	yb	= y + ty;

				if( xb >= 0 && xb < b.w && yb >= 0 && yb < b.h ) {

					double	dif = a.V[i] -
					InterpolatePixel( xb, yb, b.V, b.w );

					sum += dif * dif;
				}
			}

			if( sum < B ) {
				B = sum;
				dx = tx;
				dy = ty;
			}
		}
	}

	printf( "\nBest by sumsqdif= %f at (%f, %f).\n\n", B, dx, dy );
}


void Raw::WriteDiffImage( double dx, double dy )
{
// we work in the frame of image a
	uint8*	dif = (uint8*)RasterAlloc( a.n * sizeof(uint8) );

	for( int i = 0; i < a.n; ++i ) {

		int		y	= i / a.w;
		int		x	= i - a.w*y;
		double	xb	= x + dx;
		double	yb	= y + dy;
		int		pix	= 127;

		if( xb >= 0 && xb < b.w && yb >= 0 && yb < b.h ) {

			pix += int(a.ras[i] -
			InterpolatePixel( xb, yb, b.ras, b.w ));

			if( pix < 0 )
				pix = 0;
			else if( pix > 255 )
				pix = 255;
		}

		dif[i] = pix;
	}

	Raster8ToTif8( "diff.tif", dif, a.w, a.h );

	RasterFree( dif );
}


static void MakeThumbs(
	Thumbs				&thm,
	Raw					&R,
	int					order,
	int					reqArea,
	double				nsig )
{
	LegPolyFlatten( thm.av, R.a.P, R.a.ras,
		R.a.w, R.a.h, order );

	LegPolyFlatten( thm.bv, R.b.P, R.b.ras,
		R.b.w, R.b.h, order );

// Create thunmbnail images reduced in area by 8x8=64.
// Keeping the remaining #pixels in range [20,000; 160,000]
// works well (experimentally).

	const int linbin = 8;

	thm.ftc.clear();
	thm.reqArea	= reqArea / (linbin*linbin);
	thm.scl		= linbin;
	thm.apts	= R.a.P;
	thm.bpts	= R.b.P;

	DecimateVector( thm.apts, thm.av, R.a.w, R.a.h, linbin );
	DecimateVector( thm.bpts, thm.bv, R.b.w, R.b.h, linbin );

//	CoExcludeMiddle( thm.av, thm.bv, nsig );
//	CoNormalize( thm.av, thm.bv );
	Normalize( thm.av );
	Normalize( thm.bv );
}


static void RotatePoints(
	vector<Point>	&pts,
	TForm			&T,
	double			theta )
{
	theta *= PI/180.0;

	double	c	= cos( theta ),
			s	= sin( theta );
	TForm	t(	c, -s, 0.0,
				s,  c, 0.0 );

	T = t;

	t.Transform( pts );
}


static bool BigEnough( int sx, int sy, void *a )
{
	return (long)sx*sy > long(a);
}


static bool EnoughPoints( int count1, int count2, void *a )
{
	return count1 > long(a) && count2 > long(a);
}


static void RFromAngle(
	CorRec	&C,
	double	a,
	Thumbs	&thm )
{
	vector<Point>	ps = thm.apts;

	C.A = a;
	RotatePoints( ps, C.T, a );

	C.R = CorrImages(
		stdout, true, C.X, C.Y,
		ps, thm.av, thm.bpts, thm.bv,
		BigEnough, (void*)thm.reqArea,
		EnoughPoints, (void*)thm.reqArea,
		0.0, 0.75, thm.ftc );
}


static void BasicAScan( Raw& R )
{
	char	file[256];
	Thumbs	thm;

	sprintf( file, "angsR_%s_@_%s.log", R.a.name, R.b.name );
	FILE	*f = fopen( file, "w" );

	fprintf( f, "Deg\tR\tX\tY\tX\tY\n" );

	MakeThumbs( thm, R, 2, 25000, 0.0 );

	for( double a = -45; a <= 45; a += 0.5 ) {

		CorRec	C;

		RFromAngle( C, a, thm );

		Point	p( C.X * thm.scl, C.Y * thm.scl );

		fprintf( f, "%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
			a, C.R, C.X, C.Y, p.x, p.y );
	}

	fclose( f );
}


static void RvsSigma( Raw& R )
{
	char	file[256];
	Thumbs	thm;
	double	s0 =   0.0, sn =  2.0, sd = 0.5;
	double	a0 = -45.0, an = 45.0, ad = 0.5;
	vector<vector<double> >	r( int((sn - s0)/sd) + 1 );
	int						is = 0;

	for( double sig = s0; sig <= sn; sig += sd ) {

		r[is].resize( int((an - a0)/ad) + 1 );

		MakeThumbs( thm, R, 1, 200000, sig );

		int		ia = 0;
		for( double a = a0; a <= an; a += ad ) {

			CorRec	C;

			RFromAngle( C, a, thm );
			r[is][ia++] = C.R;
		}

		++is;
	}

	sprintf( file, "angsR_%s_@_%s.log", R.a.name, R.b.name );
	FILE	*f = fopen( file, "w" );

	for( int ia = 0; ia < r[0].size(); ++ia ) {

		fprintf( f, "%.3f\t", a0 + ia*ad );

		for( int is = 0; is < r.size(); ++is )
			fprintf( f, "%.4f\t", r[is][ia] );

		fprintf( f, "\n" );
	}

	fclose( f );
}


static void TestTemplate( Raw& R )
{
	char	buf[256];
	sprintf( buf, "out.txt" );
	freopen( buf, "w", stdout );

	PicBase		pa, pb;
	pa.SetExternal( R.a.ras, R.a.w, R.a.h );
	pb.SetExternal( R.b.ras, R.b.w, R.b.h );

	Template	M( pb, 0, 0, (int)R.b.w - 1, 0, (int)R.b.h - 1 );
	M.Match( pa, 1, 0, 0, (int)R.a.w - 1, (int)R.a.h - 1 );
}


static void JustDoOnePair( Raw& R )
{
	CorRec	C;
	Thumbs	thm;

	MakeThumbs( thm, R, 2, 25000, 0 );

	RFromAngle( C, 0, thm );
}


static void TestLoad( const char *patha, const char *pathb )
{
	uint32			w, h, n;
	uint8			*rasa, *rasb;

	rasa	= Raster8FromTif( patha, w, h, stdout );
	rasb	= Raster8FromTif( pathb, w, h, stdout );

	vector<double>	av, bv;
	vector<Point>	ap, bp;

	MakeZeroBasedPoints( ap, w, h );
	bp = ap;

	LegPolyFlatten( av, ap, rasa, w, h, 2 );
	LegPolyFlatten( bv, bp, rasb, w, h, 2 );

	RasterFree( rasa );
	RasterFree( rasb );
}


static void ReduceLeginonCoords()
{
	FILE	*fi = fopen( "120306201357_10x24_layout.txt", "r" );
	FILE	*fo = fopen( "120306201357_10x24_layout2.txt", "w" );

	char	path[2048];
	int		x, y, z;

	while( fscanf( fi, "%s%d%d%d", path, &x, &y, &z ) == 4 ) {
		fprintf( fo, "%s\t%d\t%d\t%d\n", path, x/4, y/4, z );
	}

	fclose( fi );
	fclose( fo );
}


int main( int argc, char* argv[] )
{
	Raw		R;
	double	dx, dy;

	R.PrepImages( argv[1], argv[2] );


clock_t	t0 = StartTiming();
//	BasicAScan( R );
TestLoad( argv[1], argv[2] );
StopTiming( stdout, "exper", t0 );


//	TestTemplate( R );


//	RvsSigma( R );

//	JustDoOnePair( R );



	goto exit;

// basic tests --------------------------------------------
//	R.TestPatchToImage( dx, dy );

//	R.RefineDxyByXCor( dx, dy );

//	R.RefineDxyBySumSqDif( dx, dy );

//	R.WriteDiffImage( dx, dy );
// basic tests --------------------------------------------


exit:
	R.Release();
	return 0;
}


