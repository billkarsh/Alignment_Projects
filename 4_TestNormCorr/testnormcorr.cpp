

#include	"Correlation.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* SomeArea ------------------------------------------------------ */
/* --------------------------------------------------------------- */

static bool SomeArea( int sx, int sy, void *a )
{
	return sx > 1 && sy > 1;
}

/* --------------------------------------------------------------- */
/* TestNormCorr -------------------------------------------------- */
/* --------------------------------------------------------------- */

static void TestNormCorr()
{
	vector<Point>	pts;
	vector<double>	vals;
	vector<Point>	pb;
	vector<double>	vb;

	pb.push_back(Point(0, 3)); vb.push_back(1);
	pb.push_back(Point(2, 3)); vb.push_back(2);
	pb.push_back(Point(3, 2)); vb.push_back(1);
	pb.push_back(Point(3, 3)); vb.push_back(3);
	pb.push_back(Point(3, 4)); vb.push_back(1);
	pb.push_back(Point(4, 3)); vb.push_back(2);
	pb.push_back(Point(5, 4)); vb.push_back(1);

	// for our example, want w!=h, to test code
	// exact match with delta y = 2, delta x = -4

	pts.push_back(Point(7.0, 2.0));
	pts.push_back(Point(8.0, 2.0));
	pts.push_back(Point(9.0, 2.0));
	pts.push_back(Point(7.0, 1.0));
	pts.push_back(Point(8.0, 1.0));

	vals.push_back(1.0);
	vals.push_back(0.0);
	vals.push_back(1.0);
	vals.push_back(3.0);
	vals.push_back(2.0);

	double		dx, dy;
	vector<CD>	ftc;

	clock_t		t1 = StartTiming();

	double	c = CorrPatches(
					stdout, false, dx, dy,
					pts, vals, pb, vb, 0, 0, 4096,
					SomeArea, NULL, NULL, NULL, ftc );

	printf( "TestNormCorr: corr=%f, dx=%f dy=%f\n", c, dx, dy );

	StopTiming( stdout, "TestNormCorr", t1 );
}

/* --------------------------------------------------------------- */
/* main ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

int main(int argc, char* argv[])
{
	TestNormCorr();
	return 0;
}


