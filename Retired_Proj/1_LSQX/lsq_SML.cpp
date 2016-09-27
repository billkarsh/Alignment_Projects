

#include	"lsq_Types.h"
#include	"lsq_SML.h"

#include	"LinEqu.h"

#include	<math.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* CanAlign ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return RMS error assuming that a similarity transformation
// (rot + trans) takes points {p1} into corr. points {p2}:
//
//		a x  -  b y  +  c  =  x'
//		b x  +  a y  +  d  =  y'
//
double SML::CanAlign(
	const vector<Point>	&p1,
	const vector<Point>	&p2,
	bool				print )
{
// Create system of normal equations

	double	RHS[4], *X = RHS;
	double	LHS[4*4];
	int		np    = p1.size(),
			i1[3] = { 0, 1, 2 },
			i2[3] = { 1, 0, 3 };

	Zero_Quick( LHS, RHS, 4 );

	for( int i = 0; i < np; ++i ) {

		double	v1[3] = { p1[i].x, -p1[i].y, 1.0 };
		double	v2[3] = { p1[i].x,  p1[i].y, 1.0 };

		AddConstraint_Quick( LHS, RHS, 4, 3, i1, v1, p2[i].x );
		AddConstraint_Quick( LHS, RHS, 4, 3, i2, v2, p2[i].y );
	}

// Solve

	Solve_Quick( LHS, RHS, 4 );

	if( print ) {

		double	mag = sqrt( X[0]*X[0] + X[1]*X[1] );

		printf( "  a=%f b=%f x0=%f y0=%f mag=%f\n",
		X[0], X[1], X[2], X[3], mag );
	}

// Calc error

	double	sum = 0.0;

	for( int i = 0; i < np; ++i ) {

		double px = p1[i].x * X[0] - p1[i].y * X[1] + X[2];
		double py = p1[i].x * X[1] + p1[i].y * X[0] + X[3];
		double dx = px - p2[i].x;
		double dy = py - p2[i].y;

		sum += dx*dx + dy*dy;

		if( print ) {
			printf(
			"  (%8.2f %8.2f) -> (%8.2f %8.2f) =?"
			" (%8.2f %8.2f) d=%8.2f\n",
			p1[i].x, p1[i].y, px, py,
			p2[i].x, p2[i].y, sqrt(dx*dx + dy*dy) );
		}
	}

	double	rms = sqrt( sum / np );

	if( print )
		printf( "  RMS error %f pixels.\n", rms );

	return rms;
}

/* --------------------------------------------------------------- */
/* AddPOINTPair -------------------------------------------------- */
/* --------------------------------------------------------------- */

// Add new pair of correspondence points to r12Pts table.
//
void SML::AddPOINTPair(
	int			r1,
	const Point	&p1,
	int			r2,
	const Point	&p2 )
{
	CRPair						r12( r1, r2 );
	map<CRPair,int>::iterator	mi = r12Idx.find( r12 );
	int							loc;

	if( mi == r12Idx.end() ) {

		// did not exist; add it
		r12Pts.push_back( CorrPts() );
		r12Idx[r12] = loc = r12Pts.size() - 1;
	}
	else
		loc = mi->second;

// Add points in correct order

	if( r1 < r2 ) {
		r12Pts[loc].pa.push_back( p1 );
		r12Pts[loc].pb.push_back( p2 );
	}
	else {
		r12Pts[loc].pa.push_back( p2 );
		r12Pts[loc].pb.push_back( p1 );
	}
}

/* --------------------------------------------------------------- */
/* TestPairAlignments -------------------------------------------- */
/* --------------------------------------------------------------- */

// For each pair of regions that we can (5 or more corr points)
// solve for a similarity transform. If the rms error exceeds tol.
// write a nasty report.
//
// There are no other material consequences to this step.
//
void SML::TestPairAlignments()
{
	printf(
	"---- Lyr.til-rgn pairs with simlr align err > 70 pix ----\n" );

	map<CRPair,int>::iterator	pi;	// iterator over pairs

	for( pi = r12Idx.begin(); pi != r12Idx.end(); ++pi ) {

		if( r12Pts[pi->second].pa.size() >= 5 ) {

			double	RMS = CanAlign(
							r12Pts[pi->second].pa,
							r12Pts[pi->second].pb, false );

			if( RMS > 70.0 ) {	// rough limit

				const RGN	&A = vRgn[(pi->first).a];
				const RGN	&B = vRgn[(pi->first).b];

				printf( "%d.%d-%d ^ %d.%d-%d\n",
				A.z, A.id, A.rgn, B.z, B.id, B.rgn );

				// repeat for report
				CanAlign(
					r12Pts[pi->second].pa,
					r12Pts[pi->second].pb, true );
			}
		}
	}

	printf( "\n" );
}


