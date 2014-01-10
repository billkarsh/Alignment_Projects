

#include	"lsq_Globals.h"
#include	"lsq_MPI.h"
#include	"lsq_Untwist.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"TAffine.h"
#include	"Timer.h"


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class RgdSums {
public:
	TAffine	A;
	double	Xa, Ya, Xb, Yb, XaXb, YaYb, XaYb, YaXb;
	int		N;
};

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */

static XArray			*gX;
static vector<RgdSums>	vS;
static int				nthr;






/* --------------------------------------------------------------- */
/* _RgdSums ------------------------------------------------------ */
/* --------------------------------------------------------------- */

void* _RgdSums( void* ithr )
{
	int	nb = vS.size();

// For each B-layer...

	for( int ib = (long)ithr; ib < nb; ib += nthr ) {

		// Form the sums

		int						ia	= ib + 1;
		const Rgns&				Rb	= vR[ib];
		const Rgns&				Ra	= vR[ia];
		const vector<double>&	xb	= gX->X[ib];
		const vector<double>&	xa	= gX->X[ia];
		RgdSums&				S	= vS[ib];

		// For each A-layer rgn...

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( Ra.flag[ir] )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const TAffine*		Ta = &X_AS_AFF( xa, ir );
			const TAffine*		Tb;
			int					lastb	= -1,
								np		= P.size();

			// For each of its points...

			for( int ip = 0; ip < np; ++ip ) {

				const CorrPnt&	C = vC[P[ip]];

				// Want only zA onto zB

				if( C.z1 != ia || C.z2 != ib )
					continue;

				if( C.i2 != lastb ) {

					if( Rb.flag[C.i2] )
						continue;

					Tb = &X_AS_AFF( xb, C.i2 );
					lastb = C.i2;
				}

				Point	pa = C.p1,
						pb = C.p2;

				Ta->Transform( pa );
				Tb->Transform( pb );

				S.Xa += pa.x;
				S.Ya += pa.y;
				S.Xb += pb.x;
				S.Yb += pb.y;

				S.XaXb += pa.x * pb.x;
				S.YaYb += pa.y * pb.y;
				S.XaYb += pa.x * pb.y;
				S.YaXb += pa.y * pb.x;

				++S.N;
			}
		}

		// Set transform

		if( S.N < 2 )
			continue;

		double
		theta = atan(
			(S.YaXb - S.XaYb + (S.Xa*S.Yb - S.Ya*S.Xb)/S.N) /
			((S.Xa*S.Xb + S.Ya*S.Yb)/S.N - S.XaXb - S.YaYb)
		),
		c  = cos( theta ),
		s  = sin( theta ),
		kx = (S.Xb - c*S.Xa + s*S.Ya) / S.N,
		ky = (S.Yb - s*S.Xa - c*S.Ya) / S.N;

		S.A.t[0] = c; S.A.t[1] = -s; S.A.t[2] = kx;
		S.A.t[3] = s; S.A.t[4] =  c; S.A.t[5] = ky;
	}

	return NULL;
}

/* --------------------------------------------------------------- */
/* CalcMyPairwiseTForms ------------------------------------------ */
/* --------------------------------------------------------------- */

// Pairwise layer sums and tforms can be done in parallel.
//
static void CalcMyPairwiseTForms( XArray &X )
{
	gX = &X;

	int	nb = vR.size() - 1;	// this many b-layers

	vS.resize( nb );

	nthr = maxthreads;

	if( nthr > nb )
		nthr = nb;

	if( !EZThreads( _RgdSums, nthr, 1, "_RgdSums" ) )
		exit( 42 );
}

/* --------------------------------------------------------------- */
/* WriteMyTForms ------------------------------------------------- */
/* --------------------------------------------------------------- */

// The adjustment tform needed for a given layer is the cummulative
// product of all adjustments to layers below. At this point, all
// worker nodes except the last must write their pairwise results
// to files "Untwist/id.txt".
//
// Note that printf with %.21f writes doubles to full precision
// such that scanf with %lf recovers the identical value.
//
// We bother only because of long product chains. With 21 digits
// there is no difference between many small subblocks and doing
// all together. With standard precision (%f = %.6f), differences
// appear in the 5th-most significant digit over 200 layers.
//
static void WriteMyTForms()
{
// If I'm not the last worker...

	if( wkid >= nwks - 1 )
		return;

	DskCreateDir( "Untwist", stdout );

	char	buf[32];
	sprintf( buf, "Untwist/%d.txt", wkid );
	FILE	*f = FileOpenOrDie( buf, "w" );

// Each entry affects layer zA or higher, so write lines:
// 'zA A'

	int	nb = vS.size();

	for( int ib = 0; ib < nb; ++ib ) {

		const TAffine&	A = vS[ib].A;

		fprintf( f, "%d %.21f %.21f %.21f %.21f %.21f %.21f\n",
		vR[ib+1].z,
		A.t[0], A.t[1], A.t[2], A.t[3], A.t[4], A.t[5] );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* AccumulateBefores --------------------------------------------- */
/* --------------------------------------------------------------- */

// Form product of all affines with z <= my stating z.
//
// We must be careful here because workers have overlapping
// zo ranges, while our product must include each z in order
// and once only. This is fixed just by requiring monotonic
// increase in z.
//
static void AccumulateBefores( TAffine &A0 )
{
	TAffine	A;
	int		z0		= vR[zolo].z,
			zlast	= -1,
			z;

	for( int iw = 0; iw < wkid; ++iw ) {

		char	buf[32];
		sprintf( buf, "Untwist/%d.txt", iw );
		FILE	*f = FileOpenOrDie( buf, "r" );

		while( 7 == fscanf( f, "%d%lf%lf%lf%lf%lf%lf\n",
			&z,
			&A.t[0], &A.t[1], &A.t[2],
			&A.t[3], &A.t[4], &A.t[5] ) ) {

			if( z > z0 ) {
				fclose( f );
				return;
			}

			// monotonic z rule

			if( z <= zlast )
				continue;

			A0		= A * A0;
			zlast	= z;
		}

		fclose( f );
	}
}

/* --------------------------------------------------------------- */
/* Apply --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Apply A0 and pairwise affines.
//
static void Apply( TAffine &A0 )
{
	int	nz = vR.size();

	for( int iz = (wkid ? 0 : 1); iz < nz; ++iz ) {

		if( iz )
			A0 = vS[iz - 1].A * A0;

		const Rgns&		R = vR[iz];
		vector<double>&	x = gX->X[iz];

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( R.flag[ir] )
				continue;

			TAffine& T = X_AS_AFF( x, ir );

			T = A0 * T;
		}
	}
}

/* --------------------------------------------------------------- */
/* UntwistAffines ------------------------------------------------ */
/* --------------------------------------------------------------- */

// The standard way to obtain an external scaffold for use with
// the -prior option is to use XMLGetTF on HiRes.xml; the result
// of aligning low-res strips from very good montages. However,
// better angle calculations can be made later, after the down
// correspondence points are formed. Adjustments are calculated
// here as rigid transforms T{theta, kx, ky) from layer a to b.
//
// Let c=cos(theta), s=sin(theta), Sum over all point-pairs:
//
// E = Sum[Xb - cXa + sYa - kx]^2 + [Yb - sXa - cYa - ky]^2
//
// The params u = {theta,kx,ky) are determined by dE/du = 0.
// If we use notation [arg] => Sum[argi] over point-pairs,
//
// kx = ([Xb] - c[Xa] + s[Ya]) / N
// ky = ([Yb] - s[Xa] - c[Ya]) / N
//
//				[YaXb] - [XaYb] + ([Xa][Yb] - [Ya][Xb]) / N
// tan(theta) = --------------------------------------------
//				([Xa][Xb] + [Ya][Yb]) / N - [XaXb] - [YaYb]
//
void UntwistAffines( XArray &X )
{
	if( X.NE != 6 || zolo >= zohi )
		return;

	clock_t	t0 = StartTiming();

	CalcMyPairwiseTForms( X );
	WriteMyTForms();

	MPIWaitForOthers();

	TAffine	A0;
	AccumulateBefores( A0 );
	Apply( A0 );
	vS.clear();

	StopTiming( stdout, "Untwist", t0 );
}


