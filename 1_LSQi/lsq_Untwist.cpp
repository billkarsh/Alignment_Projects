

#include	"lsq_Globals.h"
#include	"lsq_Msg.h"
#include	"lsq_Untwist.h"

#include	"EZThreads.h"
#include	"Disk.h"
#include	"File.h"
#include	"TAffine.h"
#include	"Timer.h"

#include	<stdlib.h>


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
/* CommandSharing ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void CommandSharing()
{
	if( !wkid && nwks > 1 )
		MsgSend( "untwist-share" );
}

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

int npt = 0;

		for( int ir = 0; ir < Ra.nr; ++ir ) {

			if( !Ra.used[ir] )
				continue;

			const vector<int>&	P  = Ra.pts[ir];
			const TAffine*		Ta = &X_AS_AFF( xa, ir );
			const TAffine*		Tb;
			int					lastb	= -1;
			int					np		= P.size();

			// For each of its points...

			for( int ip = 0; ip < np; ++ip ) {

				const CorrPnt&	C = vC[P[ip]];

				// Want only zA onto zB

				if( C.z1 != ia || C.z2 != ib )
					continue;

				if( C.i2 != lastb ) {

					if( !Rb.used[C.i2] )
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

	nthr = 16;

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

		fprintf( f, "%d %f %f %f %f %f %f\n",
		vR[ib+1].z,
		A.t[0], A.t[1], A.t[2], A.t[3], A.t[4], A.t[5] );
	}

	fclose( f );
}

/* --------------------------------------------------------------- */
/* AckMineShared ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void AckMineShared()
{
	if( wkid > 0 ) {

		if( !MsgWaitMatch( "untwist-share" ) )
			exit( 32 );

		MsgAck( wkid );
	}
}

/* --------------------------------------------------------------- */
/* WaitAllShared ------------------------------------------------- */
/* --------------------------------------------------------------- */

// Zero syncs by waiting for all acks, then issues "untwist".
// All others wait for "untwist".
//
static void WaitAllShared()
{
	if( nwks <= 1 )
		return;

	if( !wkid ) {

		if( !MsgWaitAck( nwks ) )
			exit( 32 );

		MsgSend( "untwist" );
	}
	else {
		if( !MsgWaitMatch( "untwist" ) )
			exit( 32 );
	}
}

/* --------------------------------------------------------------- */
/* AccumulateBefores --------------------------------------------- */
/* --------------------------------------------------------------- */

// Form product of all affines with z <= my stating z.
//
static void AccumulateBefores( TAffine &A0 )
{
	TAffine	A;
	int		z0 = vR[zolo].z,
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

			A0 = A * A0;
		}

		fclose( f );
	}
}

/* --------------------------------------------------------------- */
/* Untwist ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Apply A0 and pairwise affines.
//
static void Untwist( TAffine &A0 )
{
	int	nz = vR.size();

	for( int iz = (wkid ? 0 : 1); iz < nz; ++iz ) {

		if( iz )
			A0 = vS[iz - 1].A * A0;

		const Rgns&		R = vR[iz];
		vector<double>&	x = gX->X[iz];

		for( int ir = 0; ir < R.nr; ++ir ) {

			if( !R.used[ir] )
				continue;

			TAffine& T = X_AS_AFF( x, ir );

			T = A0 * T;
		}
	}
}

/* --------------------------------------------------------------- */
/* Finish -------------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Finish()
{
	if( nwks <= 1 )
		return;

	if( !wkid ) {

		if( !MsgWaitAck( nwks ) )
			exit( 32 );
	}
	else
		MsgAck( wkid );
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

	CommandSharing();

	CalcMyPairwiseTForms( X );

	WriteMyTForms();

	AckMineShared();

	WaitAllShared();

	TAffine	A0;

	AccumulateBefores( A0 );

	Untwist( A0 );

// Done

	vS.clear();
	Finish();

	StopTiming( stdout, "Untwist", t0 );
}


