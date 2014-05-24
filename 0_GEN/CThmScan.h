

#pragma once


#include	"GenDefs.h"
#include	"TAffine.h"

#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

enum thmerrs {
	errOK			= 0,
	errLowRDenov	= 1,
	errLowRPrior	= 2
};

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	vector<double>	av, bv;
	vector<Point>	ap, bp;
	vector<CD>		ftc;		// fourier transform cache
	long			reqArea;
	int				olap1D,
					scl;		// for caller convenience
} ThmRec;

class CorRec {
public:
	TAffine			T;
	double			X, Y,
					A, R;
public:
	CorRec()	{};
	CorRec( double A ) : A(A) {};
};

typedef void (*T_NewAngProc)(
	int		&Ox,
	int		&Oy,
	int		&Rx,
	int		&Ry,
	double	deg );

/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CThmScan {

	friend void* _TCDGet( void* icstart );
	friend bool BigEnough( int sx, int sy, void *a );

private:
	class PTWRec {
	public:
		double	a;
		int		sel;
	public:
		PTWRec( int sel, double a ) : a(a), sel(sel) {};
	};

	class ThrCorDat {
	public:
		// used only by TCDGet
		int				nthr;
	public:
		// caller I/O
		ThmRec			*thm;
		vector<CorRec>	vC;
		vector<PTWRec>	vT;
	};

private:
	ThrCorDat		TCD;
	T_NewAngProc	newAngProc;
	FILE			*flog;
	TAffine			Tdfm,
					Tptwk;
	double			rthresh,
					nbmaxht;
	int				err,
					swpConstXY,
					swpPretweak,
					swpNThreads,
					useCorrR,
					Ox, Oy, Rx, Ry,
					olap1D;

private:
	void _TCDDo1( int ic );
	void TCDGet( int nthr );

	double PTWInterp(
		double	&ynew,
		int		sel,
		double	x1,
		double	y1,
		double	d,
		double	deg,
		ThmRec	&thm );

public:
	CThmScan();

	void Initialize( FILE* flog, CorRec &best );

	void SetFlog( FILE* flog )
		{this->flog = flog;};

	void SetTdfm( const TAffine &Tdfm )
		{this->Tdfm = Tdfm;};

	void SetRThresh( double rthresh )
		{this->rthresh = rthresh;};

	void SetNbMaxHt( double nbmaxht )
		{this->nbmaxht = nbmaxht;};

	void SetSweepConstXY( int swpConstXY )
		{this->swpConstXY = swpConstXY;};

	void SetSweepPretweak( int swpPretweak )
		{this->swpPretweak = swpPretweak;};

	void SetSweepNThreads( int swpNThreads )
		{this->swpNThreads = swpNThreads;};

	void SetUseCorrR( int useCorrR )
		{this->useCorrR = useCorrR;};

	void SetDisc( int Ox, int Oy, int Rx, int Ry )
		{this->Ox = Ox; this->Oy = Oy; this->Rx = Rx; this->Ry = Ry;};

	void SetNewAngProc( T_NewAngProc proc )
		{newAngProc = proc;};

	void GetTptwk( TAffine &Tptwk )
		{Tptwk = this->Tptwk;};

	void SetTptwk( const TAffine &Tptwk )
		{this->Tptwk = Tptwk;};

	int GetErr()
		{return err;};

	bool Pretweaks( double bestR, double deg, ThmRec &thm );

	void RFromAngle(
		CorRec	&C,
		double	deg,
		ThmRec	&thm,
		int		iaux = -1 );

	void DebugAngs(
		int		alyr,
		int		atil,
		int		blyr,
		int		btil,
		double	center,
		double	hlfwid,
		double	step,
		ThmRec	&thm );

	double AngleScanMaxR(
		CorRec	&best,
		double	center,
		double	hlfwid,
		double	step,
		bool	favorCenter,
		ThmRec	&thm );

	double AngleScanConstXY(
		CorRec	&best,
		double	center,
		double	hlfwid,
		double	step,
		ThmRec	&thm );

	double AngleScanSel(
		CorRec	&best,
		double	center,
		double	hlfwid,
		double	step,
		bool	favorCenter,
		ThmRec	&thm );

	double AngleScanWithTweaks(
		CorRec	&best,
		double	center,
		double	hlfwid,
		double	step,
		ThmRec	&thm );

	double PeakHunt( CorRec &best, double hlfwid, ThmRec &thm );

	bool UsePriorAngles(
		CorRec	&best,
		int		nprior,
		double	ang0,
		double	hfangpr,
		ThmRec	&thm );

	bool DenovoBestAngle(
		CorRec	&best,
		double	ang0,
		double	hfangdn,
		double	step,
		ThmRec	&thm,
		bool	failmsg );

	void PostTweaks( CorRec &best, ThmRec &thm );

	void FinishAtFullRes( CorRec &best, ThmRec &thm );
};


