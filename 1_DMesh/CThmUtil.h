

#pragma once


#include	"FoldMask.h"
#include	"CPixPair.h"
#include	"CThmScan.h"


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef struct {
	vector<double>	v;
	vector<Point>	p;
	Point			O;		// subimage ref. point
	int				w, h;	// sub-img dims
} SubI;

typedef struct {
	SubI		a, b;
} OlapRec;

/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CThmUtil {

private:
	const PicSpecs	&A, &B;
	const int		acr, bcr;
	const PixPair	&px;
	long			&OLAP2D;
	FILE			*flog;

public:
	TAffine	&Tab;
	double	ang0,
			HFANGDN,
			HFANGPR,
			RTRSH;
	int		OLAP1D,
			MODE,
			LIMXY;

public:
	CThmUtil(
		const PicSpecs	&A,
		int				acr,
		const PicSpecs	&B,
		int				bcr,
		const PixPair	&px,
		TAffine			&Tab,
		long			&OLAP2D,
		FILE			*flog )
		:
			A(A), B(B), acr(acr), bcr(bcr),
			px(px), Tab(Tab), OLAP2D(OLAP2D),
			flog(flog), ang0(0.0)
		{};

	void SetParams(
		double	HFANGDN,
		double	HFANGPR,
		double	RTRSH,
		int		OLAP1D,
		int		MODE,
		int		LIMXY )
		{
			this->HFANGDN=HFANGDN; this->HFANGPR=HFANGPR;
			this->RTRSH=RTRSH; this->OLAP1D=OLAP1D;
			this->MODE=MODE; this->LIMXY=LIMXY;
		};

	bool Echo( vector<TAffine> &guesses );
	bool FromLog( vector<TAffine> &guesses );

	int SetStartingAngle( const TAffine &Tdfm, double CTR );

	void SubI_ThesePoints(
		SubI					&S,
		const vector<double>	&v,
		const vector<Point>		&p );

	bool SubI_ThisBox(
		SubI					&S,
		const vector<double>	&v,
		const vector<Point>		&p,
		const IBox				&Bolap );

	bool Olap_WholeImage(
		OlapRec				&olp,
		const vector<Point>	&apts,
		const vector<Point>	&bpts );

	bool Olap_TheseBoxes_NoCR(
		OlapRec		&olp,
		const IBox	&Ba,
		const IBox	&Bb );

	bool Olap_WholeImage_NoCR( OlapRec &olp );

	void GetOlapBoxes( IBox &Ba, IBox &Bb, double XYCONF );

	bool Crop(
		OlapRec				&olp,
		const ConnRegion	&acr,
		const ConnRegion	&bcr,
		double				XYCONF );

	bool Crop_NoCR( OlapRec &olp, double XYCONF );

	bool MakeThumbs(
		ThmRec			&thm,
		const OlapRec	&olp,
		int				decfactor );

	bool Disc(
		CorRec			&best,
		CThmScan		&S,
		ThmRec			&thm,
		const OlapRec	&olp,
		int				PRETWEAK );

	void DebugSweepKill( CThmScan &S, ThmRec thm );

	bool Sweep(
		CorRec		&best,
		CThmScan	&S,
		ThmRec		&thm,
		int			nPrior );

	void IsectToImageCoords(
		CorRec		&best,
		Point		aO,
		const Point	&bO );

	void RecordSumSqDif( const TAffine &T );

	void FullScaleReportToLog( CorRec &best );

	double XYChange(
		CorRec			CR,
		ThmRec			&thm,
		const OlapRec	&olp );

	bool Check_LIMXY( const TAffine &Tbest );

	void TabulateResult( const CorRec &best, int err );
	bool Failure( const CorRec &best, int err );

	bool Finish(
		CorRec			&best,
		CThmScan		&S,
		ThmRec			&thm,
		const OlapRec	&olp,
		int				TWEAKS );
};


