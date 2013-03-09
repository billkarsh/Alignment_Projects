

#pragma once


#include	"lsq_Types.h"

#include	"LinEqu.h"


/* --------------------------------------------------------------- */
/* MDL base classs ----------------------------------------------- */
/* --------------------------------------------------------------- */

class MDL {

protected:
	const int	NT, NX;

protected:
	MDL( int NT, int NX ) : NT(NT), NX(NX) {};
	void PrintMagnitude( const vector<double> &X );

	void NewAffine(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength,
		double			square_strength,
		int				nTr,
		int				itr );

	void NewHmgphy(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength,
		double			square_strength,
		int				nTr,
		int				itr );

private:
	virtual void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength ) = 0;

	virtual void SetIdentityTForm(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				itr ) = 0;

	virtual void SetUniteLayer(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		int				unite_layer,
		const char		*tfm_file ) = 0;

	virtual void SolveWithSquareness(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTr,
		double			square_strength ) = 0;

	virtual void SolveWithUnitMag(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTR,
		double			scale_strength ) = 0;

	virtual void RescaleAll(
		vector<double>	&X,
		double			sc ) = 0;

	virtual void RotateAll(
		vector<double>	&X,
		double			degcw ) = 0;

	virtual void NewOriginAll(
		vector<double>	&X,
		double			xorg,
		double			yorg ) = 0;

public:
	int MinLinks()	{return NX/2;};

	void SolveSystem(
		vector<double>	&X,
		int				nTr,
		int				gW,
		int				gH,
		double			same_strength,
		double			square_strength,
		double			scale_strength,
		int				unite_layer,
		const char		*tfm_file );

	void Bounds(
		double					&xbnd,
		double					&ybnd,
		vector<double>			&X,
		int						gW,
		int						gH,
		const vector<double>	&lrbt,
		double					degcw,
		FILE					*FOUT );

	virtual void WriteTransforms(
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						bstrings,
		FILE					*FOUT ) = 0;

	virtual void WriteTrakEM(
		double					xmax,
		double					ymax,
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						gW,
		int						gH,
		double					trim,
		int						xml_type,
		int						xml_min,
		int						xml_max ) = 0;

	virtual void L2GPoint(
		Point			&p,
		vector<double>	&X,
		int				itr ) = 0;

	virtual void L2GPoint(
		vector<Point>	&p,
		vector<double>	&X,
		int				itr ) = 0;
};


