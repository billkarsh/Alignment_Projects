

#pragma once


#include	"lsq_MDL.h"


/* --------------------------------------------------------------- */
/* Affine MDL ---------------------------------------------------- */
/* --------------------------------------------------------------- */

class MAffine : public MDL {

public:
	MAffine() : MDL( 6, 6 ) {};

private:
	void SetPointPairs(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		double			same_strength );

	void SetIdentityTForm(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				itr );

	void SetUniteLayer(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc,
		int				unite_layer,
		const char		*tfm_file );

	void SolveWithSquareness(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTr,
		double			square_strength );

	void SolveWithUnitMag(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTR,
		double			scale_strength );

	void RescaleAll(
		vector<double>	&X,
		double			sc );

	void RotateAll(
		vector<double>	&X,
		double			degcw );

	void NewOriginAll(
		vector<double>	&X,
		double			xorg,
		double			yorg );

	void AffineEquTransWt(
		vector<double>	&X,
		int				nTr,
		int				gW,
		int				gH,
		double			same_strength,
		double			square_strength );

	void AffineEquTrans(
		vector<double>	&X,
		int				nTr,
		int				gW,
		int				gH,
		double			square_strength );

	void SolveSystemStandard(
		vector<double>	&X,
		int				nTr,
		int				gW,
		int				gH,
		double			same_strength,
		double			square_strength,
		double			scale_strength,
		int				unite_layer,
		const char		*tfm_file );

public:
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

	void WriteTransforms(
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						bstrings,
		FILE					*FOUT );

	void WriteTrakEM(
		double					xmax,
		double					ymax,
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						gW,
		int						gH,
		double					trim,
		int						xml_type,
		int						xml_min,
		int						xml_max );

	void WriteJython(
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						gW,
		int						gH,
		double					trim,
		int						Ntr );

	void G2LPoint(
		Point					&p,
		const vector<double>	&X,
		int						itr );

	void L2GPoint(
		Point					&p,
		const vector<double>	&X,
		int						itr );

	void L2GPoint(
		vector<Point>			&p,
		const vector<double>	&X,
		int						itr );
};


