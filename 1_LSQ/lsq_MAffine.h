

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
		double			sc );

	void SetIdentityTForm(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				itr );

	void SetUniteLayer(
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		double			sc );

	void SolveWithSquareness(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTr );

	void SolveWithUnitMag(
		vector<double>	&X,
		vector<LHSCol>	&LHS,
		vector<double>	&RHS,
		int				nTR );

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

	void DeviantAffines(
		const vector<double>	&T,
		const vector<double>	&X );

	void AffineFromTransWt( vector<double> &X, int nTr );

	void AffineFromAffine2( vector<double> &X, int nTr );

	void AffineFromAffine( vector<double> &X, int nTr );

	void AffineFromTrans( vector<double> &X, int nTr );

	void SolveSystemStandard( vector<double> &X, int nTr );

public:
	void SolveSystem( vector<double> &X, int nTr );

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
		double					trim,
		int						xml_type,
		int						xml_min,
		int						xml_max );

	void WriteJython(
		const vector<zsort>		&zs,
		const vector<double>	&X,
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


