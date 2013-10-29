

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

	void DevFromTrans(
		const vector<double>	&T,
		const vector<double>	&X );

	void DevFromPrior(
		const vector<double>	&A,
		const vector<double>	&X );

	void LoadAffTable(
		vector<double>	&X,
		int				&z0,
		int				&nz,
		int				nTr );

	void UntwistScaffold(
		vector<double>	&X,
		int				z0,
		int				nz );

	void AffineFromFile( vector<double> &X, int nTr );

	void AffineFromFile2( vector<double> &X, int nTr );

	void AffineFromTransWt( vector<double> &X, int nTr );

	// ----------------------------------------------------------
	// Iterative experiments
	void Fill_myc( const vector<double> &X );
	void GetStageT( vector<double> &X, int nTr );
	void GetTableT( vector<double> &X, int nTr );
	void GetScaffT( vector<double> &X, int nTr );
	void OnePass(
		vector<double>	&Xout,
		vector<double>	&Xin,
		vector<double>	&S,
		int				nTr,
		double			w );
	// ----------------------------------------------------------

	void SolveSystemStandard( vector<double> &X, int nTr );

public:
	void SolveSystem( vector<double> &X, int nTr );

	void UpdateScaffold( vector<double> &X, int nTr );

	void WriteTransforms(
		const vector<double>	&X,
		int						bstrings,
		FILE					*FOUT );

	void WriteTrakEM(
		double					xmax,
		double					ymax,
		const vector<double>	&X,
		double					trim,
		int						xml_type,
		int						xml_min,
		int						xml_max );

	void WriteJython(
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


