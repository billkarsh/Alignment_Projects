

#pragma once


#include	"lsq_Types.h"

#include	"TAffine.h"
#include	"LinEqu.h"


/* --------------------------------------------------------------- */
/* MDL base classs ----------------------------------------------- */
/* --------------------------------------------------------------- */

class MDL {

protected:
	const int	NT, NX;

protected:
	MDL( int NT, int NX ) : NT(NT), NX(NX) {};

	double Magnitude( const vector<double> &X, int itr );

	void PrintMagnitude( const vector<double> &X );

private:
	virtual void RotateAll(
		vector<double>	&X,
		double			degcw ) = 0;

	virtual void NewOriginAll(
		vector<double>	&X,
		double			xorg,
		double			yorg ) = 0;

public:
	int MinLinks()	{return NX/2;};

	virtual void SolveSystem(
		vector<double>	&X,
		int				nTr,
		int				gW,
		int				gH,
		double			same_strength,
		double			square_strength,
		double			scale_strength,
		int				unite_layer,
		const char		*tfm_file ) = 0;

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

	virtual void WriteJython(
		const vector<zsort>		&zs,
		const vector<double>	&X,
		int						gW,
		int						gH,
		double					trim,
		int						Ntr ) = 0;

	virtual void G2LPoint(
		Point					&p,
		const vector<double>	&X,
		int						itr ) = 0;

	virtual void L2GPoint(
		Point					&p,
		const vector<double>	&X,
		int						itr ) = 0;

	virtual void L2GPoint(
		vector<Point>			&p,
		const vector<double>	&X,
		int						itr ) = 0;

	virtual TAffine EqvAffine(
		const vector<double>	&X,
		int						itr );
};


