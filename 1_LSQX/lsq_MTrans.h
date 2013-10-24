

#pragma once


#include	"lsq_MDL.h"


/* --------------------------------------------------------------- */
/* Trans MDL ----------------------------------------------------- */
/* --------------------------------------------------------------- */

class MTrans : public MDL {

public:
	MTrans() : MDL( 6, 2 ) {};

private:
	void RotateAll(
		vector<double>	&X,
		double			degcw );

	void NewOriginAll(
		vector<double>	&X,
		double			xorg,
		double			yorg );

public:
	void SolveSystem( vector<double> &X, int nTr );

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


