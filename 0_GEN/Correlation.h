

#pragma once


#include	"GenDefs.h"
#include	"CPoint.h"

#include	<stdio.h>

#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* FFT support --------------------------------------------------- */
/* --------------------------------------------------------------- */

int FFTSize( int n1, int n2 );

int FFTSizeSQR( int w1, int h1, int w2, int h2 );

int FFT_2D(
	vector<CD>				&out,
	const vector<double>	&in,
	int						Nfast,
	int						Nslow );

void IFT_2D(
	vector<double>			&out,
	const vector<CD>		&in,
	int						Nfast,
	int						Nslow );

/* --------------------------------------------------------------- */
/* Correlation support ------------------------------------------- */
/* --------------------------------------------------------------- */

void ParabPeakFFT(
	double			&xpk,
	double			&ypk,
	int				d,
	const double	*I,
	int				nx,
	int				ny );

void PrintCorLandscape(
	double			biggest,
	int				bigx,
	int				bigy,
	int				Ox,
	int				Oy,
	int				radius,
	int				lim,
	int				step,
	const double	*I,
	int				nx,
	int				ny,
	double			norm,
	FILE			*flog );

/* --------------------------------------------------------------- */
/* RVectors ------------------------------------------------------ */
/* --------------------------------------------------------------- */

double RVectors(
	const vector<double>	&a,
	const vector<double>	&b );

/* --------------------------------------------------------------- */
/* CorrPatches --------------------------------------------------- */
/* --------------------------------------------------------------- */

typedef bool (*EvalType)( int sx, int sy, void *v );

double CorrPatches(
	FILE					*flog,
	int						verbose,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						Ox,
	int						Oy,
	int						radius,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 );

double CorrPatchToImage(
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<double>	&i2,
	int						Ox,
	int						Oy,
	int						radius,
	bool					bFilter );

/* --------------------------------------------------------------- */
/* Mesh Optimization --------------------------------------------- */
/* --------------------------------------------------------------- */

void GradDescStep(
	vector<double>					&bnew,
	vector<Point>					&dRdc,
	const vector<Point>				&ac,
	const vector<vector<double> >	&am,
	const vector<double>			&av,
	const vector<double>			&bimg,
	int								w,
	int								h );

double CorrVectors(
	FILE					*flog,
	const vector<double>	&a,
	const vector<double>	&b,
	int						*nnz = NULL );

double ImproveControlPts(
	vector<Point>					&ac,
	const vector<vector<double> >	&am,
	const vector<double>			&av,
	const vector<double>			&bimg,
	int								w,
	int								h,
	FILE							*flog,
	const char						*describe,
	double							iniThresh,
	double							finThresh );

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

double CorrPatchesRQ(
	FILE					*flog,
	int						verbose,
	double					&Q,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 );

double CorrPatchesMaxQ(
	FILE					*flog,
	int						verbose,
	double					&Q,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	vector<CD>				&fft2 );

double CorrPatchesMaxR(
	FILE					*flog,
	int						verbose,
	double					&dx,
	double					&dy,
	const vector<Point>		&ip1,
	const vector<double>	&iv1,
	const vector<Point>		&ip2,
	const vector<double>	&iv2,
	int						nnegx,
	int						nposx,
	int						nnegy,
	int						nposy,
	EvalType				LegalRgn,
	void*					arglr,
	EvalType				LegalCnt,
	void*					arglc,
	double					qtol,
	vector<CD>				&fft2 );


