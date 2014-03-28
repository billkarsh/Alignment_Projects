

#pragma once


#include	"PipeFiles.h"
#include	"CCropMask.h"


/* --------------------------------------------------------------- */
/* class ConnRegion ---------------------------------------------- */
/* --------------------------------------------------------------- */

class ConnRegion {

// Note: All coordinates here are scaled.

public:
	vector<Point>	pts;	// pixels within the region
	double			dx, dy;	// deltas to line up with image2
	IBox			B;		// bounding box in original image
	int				id;		// this region's mask value

public:
	ConnRegion()	{B.L = B.B = BIG; B.R = B.T = -BIG;};
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */

uint8* GetFoldMask(
	const string		&idb,
	const PicSpecs		&P,
	const char			*forcepath,
	const vector<uint8>	&resmsk,
	CCropMask			*CM,
	int					wf,
	int					hf,
	bool				nofile,
	bool				transpose,
	bool				force1rgn );

void SetWithinSectionBorders( uint8* foldMask, int wf, int hf );

void SetBoundsAndColors(
	vector<ConnRegion>	&cr,
	uint8*				foldMask,
	int					wf,
	int					hf );

void SetBoundsAndColors(
	vector<ConnRegion>		&cr,
	uint8*					foldMask,
	const vector<double>	&valid,
	int						wf,
	int						hf,
	double					thresh,
	int						D );

void ConnRgnsFromFoldMask(
	vector<ConnRegion>	&cr,
	const uint8*		foldMask,
	int					wf,
	int					hf,
	int					scale,
	uint32				minpts,
	FILE				*flog );

void ConnRgnForce1( vector<ConnRegion> &cr, int ws, int hs );


