

#pragma once


#include	"BK_DEFS_GEN.h"

#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Mask ---------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CMask {

private:
	enum MaskConsts {
		kMaxHPix		= 2048,
		kMaxVPix		= 2048,
		kMaxPix			= kMaxHPix * kMaxVPix,
		kMapWords		= kMaxPix / WordBits,
		kBucketBytes	= kMaxPix * 4
	};

	enum StdDefConsts {
	// global thresholding methods
		stdDefMthdFixed			= 0,
		stdDefMthdAddFixed		= 1,
		stdDefMthdAddRel		= 2,
		stdDefMthdMulFixed		= 3,

	// local rethreholding methods
		stdLRTMthdNone			= 0,
		stdLRTMthdPercentile	= 1,
		stdLRTMthdPercentPk		= 2
	};

private:
	FILE	*flog;

// StdDef params
	UInt8	threshMethod,	// enum threshold method
			lrtMethod,		// enum loc rethresh method
			fillHoles,		// checkBox: true or false
			reserved;
	SInt32	erodePix;		// pixels
	FP32	threshConst,	// range: method-dependent
			rethreshFrac;	// range:[0,1] inclusive

// StdFlt params
	UInt8	checkBBox,		// BBox flt applied?
			check_I,		// AveI flt applied?
			reserved1,
			reserved2;
	SInt32	minBBox,		// pixels
			maxBBox,		// pixels
			min_I,			// counts
			max_I;			// counts

// image
	const UInt16	*m_mrkPix;
	UInt32			m_hPix,
					m_vPix;
	U16Box			m_roi;

// working...
	FP32			m_peak,
					m_thresh;

// bitmaps
	UInt8	m_scratch[kBucketBytes];
	UInt32	m_iniMap[kMapWords],
			m_fltMap[kMapWords],
			m_scnMap[kMapWords],
			m_trkMap[kMapWords],
			m_isoMap[kMapWords],
			m_tmpMap0[kMapWords],
			m_tmpMap1[kMapWords];

private:
	void   ThreshFromPeak();
	void   GetMarkerPeak();
	UInt32 GetThresh();

	void   LocalRethreshold(
		UInt32	*dstMap,
		UInt32	*trackingMap,
		UInt32	*scanningMap,
		UInt32	*tmpMap0,
		UInt32	*tmpMap1,
		UInt32	T );

	UInt32 InitObjectDetectionMaps(
		UInt32	*auxObjectMap,
		UInt32	*trackingMap,
		UInt32	*scanningMap,
		UInt32	*tmpMap0,
		UInt32	*tmpMap1 );

	int    BBoxPass( const U16BoxPtr rBlob );
	int    Marker_I_Pass( const UInt32 *map, const U16BoxPtr rBlob );
	void   ApplyFilters();

public:
	void ReadParamFile( FILE *flog, const char *path );

	void MaskFromImage(
		const char		*outpath,
		const UInt16	*src,
		int				hPix,
		int				vPix );
};


