/* BK_HST ------------------------------------------------------------
 *
 * Histogramming.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_HST.h"
#include	"BK_MEM.h"






/* HSTAltImageRowsUInt16 ---------------------------------------------
 *
 * Histogram every other row in 16-bit image data using fixed
 * bin width of one, i.e., values must be in range [0, nBins-1].
 * Values greater than or equal to nBins are collected in the
 * oFlow parameters (if non-NULL).
 *
 * Return pixel count, excluding overflow.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * hPix must be an even number.
 *
 * -----
 * Notes
 * -----
 *
 * 16-bit pixel pairs are read into a 32-bit register and extracted
 * separately. On X86 family processors, the 16-bit number at the
 * lower RAM address occupies the lower 16 bits of the register.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 HSTAltImageRowsUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	int				hPix,
	int				vPix )
{
	UInt32	*pPair;
	UInt32	pair, v0, v1, oflowC, oflowS;
	int		h, v, rowBytes;

/* ----------------- */
/* Initialize output */
/* ----------------- */

	MEMZeroBytes( binArray, nBins * sizeof(UInt32) );
	oflowC	= 0;
	oflowS	= 0;

/* ---------- */
/* Accumulate */
/* ---------- */

	rowBytes	= hPix * sizeof(UInt16);
	hPix		>>= 1;
	pPair		= (UInt32*)data16Bit;

	for( v = 0; v < vPix;
		v += 2,
		pPair = (UInt32*)((char*)pPair + rowBytes) ) {

		for( h = 0; h < hPix; ++h ) {

			pair	= *pPair++;
			v0		= pair & LoMask;
			v1		= pair >> HalfBits;

			if( v0 < nBins )
				++binArray[v0];
			else {
				oflowS += v0;
				++oflowC;
			}

			if( v1 < nBins )
				++binArray[v1];
			else {
				oflowS += v1;
				++oflowC;
			}
		}
	}

	if( oFlowCnt )
		*oFlowCnt = oflowC;

	if( oFlowSum )
		*oFlowSum = oflowS;

	return (vPix + (vPix & 1)) / 2 * hPix - oflowC;
}


/* HSTImageBoxUInt16 -------------------------------------------------
 *
 * Histogram all boxed pixels in 16-bit image data using fixed
 * bin width of one, i.e., values must be in range [0, nBins-1].
 * Values greater than or equal to nBins are collected in the
 * oFlow parameters (if non-NULL).
 *
 * Return pixel count, excluding overflow.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 HSTImageBoxUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	int				hPix,
	int				vPix,
	const U16BoxPtr	box )
{
	UInt32	imgRowBytes, N, x, oflowC, oflowS;
	int		v, h, top, left, bot, right;

/* --------------- */
/* Initializations */
/* --------------- */

	MEMZeroBytes( binArray, nBins * sizeof(UInt32) );
	N			= 0;
	oflowC		= 0;
	oflowS		= 0;

	imgRowBytes	= hPix * sizeof(UInt16);

/* box limits */

	if( box ) {

		top		= box->top;
		left	= box->left;
		bot		= box->bottom;
		right	= box->right;

		if( bot > (int)vPix )
			bot = vPix;

		if( right > (int)hPix )
			right = hPix;

		if( top >= bot )
			goto exit;

		if( left >= right )
			goto exit;
	}
	else {
		top		= 0;
		left	= 0;
		bot		= vPix;
		right	= hPix;
	}

	data16Bit = (UInt16*)((char*)data16Bit + imgRowBytes * top);

/* ---------------------- */
/* Loop over boxed pixels */
/* ---------------------- */

	for( v = top; v < bot; ++v,
		data16Bit = (UInt16*)((char*)data16Bit + imgRowBytes) ) {

		for( h = left; h < right; ++h ) {

			if( (x = data16Bit[h]) < nBins ) {
				++binArray[x];
				++N;
			}
			else {
				oflowS += x;
				++oflowC;
			}
		}
	}

exit:
	if( oFlowCnt )
		*oFlowCnt = oflowC;

	if( oFlowSum )
		*oFlowSum = oflowS;

	return N;
}


/* HSTUnitWidthUInt16 ------------------------------------------------
 *
 * Histogram 16-bit data array using fixed bin width of one,
 * i.e., values must be in range [0, nBins-1]. Values greater
 * than or equal to nBins are collected in the oFlow parameters
 * (if non-NULL).
 *
 * Return pixel count, excluding overflow.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * nData must be an even number.
 *
 * -----
 * Notes
 * -----
 *
 * 16-bit pixel pairs are read into a 32-bit register and extracted
 * separately. On X86 family processors, the 16-bit number at the
 * lower RAM address occupies the lower 16 bits of the register.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 HSTUnitWidthUInt16(
	UInt32			*oFlowCnt,
	UInt32			*oFlowSum,
	UInt32			*binArray,
	UInt32			nBins,
	const UInt16	*data16Bit,
	UInt32			nData )
{
	UInt32	*pPair;
	UInt32	i, pair, v0, v1, oflowC, oflowS;

/* ----------------- */
/* Initialize output */
/* ----------------- */

	MEMZeroBytes( binArray, nBins * sizeof(UInt32) );
	oflowC	= 0;
	oflowS	= 0;

/* ---------- */
/* Accumulate */
/* ---------- */

	nData >>= 1;
	pPair   = (UInt32*)data16Bit;

	for( i = 0; i < nData; ++i ) {

		pair	= *pPair++;
		v0		= pair & LoMask;
		v1		= pair >> HalfBits;

		if( v0 < nBins )
			++binArray[v0];
		else {
			oflowS += v0;
			++oflowC;
		}

		if( v1 < nBins )
			++binArray[v1];
		else {
			oflowS += v1;
			++oflowC;
		}
	}

	if( oFlowCnt )
		*oFlowCnt = oflowC;

	if( oFlowSum )
		*oFlowSum = oflowS;

	return nData - oflowC;
}


/* HSTGeneralFP32 ----------------------------------------------------
 *
 * Histogram FP32 data.
 *
 * There are two modes of operation:
 *
 * (1) Simple and fast mode
 * ------------------------
 * Set both uFlow and oFlow to NULL. In this mode [minH, maxH]
 * is interpreted as both the range of the bins and the range
 * of the data. No underflow or overflow checking is done. The
 * caller must guarantee that all data fall in this range.
 *
 * (2) Complex windowing mode
 * --------------------------
 * uFlow and oFlow must not be NULL. In this mode [minH, maxH]
 * is interpreted as the range of the bins, that is, only data
 * within the acceptance window [minH, maxH] are binned. Data
 * outside the range increase the underflow and overflow counters.
 *
 * Return highest bin count.
 *
 * ------------------
 * Data storage notes
 * ------------------
 *
 * The input values are assumed to reside at fixed offsets within
 * an array of records. The recBytes parameter should be set to
 * the record size. Note that pFirstVal should point to the VALUE
 * WITHIN the first record--this is NOT necessarily the start of
 * the record, itself.
 *
 * If the data are in a simple array of FP32, then recBytes should
 * be set to sizeof(FP32) and pFirstVal should point to the base
 * address of the array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int HSTGeneralFP32(
	FP32			*uFlow,
	FP32			*oFlow,
	UInt16			*binArray,
	UInt32			nBins,
	const FP32		*pFirstVal,
	UInt32			recBytes,
	UInt32			nData,
	FP32			minH,
	FP32			maxH )
{
	FP32	rng;
	int		i, count, maxV;

/* ----------------- */
/* Initialize output */
/* ----------------- */

	MEMZeroBytes( binArray, nBins * sizeof(UInt16) );

	rng = maxH - minH;

	if( !nData || rng <= 0.0F ) {
		binArray[0] = (int)nData;
		maxV		= nData;
		goto exit;
	}

/* ---------- */
/* Accumulate */
/* ---------- */

	--nBins;

	if( !uFlow ) {

		for( i = 0; i < (int)nData;
			++i,
			pFirstVal = (const FP32*)((char*)pFirstVal + recBytes) ) {

			++binArray[(int)((*pFirstVal - minH) * nBins / rng)];
		}
	}
	else {

		FP32	v;
		UInt32	uflow = 0, oflow = 0;

		for( i = 0; i < (int)nData;
			++i,
			pFirstVal = (const FP32*)((char*)pFirstVal + recBytes) ) {

			if( (v = *pFirstVal) < minH )
				++uflow;
			else if( v > maxH )
				++oflow;
			else
				++binArray[(int)((v - minH) * nBins / rng)];
		}

		*uFlow	= (FP32)uflow;
		*oFlow	= (FP32)oflow;
	}

	++nBins;

/* ---------------- */
/* Find highest bin */
/* ---------------- */

	maxV = 0;

	for( i = 0; i < (int)nBins; ++i ) {

		if( (count = binArray[i]) > maxV )
			maxV = count;
	}

exit:
	return maxV;
}


