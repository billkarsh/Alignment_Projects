/* BK_STAT -----------------------------------------------------------
 *
 * Common statistics calculations.
 *
 * ---------------------
 * Propagation of errors
 * ---------------------
 *
 * The following recounts how to combine errors in common
 * expressions...
 *
 * Notation
 * --------
 * p, q	= independent variables
 * a, b = unsigned constants
 * x	= dependent variable
 * V	= variance
 *
 * Addition and subtraction
 * ------------------------
 *
 * x = ap +/- bq
 *
 * Vx = [aaVp + bbVq] +/- 2abVpq
 *
 * Multiplication
 * --------------
 *
 * x = +/- apq
 *
 * Vx/(xx) = Vp/(pp) + Vq/(qq) + 2Vpq/(pq)
 *
 * This is calculated efficiently as:
 *
 * Vx = aa[qqVp + ppVq + 2pqVpq]
 *
 * Division
 * --------
 *
 * x = +/- ap/q
 *
 * Vx/(xx) = Vp/(pp) + Vq/(qq) - 2Vpq/(pq)
 *
 * This is calculated efficiently as:
 *
 * Vx = [aaVp + xxVq - 2axVpq] / qq
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_MEM.h"
#include	"BK_SORT_FP32.h"
#include	"BK_STAT.h"

#include	<math.h>






/* STATAverageUInt16 -------------------------------------------------
 *
 * Return average of array of N 16-bit values.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * N must be an even number.
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

#if defined(WIN32)

int STATAverageUInt16( const UInt16 *array, UInt32 N )
{
    UInt32	*pPair;
    UInt32	i, sum, pair;

    sum		= 0;
    pPair	= (UInt32*)array;

    for( i = 0; i < N; i += 2 ) {

        pair = *pPair++;
        sum += pair & LoMask;
        sum += pair >> HalfBits;
    }

    return sum / N;
}

#endif


/* STATIndexOfMaxValueUInt16 -----------------------------------------
 *
 * Return index of maximum value in array of 16-bit numbers.
 *
 * ---------
 * IMPORTANT
 * ---------
 *
 * N must be an even number.
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

UInt32 STATIndexOfMaxValueUInt16( const UInt16 *array, UInt32 N )
{
    UInt32	*pPair;
    UInt32	i, idx, pair;
    int		t, max;

    max		= 0;
    idx		= 0;
    pPair	= (UInt32*)array;

    for( i = 0; i < N; i += 2 ) {

        pair = *pPair++;

        t = pair & LoMask;

        if( t > max ) {
            max	= t;
            idx	= i;
        }

        t = pair >> HalfBits;

        if( t > max ) {
            max	= t;
            idx	= i + 1;
        }
    }

    return idx;
}


/* STATIndexOfMaxValueUInt32 -----------------------------------------
 *
 * Return index of maximum value in array of 32-bit numbers.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 STATIndexOfMaxValueUInt32( const UInt32 *array, UInt32 N )
{
    UInt32	i, idx, t, max;

    max		= 0;
    idx		= 0;

    for( i = 0; i < N; ++i ) {

        if( (t = array[i]) > max ) {
            max	= t;
            idx	= i;
        }
    }

    return idx;
}


/* STATPercentileBinUInt32 -------------------------------------------
 *
 * Return index of bin such that integrated count from
 * bin zero through result bin is indicated fraction of
 * total counts.
 *
 * For example, if frac = 0.5, the result is the median bin.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

UInt32 STATPercentileBinUInt32(
    const UInt32	*hist,
    UInt32			nBins,
    UInt32			nCounts,
    FP64			frac )
{
    UInt32	i, sum, v, bLast;

    nCounts	= (UInt32)(nCounts * frac);
    sum		= 0;
    bLast	= 0;

    for( i = 0; i < nBins; ++i ) {

        if( v = hist[i] ) {

            if( (sum += v) > nCounts )
                return (bLast + i) >> 1;
            else if( sum == nCounts )
                return i;

            bLast = i;
        }
    }

    return nBins - 1;
}


/* STATSigma32 -------------------------------------------------------
 *
 * Return sigma(X) = sqrt( 1/(N-1) * SUM(x - meanX)^2 ).
 *
 * If sumX2 = SUM(X^2), then
 *
 * sigma(X) = sqrt( 1/(N-1) * [sumX2 - N * meanX^2] );
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 STATSigma32( FP32 N, FP32 meanX, FP32 sumX2 )
{
    FP32	sigma = 0.0F;

    if( N > 1.0F ) {

        sigma = (sumX2 - N * meanX * meanX) / (N - 1.0F);
        sigma = (FP32)sqrt( sigma );
    }

    return sigma;
}


/* STATAveFP32 -------------------------------------------------------
 *
 * Return average of list of FP32.
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

FP32 STATAveFP32(
    const FP32		*pFirstVal,
    UInt32			recBytes,
    UInt32			N )
{
    FP32	ave = 0.0F;
    UInt32	i;

    if( N ) {

        for( i = 0; i < N;
            ++i,
            pFirstVal = (const FP32*)((char*)pFirstVal + recBytes) ) {

            ave += *pFirstVal;
        }

        ave /= N;
    }

    return ave;
}


/* STATMinMaxFP32 ----------------------------------------------------
 *
 * Return maximum and (fill in) minimum values in list of FP32.
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

FP32 STATMinMaxFP32(
    FP32			*minVal,
    const FP32		*pFirstVal,
    UInt32			recBytes,
    UInt32			N )
{
    FP32	v, lMin, lMax;
    UInt32	i;

    if( N ) {

        lMin = lMax = *pFirstVal;

        pFirstVal = (const FP32*)((char*)pFirstVal + recBytes);

        for( i = 1; i < N;
            ++i,
            pFirstVal = (const FP32*)((char*)pFirstVal + recBytes) ) {

            if( (v = *pFirstVal) > lMax )
                lMax = v;
            else if( v < lMin )
                lMin = v;
        }
    }
    else {
        lMin = 0;
        lMax = 0;
    }

    *minVal = lMin;

    return lMax;
}


/* STATAveAndSDInRangeFP32 -------------------------------------------
 *
 * Calculate the sample average and standard deviation for a
 * range of data points: [i0,...,iLim).
 *
 * Return the average, std dev is filled into stdDev parameter.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 STATAveAndSDInRangeFP32(
    FP32			*stdDev,
    const FP32		*array,
    UInt32			i0,
    UInt32			iLim )
{
    FP32	x, sumX = 0.0F, sumX2 = 0.0F;
    UInt32	i, n = 0;

    for( i = i0; i < iLim; ++i, ++n ) {

        sumX	+= (x = array[i]);
        sumX2	+= x * x;
    }

    x		= (FP32)n;
    *stdDev = STATSigma32( x, sumX /= x, sumX2 );

    return sumX;
}


/* STATWindowedMaxFP32 -----------------------------------------------
 *
 * Given an array of data with elements [i0,iLim) and a specified
 * windowSize...
 *
 * Calculate the average value of the array in every window to
 * find the maximum value window.
 *
 * This is essentially a smoothing method that can help identify
 * the peak of a sampled curve when the data contain significant
 * fluctuations.
 *
 * Return the average array value over the maximum value window;
 * that is, the approximate maximum value of the curve.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 STATWindowedMaxFP32(
    const FP32		*array,
    UInt32			i0,
    UInt32			iLim,
    UInt32			windowSize )
{
    FP32	max, win;
    UInt32	iRMost;		/* right-most point in window */

    if( iLim - i0 > windowSize ) {

        /* ---------------------------------- */
        /* Construct first window and set max */
        /* ---------------------------------- */

        for( iRMost = i0 + 1, win = array[i0];
             iRMost < i0 + windowSize; ++iRMost ) {

             win += array[iRMost];
        }

        max = win;

        /* ------------------------------------------------------ */
        /* Slide window, remove old left-most, add new right-most */
        /* ------------------------------------------------------ */

        for( ; iRMost < iLim; ++iRMost ) {

            win -= array[iRMost - windowSize];
            win += array[iRMost];

            if( win > max )
                max = win;
        }

        /* ----------------- */
        /* Report as average */
        /* ----------------- */

        max /= windowSize;
    }
    else {

        /* --------------------------------- */
        /* Data span smaller than one window */
        /* --------------------------------- */

        max = STATMinMaxFP32( &win, array + i0,
                sizeof(FP32), iLim - i0 );
    }

    return max;
}


/* STATMedianFP32 ----------------------------------------------------
 *
 * Return the median value of an FP32 array.
 *
 * array
 * -----
 * The array of values.
 *
 * wrkspc
 * ------
 * This routine involves sorting the input values. If you set
 * wrkspc to NULL, then the array, itself, will be sorted (in
 * ascending order). If you point wrkspc to a non-NULL array
 * of size N, then wrkspc will be used to sort a copy of array.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 STATMedianFP32( FP32 *array, FP32 *wrkspc, UInt32 N )
{
    FP32	median;

/* -------------------- */
/* Copy array if needed */
/* -------------------- */

    if( wrkspc )
        MEMCopyBytes( wrkspc, array, N * sizeof(FP32) );
    else
        wrkspc = array;

/* ---- */
/* Sort */
/* ---- */

    SORT_FP32_Ascending( wrkspc, N );

/* ------------------- */
/* Select middle value */
/* ------------------- */

    if( N & 1 )
        median = wrkspc[N / 2 + 1];
    else {
        N /= 2;
        median = (wrkspc[N] + wrkspc[N + 1]) * 0.5F;
    }

    return median;
}


/* STATUnitWeightsFP32 -----------------------------------------------
 *
 * Set FP32 fitting weights to unity.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void STATUnitWeightsFP32( FP32 *wt, UInt32 N )
{
    UInt32	i;

    for( i = 0; i < N; ++i )
        wt[i] = 1.0F;
}


/* STATSigmaMADFP32 --------------------------------------------------
 *
 * Given an array (r) of signed residuals from a curve fit,
 * calculate the Median Absolute Deviation (sigma):
 *
 *			median( |Ri - median( Ri )| )
 *	sigma =	-----------------------------
 *					0.6745
 *
 * wrkspc
 * ------
 * You must set wrkspc to the non-NULL address of
 * an FP32 scratch buffer of size N.
 *
 * Return sigma.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 STATSigmaMADFP32( FP32 *r, FP32 *wrkspc, UInt32 N )
{
    FP32	median;
    UInt32	i;

/* ---------------------- */
/* Calculate median( Ri ) */
/* ---------------------- */

    median = STATMedianFP32( r, wrkspc, N );

/* -------------------------- */
/* Create array |Ri - median| */
/* -------------------------- */

    for( i = 0; i < N; ++i )
        wrkspc[i] = (FP32)fabs( wrkspc[i] - median );

/* --------------------------------- */
/* Calculate median( |Ri - median| ) */
/* --------------------------------- */

    return STATMedianFP32( wrkspc, NULL, N ) / 0.6745F;
}


/* STATTukeyWtsFP32 --------------------------------------------------
 *
 * Given an array (r) of signed residuals from a curve fit,
 * calculate a new array (w) of Tukey biweights:
 *
 *	Wi = [1 - (Ri/B.s)^2]^2		; |Ri| <  B.s
 *
 *	Wi = 0						; |Ri| >= B.s,
 *
 *	where, s = sigma(MAD) (see STATSigmaMADFP32), and B is
 *	constant tuning parameter with value 4.685.
 *
 * wrkspc
 * ------
 * You must set wrkspc to the non-NULL address of
 * an FP32 scratch buffer of size N.
 *
 * Return count of non-zero weights.
 *
 * -----
 * Notes
 * -----
 *
 * Output array w[] and input array r[] can be the same if desired.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int STATTukeyWtsFP32( FP32 *w, FP32 *r, FP32 *wrkspc, UInt32 N )
{
    FP32	Bs2;
    int		nnz;

    Bs2 = STATSigmaMADFP32( r, wrkspc, N ) * 4.685F;
    nnz = 0;

    if( Bs2 <= 0.0F )
        MEMZeroBytes( w, N * sizeof(FP32) );
    else {

        FP32	invBs2, r2, x;
        UInt32	i;

        Bs2	   *= Bs2;
        invBs2	= 1.0F / Bs2;

        for( i = 0; i < N; ++i ) {

            x = r[i];

            if( (r2 = x * x) < Bs2 ) {

                x	 = (1.0F - r2 * invBs2);
                ++nnz;
                w[i] = x * x;
            }
            else
                w[i] = 0.0F;
        }
    }

    return nnz;
}


/* STATBisquareWtsFP32 -----------------------------------------------
 *
 * Given an array (r) of signed residuals from a curve fit,
 * calculate a new array (w) of "bisquare" weights. This is
 * a common variation on Tukey biweights that only requires
 * one calculation of a median:
 *
 *	Wi = [1 - (Ri/B.s)^2]^2		; |Ri| <  B.s
 *
 *	Wi = 0						; |Ri| >= B.s,
 *
 *	where, s = median( |Ri| ), and B is constant
 *	tuning parameter with value 6.0.
 *
 * wrkspc
 * ------
 * You must set wrkspc to the non-NULL address of
 * an FP32 scratch buffer of size N.
 *
 * Return count of non-zero weights.
 *
 * -----
 * Notes
 * -----
 *
 * Output array w[] and input array r[] can be the same if desired.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int STATBisquareWtsFP32( FP32 *w, FP32 *r, FP32 *wrkspc, UInt32 N )
{
    FP32	Bs2;
    UInt32	i;
    int		nnz;

    for( i = 0; i < N; ++i )
        wrkspc[i] = (FP32)fabs( r[i] );

    Bs2	= STATMedianFP32( wrkspc, NULL, N ) * 6.0F;
    nnz	= 0;

    if( Bs2 <= 0.0F )
        MEMZeroBytes( w, N * sizeof(FP32) );
    else {

        FP32	invBs2, r2, x;

        Bs2	   *= Bs2;
        invBs2	= 1.0F / Bs2;

        for( i = 0; i < N; ++i ) {

            x = r[i];

            if( (r2 = x * x) < Bs2 ) {

                x	 = (1.0F - r2 * invBs2);
                ++nnz;
                w[i] = x * x;
            }
            else
                w[i] = 0.0F;
        }
    }

    return nnz;
}


/* STATSlopeFit32XY --------------------------------------------------
 *
 * Fit data pairs (x,y) to line through origin on interval [i0,iLim).
 *
 * Calculate: slope, linear correlation coefficient.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void STATSlopeFit32XY(
    STATLinFitPtr	L,
    const FP32		*x,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim )
{
    FP64	X, Y, r, varX, varY,
            n		= 0.0,
            sumX	= 0.0,
            sumY	= 0.0,
            sumXX	= 0.0,
            sumYY	= 0.0,
            sumXY	= 0.0;
    UInt32	i;

/* --------------------------- */
/* Collect sums over [i0,iLim) */
/* --------------------------- */

    for( i = i0; i < iLim; ++i ) {

        ++n;
        X	= x[i];
        Y	= y[i];
        sumX	+= X;
        sumY	+= Y;
        sumXX	+= X * X;
        sumYY	+= Y * Y;
        sumXY	+= X * Y;
    }

/* --------------- */
/* Compute results */
/* --------------- */

    r		= n * sumXY - sumX * sumY;
    varX	= n * sumXX - sumX * sumX;	/* a.k.a. delta */
    varY	= n * sumYY - sumY * sumY;

    L->icpt		= 0;
    L->slope	= sumXY / sumXX;
    L->linCor	= r / sqrt( varX * varY );
}


/* STATLinFit32XY ----------------------------------------------------
 *
 * Fit data pairs (x,y) to line on interval [i0,iLim).
 *
 * Calculate: intercept, slope, linear correlation coefficient.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void STATLinFit32XY(
    STATLinFitPtr	L,
    const FP32		*x,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim )
{
    FP64	X, Y, r, varX, varY,
            n		= 0.0,
            sumX	= 0.0,
            sumY	= 0.0,
            sumXX	= 0.0,
            sumYY	= 0.0,
            sumXY	= 0.0;
    UInt32	i;

/* --------------------------- */
/* Collect sums over [i0,iLim) */
/* --------------------------- */

    for( i = i0; i < iLim; ++i ) {

        ++n;
        X	= x[i];
        Y	= y[i];
        sumX	+= X;
        sumY	+= Y;
        sumXX	+= X * X;
        sumYY	+= Y * Y;
        sumXY	+= X * Y;
    }

/* --------------- */
/* Compute results */
/* --------------- */

    r		= n * sumXY - sumX * sumY;
    varX	= n * sumXX - sumX * sumX;	/* a.k.a. delta */
    varY	= n * sumYY - sumY * sumY;

    L->icpt		= (sumXX * sumY - sumX * sumXY) / varX;
    L->slope	= r / varX;
    L->linCor	= r / sqrt( varX * varY );
}


/* STATLinFit64XY ----------------------------------------------------
 *
 * Fit data pairs (x,y) to line on interval [i0,iLim).
 *
 * Calculate: intercept, slope, linear correlation coefficient.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void STATLinFit64XY(
    STATLinFitPtr	L,
    const FP64		*x,
    const FP64		*y,
    UInt32			i0,
    UInt32			iLim )
{
    FP64	X, Y, r, varX, varY,
            n		= 0.0,
            sumX	= 0.0,
            sumY	= 0.0,
            sumXX	= 0.0,
            sumYY	= 0.0,
            sumXY	= 0.0;
    UInt32	i;

/* --------------------------- */
/* Collect sums over [i0,iLim) */
/* --------------------------- */

    for( i = i0; i < iLim; ++i ) {

        ++n;
        X	= x[i];
        Y	= y[i];
        sumX	+= X;
        sumY	+= Y;
        sumXX	+= X * X;
        sumYY	+= Y * Y;
        sumXY	+= X * Y;
    }

/* --------------- */
/* Compute results */
/* --------------- */

    r		= n * sumXY - sumX * sumY;
    varX	= n * sumXX - sumX * sumX;	/* a.k.a. delta */
    varY	= n * sumYY - sumY * sumY;

    L->icpt		= (sumXX * sumY - sumX * sumXY) / varX;
    L->slope	= r / varX;
    L->linCor	= r / sqrt( varX * varY );
}


/* STATLinFit32Y -----------------------------------------------------
 *
 * Fit data pairs (x0+k*dx,y) to line on interval [i0,iLim).
 *
 * Calculate: intercept, slope, linear correlation coefficient.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void STATLinFit32Y(
    STATLinFitPtr	L,
    FP32			x0,
    FP32			dx,
    const FP32		*y,
    UInt32			i0,
    UInt32			iLim )
{
    FP64	X, Y, r, varX, varY,
            n		= 0.0,
            sumX	= 0.0,
            sumY	= 0.0,
            sumXX	= 0.0,
            sumYY	= 0.0,
            sumXY	= 0.0;
    UInt32	i;

/* --------------------------- */
/* Collect sums over [i0,iLim) */
/* --------------------------- */

    for( i = i0, X = x0; i < iLim; ++i, X += dx ) {

        ++n;
        Y	= y[i];
        sumX	+= X;
        sumY	+= Y;
        sumXX	+= X * X;
        sumYY	+= Y * Y;
        sumXY	+= X * Y;
    }

/* --------------- */
/* Compute results */
/* --------------- */

    r		= n * sumXY - sumX * sumY;
    varX	= n * sumXX - sumX * sumX;	/* a.k.a. delta */
    varY	= n * sumYY - sumY * sumY;

    L->icpt		= (sumXX * sumY - sumX * sumXY) / varX;
    L->slope	= r / varX;
    L->linCor	= r / sqrt( varX * varY );
}


