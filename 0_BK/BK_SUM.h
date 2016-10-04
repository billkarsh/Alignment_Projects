/* BK_SUM.h ----------------------------------------------------------
 *
 * Sums over masked images.
 *
 * -----------
 * Definitions
 * -----------
 *
 * Intensity (I)	= power/unit area
 *
 * IOD				= "integrated object density"
 *
 * Flux (F)			= IOD
 *					= Integral of (I . dA)
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#pragma once


#include	"BK_DEFS_GEN.h"


#if defined(__cplusplus)
extern "C" {
#endif


/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	SUM_RecIntensity( recptr )								\
    ((recptr)->A ? (FP32)(recptr)->F / (recptr)->A : 0.0F)

/* --------------------------------------------------------------- */
/* Data structures ----------------------------------------------- */
/* --------------------------------------------------------------- */

#if defined(WIN32)
#pragma pack( push, PACK_S )
#pragma pack( 1 )
#elif defined(__linux__)
#pragma pack( push, 1 )
#pragma pack( 1 )
#endif

typedef struct {

    UInt32	F,		/* flux	*/
            A;		/* area	*/

} SUM_F_A_Rec, *SUM_F_A_Ptr;

typedef struct {

    UInt32	F,		/* flux					*/
            F2,		/* sum(flux squared)	*/
            A;		/* area					*/

} SUM_F_F2_A_Rec, *SUM_F_F2_A_Ptr;

typedef struct {

    UInt32	F,		/* flux				*/
            A,		/* area				*/
            Ipk;	/* peak intensity	*/

} SUM_F_A_Ipk_Rec, *SUM_F_A_Ipk_Ptr;

typedef struct {

    SUM_F_A_Rec		obj,
                    bkg;

} SUM_Obj_Bkg_Rec, *SUM_Obj_Bkg_Ptr;

#if defined(WIN32)
#pragma pack( pop, PACK_S )
#elif defined(__linux__)
#pragma pack( pop )
#endif






/* --------------------------------------------------------------- */
/* Whole map ----------------------------------------------------- */
/* --------------------------------------------------------------- */

void SUMMap_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix );

void SUMMap_Obj_Bkg(
    SUM_Obj_Bkg_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*objMap,
    const UInt32		*bkgMap,
    UInt32				hPix,
    UInt32				vPix );

/* --------------------------------------------------------------- */
/* Blobs --------------------------------------------------------- */
/* --------------------------------------------------------------- */

void SUMBlob_F_A(
    SUM_F_A_Ptr			sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

void SUMBlob_F_A_Ipk(
    SUM_F_A_Ipk_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*map,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );

void SUMBlob_Obj_Bkg(
    SUM_Obj_Bkg_Ptr		sum,
    const UInt16		*data16Bit,
    const UInt32		*objMap,
    const UInt32		*bkgMap,
    UInt32				hPix,
    UInt32				vPix,
    const U16BoxPtr		rBlob );


#if defined(__cplusplus)
}
#endif
