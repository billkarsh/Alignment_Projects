/* BK_GEOM -----------------------------------------------------------
 *
 * Common geometry calculations.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */


#include	"BK_GEOM.h"






/* GEOMInsetU16Box1Pixel ---------------------------------------------
 *
 * Shrink box by 1 pixel on each side.
 *
 * -----
 * Notes
 * -----
 *
 * (1) srcBox must be valid. In particular, the bottom and right
 * must be non-zero.
 *
 * (2) Source and destination boxes can be same.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void GEOMInsetU16Box1Pixel( U16BoxPtr dstBox, U16BoxPtr srcBox )
{
    *(UInt32*)&dstBox->top = *(UInt32*)&srcBox->top + 0x10001;
    *(UInt32*)&dstBox->bottom = *(UInt32*)&srcBox->bottom - 0x10001;
}


/* GEOMOutsetU16Box --------------------------------------------------
 *
 * Enlarge all sides of box by given number of pixels.
 *
 * Left and right are confined to range [0,hLim].
 * Top and bottom are confined to range [0,vLim].
 *
 * hLim and vLim are assumed to be < (MAXUINT16 - outset).
 *
 * -----
 * Notes
 * -----
 *
 * Source and destination boxes can be same.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void GEOMOutsetU16Box(
    U16BoxPtr	dstBox,
    U16BoxPtr	srcBox,
    int			hLim,
    int			vLim,
    int			outset )
{
    int		temp;

    if( (temp = srcBox->top) > outset )
        dstBox->top = temp - outset;
    else
        dstBox->top = 0;

    if( (temp = srcBox->left) > outset )
        dstBox->left = temp - outset;
    else
        dstBox->left = 0;

    if( (temp = srcBox->bottom + outset) <= vLim )
        dstBox->bottom = temp;
    else
        dstBox->bottom = vLim;

    if( (temp = srcBox->right + outset) <= hLim )
        dstBox->right = temp;
    else
        dstBox->right = hLim;
}


/* GEOMU16PatchBox ---------------------------------------------------
 *
 * Calculate bounding box of patch that corresponds to srcBox.
 *
 * Left and right are confined to range [0,hLim].
 * Top and bottom are confined to range [0,vLim].
 *
 * -----
 * Notes
 * -----
 *
 * Source and destination boxes can be same.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void GEOMU16PatchBox(
    U16BoxPtr	dstBox,
    U16BoxPtr	srcBox,
    int			hLim,
    int			vLim )
{
    int		temp;

    dstBox->top    = srcBox->top;
    dstBox->left   = BitToWord( srcBox->left ) * WordBits;

    if( (temp = srcBox->bottom) <= vLim )
        dstBox->bottom = temp;
    else
        dstBox->bottom = vLim;

    if( (temp = srcBox->right) >= hLim )
        dstBox->right = hLim;
    else if( !temp )
        dstBox->right = 0;
    else
        dstBox->right = (BitToWord( temp - 1 ) + 1) * WordBits;
}


/* GEOMU16BoxIntersection --------------------------------------------
 *
 * Calculate intersection of two boxes.
 *
 * Return true if intersection non-empty.
 *
 * -----
 * Notes
 * -----
 *
 * outBox and box1 can be same.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

int GEOMU16BoxIntersection(
    U16BoxPtr		outBox,
    const U16BoxPtr	box1,
    const U16BoxPtr	box2 )
{
    int		temp;

    *outBox = *box1;

    if( (temp = box2->top) > box1->top )
        outBox->top = temp;

    if( (temp = box2->left) > box1->left )
        outBox->left = temp;

    if( (temp = box2->bottom) < box1->bottom )
        outBox->bottom = temp;

    if( (temp = box2->right) < box1->right )
        outBox->right = temp;

    return (outBox->top  < outBox->bottom &&
            outBox->left < outBox->right);
}


/* GEOMU16BoxUnion --------------------------------------------------
 *
 * Calculate union of two boxes.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

void GEOMU16BoxUnion(
    U16BoxPtr		outBox,
    const U16BoxPtr	box1,
    const U16BoxPtr	box2 )
{
    int		temp;

    *outBox = *box1;

    if( (temp = box2->top) < box1->top )
        outBox->top = temp;

    if( (temp = box2->left) < box1->left )
        outBox->left = temp;

    if( (temp = box2->bottom) > box1->bottom )
        outBox->bottom = temp;

    if( (temp = box2->right) > box1->right )
        outBox->right = temp;
}


/* GEOMPinSubrect ----------------------------------------------------
 *
 * Pin smaller RECT (r) to bounds of larger RECT (R).
 *
 * Return result as follows:
 *
 *	0:			No adjustments made.
 *	bit-0 set:	X adjusted.
 *	bit-1 set:	Y adjusted.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

#if defined(WIN32)

int GEOMPinSubrect( RECT *r, const RECT *R )
{
    int		dx, dy, moved = 0;

    dx	= r->right - r->left;
    dy	= r->bottom - r->top;

    if( r->left < R->left ) {
        r->left		= R->left;
        r->right	= R->left + dx;
        moved		= 1;
    }
    else if( r->right > R->right ) {
        r->right	= R->right;
        r->left		= R->right - dx;
        moved		= 1;
    }

    if( r->top < R->top ) {
        r->top		= R->top;
        r->bottom	= R->top + dy;
        moved		+= 2;
    }
    else if( r->bottom > R->bottom ) {
        r->bottom	= R->bottom;
        r->top		= R->bottom - dy;
        moved		+= 2;
    }

    return moved;
}

#endif	/* WIN32 */


/* GEOMSegLen --------------------------------------------------------
 *
 * Return binomial expansion of sqrt( dx^2 + dy^2 ).
 *
 * -----
 * Notes
 * -----
 *
 * Maximum error is about 5.6% for {dx,dy} <= 200.
 *
 * Copyright (c) 1999 Bill Karsh.
 * All rights reserved.
 *
 */

FP32 GEOMSegLen( int dx, int dy )
{
    if( dx < 0 )
        dx = -dx;

    if( dy < 0 )
        dy = -dy;

    if( dx < dy )
        return (dy + (dx * dx) / (2.0F * dy));
    else if( dy < dx )
        return (dx + (dy * dy) / (2.0F * dx));
    else
        return dx * 1.4142F;
}


