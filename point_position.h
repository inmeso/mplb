#ifndef POINT_POSITION_H
#define POINT_POSITION_H
#include "type.h"
/*
 * Utilties for comparing real numbers
 * might be moved to a source file specifically for utilites in the future
 */

/*!
 * @fn judging if a point is inside a polygon or not
 * @param point the coordinate of the point
 * @param polygon the coordinates of the polygon vertex
 * @param polyVertexNum the total number of polygon vertex
 */
PointPosition IfPointInPoly(const Real* point, const Real* polygon,
                            const long long polyVertexNum);
#endif  // POINT_POSITION_H
