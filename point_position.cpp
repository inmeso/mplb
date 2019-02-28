// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief utilities for judging the position of a point
 *  @author Jianping Meng
 *  @details judging the position of a point relative to a polygon or polyhedron
 *  using ray-crossing scheme
 */
#include "point_position.h"

PointPosition IfPointInPoly(const Real* point, const Real* polygon,
                            const long long polyVertexNum) {
    long long currentVertex, previousVertex;
    Real x;
    long long rightRayCross = 0;
    long long leftRayCross = 0;
    const int DIM = 2;
    bool rStrad, lStrad;
    /* For each edge e=(i-1,i), see if crosses ray. */
    for (currentVertex = 0; currentVertex < polyVertexNum; currentVertex++) {
        // First see if the point is a vertex.
        const Real xDiffCurrent{polygon[DIM * currentVertex] - point[0]};
        const Real yDiffCurrent{polygon[DIM * currentVertex + 1] - point[1]};
        if (EssentiallyEqual(&xDiffCurrent, &ZERO, EPS) &&
            EssentiallyEqual(&yDiffCurrent, &ZERO, EPS)) {
            return IsVertex;
        }
        // calculate the index of previous point
        previousVertex = (currentVertex + polyVertexNum - 1) % polyVertexNum;
        const Real xDiffPrevious{polygon[DIM * previousVertex] - point[0]};
        const Real yDiffPrevious{polygon[DIM * previousVertex + 1] - point[1]};
        rStrad = (DefinitelyGreaterThan(&yDiffCurrent, &ZERO, EPS) !=
                  DefinitelyGreaterThan(&yDiffPrevious, &ZERO, EPS));
        lStrad = (DefinitelyLessThan(&yDiffCurrent, &ZERO, EPS) !=
                  DefinitelyLessThan(&yDiffPrevious, &ZERO, EPS));
        if (rStrad || lStrad) {
            /* e straddles ray, so compute intersection with ray. */
            x = (xDiffCurrent * yDiffPrevious - xDiffPrevious * yDiffCurrent) /
                (yDiffPrevious - yDiffCurrent);
            if (rStrad && DefinitelyGreaterThan(&x, &ZERO, EPS)) {
                rightRayCross++;
            }
            if (lStrad && DefinitelyLessThan(&x, &ZERO, EPS)) {
                leftRayCross++;
            }
        }
    }
    // point on the edge if left and right cross are not the same parity.
    if ((rightRayCross % 2) != (leftRayCross % 2)) {
        return RelativelyInteriorToEdge;
    }
    // point inside iff an odd number of crossings.
    if ((rightRayCross % 2) == 1) {
        return StrictlyInterior;
    } else {
        return StrictlyExterior;
    }
}
