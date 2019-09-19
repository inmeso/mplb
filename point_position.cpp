/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
*/

/*! @brief utilities for judging the position of a point
 *  @author Jianping Meng
 *  @details judging the position of a point relative to a polygon or polyhedron
 *  using ray-crossing scheme
 */
#include "point_position.h"

PointPosition IfPointInPoly(const Real* point, const Real* polygon,
                            const long long polyVertexNum) {
    /*This routine implements the algorithm discussed in Joseph O' Rourke,
     *Computational Geometry in C, 1998, Cambridge University Press.
     *Two dimensional version
    */
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
