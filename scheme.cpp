// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief Define functions for numerical schemes.
 *  @author Jianping Meng
 *  @details Mainly define the functions that facilitating
 *  kernels for implementing numerical schemes
 */
#include "scheme.h"
ops_stencil LOCALSTENCIL;
ops_stencil ONEPTREGULARSTENCIL;
ops_stencil ONEPTLATTICESTENCIL;
ops_stencil TWOPTREGULARSTENCIL;
/*!
 * A numerical scheme may need two kinds of halo points:
 * 1. Halo points between two blocks, we may call them connection boundary.
 * 2. Halo points for implementing high order schemes, which needs extra halo
 *    points. For example, a second order scheme needs one extra point, which
 *    has the information for implementing a actual first order scheme, we may
 *    call them virtual boundary.
 * On the other hand, there is some difference between two kinds of boundary
 * conditions, for a second order scheme, 1 needs two halo points but 2 only
 * need one. We will use the larger one here.
 */
int schemeHaloPt = 1;
SchemeType schemeType{Scheme_StreamCollision};
const SchemeType Scheme() { return schemeType; }
void SetupCommonStencils() {
#ifdef OPS_2D
    int currentNode[] = {0, 0};
    LOCALSTENCIL = ops_decl_stencil(2, 1, currentNode, "00");
    int rectangle[] = {0, 0, 1, 0, -1, 0, 0, 1, 0, -1};
    ONEPTREGULARSTENCIL = ops_decl_stencil(2, 5, rectangle, "00:10:-10:01:0-1");
    int twoPtRectangle[] = {0,  0, 1, 0,  -1, 0, 0, 1, 0,
                            -1, 2, 0, -2, 0,  0, 2, 0, -2};
    TWOPTREGULARSTENCIL = ops_decl_stencil(2, 9, twoPtRectangle,
                                           "00:10:-10:01:0-1:20:-20:02:0-2");

    int d2q9[] = {0, 0, 1, 0, -1, 0, 0, 1, 0, -1, 1, 1, -1, -1, 1, -1, -1, 1};
    ONEPTLATTICESTENCIL =
        ops_decl_stencil(2, 9, d2q9, "00:10:-10:01:0-1:11:-1-1:1-1:-11");
#endif /* OPS_2D */
#ifdef OPS_3D
    int currentNode[]{0, 0, 0};
    LOCALSTENCIL = ops_decl_stencil(3, 1, currentNode, "000");
    int rectangle[]{0, 0, 0,  1, 0, 0, -1, 0, 0, 0, 1,
                    0, 0, -1, 0, 0, 0, 1,  0, 0, -1};
    ONEPTREGULARSTENCIL =
        ops_decl_stencil(3, 7, rectangle, "000:100:-100:010:0-10:001:00-1");
    // ops_printf("%s\n", "ONEPTREGULAR finished!");
    int twoPtRectangle[]{0,  0, 0, 1, 0, 0,  -1, 0,  0, 0, 1, 0,  0,
                         -1, 0, 0, 0, 1, 0,  0,  -1, 2, 0, 0, -2, 0,
                         0,  0, 2, 0, 0, -2, 0,  0,  0, 2, 0, 0,  -2};
    TWOPTREGULARSTENCIL = ops_decl_stencil(
        3, 13, twoPtRectangle,
        "000:100:-100:010:0-10:001:00-1:200:-200:020:0-20:002:00-2");
    int d3q27[] = {0,  0,  0,  -1, -1, -1, -1, -1, 0,  -1, -1, 1, -1, 0,
                   -1, -1, 0,  0,  -1, 0,  1,  -1, 1,  -1, -1, 1, 0,  -1,
                   1,  1,  0,  -1, -1, 0,  -1, 0,  0,  -1, 1,  0, 0,  -1,
                   0,  0,  0,  0,  0,  1,  0,  1,  -1, 0,  1,  0, 0,  1,
                   1,  1,  -1, -1, 1,  -1, 0,  1,  -1, 1,  1,  0, -1, 1,
                   0,  0,  1,  0,  1,  1,  1,  -1, 1,  1,  0,  1, 1,  1};
    ONEPTLATTICESTENCIL = ops_decl_stencil(3, 28, d3q27, "D3Q27");
#endif /* OPS_3D */
}

void DefineScheme(const SchemeType scheme) {
    schemeType = scheme;
    SetupCommonStencils();
    switch (schemeType) {
        case Scheme_StreamCollision: {
            SetSchemeHaloNum(1);
            ops_printf("The stream-collision scheme is chosen!\n");
        } break;
        default:
            break;
    }
}
const int SchemeHaloNum() { return schemeHaloPt; }
void SetSchemeHaloNum(const int schemeHaloNum) { schemeHaloPt = schemeHaloNum; }
#include "scheme_kernel.h"
