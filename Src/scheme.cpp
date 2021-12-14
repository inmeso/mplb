
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

/*! @brief Define functions for numerical schemes.
 *  @author Jianping Meng
 *  @details Mainly define the functions that facilitating
 *  kernels for implementing numerical schemes
 */
#include "scheme.h"
#ifdef OPS_2D
int currentNode[] = {0, 0};
ops_stencil LOCALSTENCIL{ops_decl_stencil(2, 1, currentNode, "00")};
int rectangle[] = {0, 0, 1, 0, -1, 0, 0, 1, 0, -1};
ops_stencil ONEPTREGULARSTENCIL{
    ops_decl_stencil(2, 5, rectangle, "00:10:-10:01:0-1")};
int twoPtRectangle[] = {0,  0, 1, 0,  -1, 0, 0, 1, 0,
                        -1, 2, 0, -2, 0,  0, 2, 0, -2};
ops_stencil TWOPTREGULARSTENCIL{
    ops_decl_stencil(2, 9, twoPtRectangle, "00:10:-10:01:0-1:20:-20:02:0-2")};

int d2q9pts[] = {0, 0, 1, 0, -1, 0, 0, 1, 0, -1, 1, 1, -1, -1, 1, -1, -1, 1};
ops_stencil ONEPTLATTICESTENCIL{
    ops_decl_stencil(2, 9, d2q9pts, "00:10:-10:01:0-1:11:-1-1:1-1:-11")};
#endif /* OPS_2D */
#ifdef OPS_3D
int currentNode[]{0, 0, 0};
ops_stencil LOCALSTENCIL{ops_decl_stencil(3, 1, currentNode, "000")};
int rectangle[]{0, 0, 0,  1, 0, 0, -1, 0, 0, 0, 1,
                0, 0, -1, 0, 0, 0, 1,  0, 0, -1};
ops_stencil ONEPTREGULARSTENCIL{
    ops_decl_stencil(3, 7, rectangle, "000:100:-100:010:0-10:001:00-1")};
// ops_printf("%s\n", "ONEPTREGULAR finished!");
int twoPtRectangle[]{0,  0, 0, 1, 0, 0,  -1, 0,  0, 0, 1, 0,  0,
                     -1, 0, 0, 0, 1, 0,  0,  -1, 2, 0, 0, -2, 0,
                     0,  0, 2, 0, 0, -2, 0,  0,  0, 2, 0, 0,  -2};
ops_stencil TWOPTREGULARSTENCIL{ops_decl_stencil(
    3, 13, twoPtRectangle,
    "000:100:-100:010:0-10:001:00-1:200:-200:020:0-20:002:00-2")};
int d3q27[] = {-1, -1, -1, -1, -1, 0,  -1, -1, 1, -1, 0,  -1, -1, 0,  0,  -1, 0,
               1,  -1, 1,  -1, -1, 1,  0,  -1, 1, 1,  0,  -1, -1, 0,  -1, 0,  0,
               -1, 1,  0,  0,  -1, 0,  0,  0,  0, 0,  1,  0,  1,  -1, 0,  1,  0,
               0,  1,  1,  1,  -1, -1, 1,  -1, 0, 1,  -1, 1,  1,  0,  -1, 1,  0,
               0,  1,  0,  1,  1,  1,  -1, 1,  1, 0,  1,  1,  1};
ops_stencil ONEPTLATTICESTENCIL{ops_decl_stencil(3, 27, d3q27, "D3Q27")};
#endif /* OPS_3D */

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
int schemeHaloPt{1};
SchemeType schemeType{Scheme_StreamCollision};
SchemeType Scheme() { return schemeType; }


void DefineScheme(const SchemeType scheme, const SizeType timeStep) {
    schemeType = scheme;
    switch (schemeType) {
        case Scheme_StreamCollision: {
            SetSchemeHaloNum(1);
            if (timeStep == 0) {
                g_f().CreateFieldFromScratch(g_Block());
            } else {
                g_f().CreateFieldFromFile(CaseName(), g_Block(), timeStep);
            }
            g_fStage().CreateFieldFromScratch(g_Block());
            RegisterFieldNeedHalo(g_fStage());
            ops_printf("The stream-collision scheme is chosen!\n");
        } break;
        case Scheme_StreamCollision_Swap: {
            SetSchemeHaloNum(1);
            if (timeStep == 0) {
                g_f().CreateFieldFromScratch(g_Block());
            } else {
                g_f().CreateFieldFromFile(CaseName(), g_Block(), timeStep);
            }
            RegisterFieldNeedHalo(g_f());
            ops_printf("The stream-collision_swap scheme is chosen!\n");
        } break;
        default:
            break;
    }
}
int SchemeHaloNum() { return schemeHaloPt; }
void SetSchemeHaloNum(const int schemeHaloNum) { schemeHaloPt = schemeHaloNum; }
