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

/*! @brief Declare functions for numerical schemes.
 *  @author Jianping Meng
 *  @details Declare functions for implementing various numerical
 *  schemes, including the kernels
 **/
#ifndef SCHEME_H
#define SCHEME_H
#include "flowfield.h"
#include "model.h"
#include "type.h"

// Define common stencils for implementing numerical schemes
/*!
 * LOCALSTENCIL: the current node, suitable for mainly the collision process
 */
extern ops_stencil LOCALSTENCIL;

/*!
 * ONEPTREGULARSTENCIL: the standard rectangular stencil, mainly suitable for a
 * first order finite difference scheme
 */
extern ops_stencil ONEPTREGULARSTENCIL;

/*!
 * TWOPTREGULARSTENCIL: the standard rectangular stencil, mainly suitable for a
 * second order finite difference scheme
 */
extern ops_stencil TWOPTREGULARSTENCIL;

/*!
 * ONEPTLATTICESTENCIL: the standard stencil for the stream scheme
 */
extern ops_stencil ONEPTLATTICESTENCIL;

void SetupCommonStencils();
// End Define common stencils

//HiLeMMS interface see https://gitlab.com/jpmeng/hilemms
void DefineScheme(const SchemeType scheme);
const int SchemeHaloNum();
void SetSchemeHaloNum(const int schemeHaloNum);
const SchemeType Scheme();
#ifdef OPS_3D
void Stream3D();
#endif //OPS_3D

#ifdef OPS_2D
void Stream();
#endif //OPS_2D
#endif
