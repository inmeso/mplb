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

/*! @brief   Head files for boundary conditions
 * @author  Jianping Meng
 * @details Declaring functions related to boundary conditions.
 */
#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "type.h"
#include <vector>
#include "block.h"
#include "flowfield_host_device.h"
#include "boundary_host_device.h"
#include "model.h"
#include "model_host_device.h"

const std::vector<BoundarySurface> AllBoundarySurface{BoundarySurface::Left,
                                   BoundarySurface::Right,
                                   BoundarySurface::Top,
                                   BoundarySurface::Bottom,
#ifdef OPS_3D
                                   BoundarySurface::Front,
                                   BoundarySurface::Back,
                                   BoundarySurface::LeftBack,
                                   BoundarySurface::LeftFront,
                                   BoundarySurface::RightBack,
                                   BoundarySurface::RightFront,
                                   BoundarySurface::TopBack,
                                   BoundarySurface::TopFront,
                                   BoundarySurface::BottomBack,
                                   BoundarySurface::BottomFront,
#endif
                                   BoundarySurface::LeftTop,
                                   BoundarySurface::LeftBottom,
                                   BoundarySurface::RightTop,
                                   BoundarySurface::RightBottom

};

enum class BoundaryScheme {
    KineticDiffuseWall = 11,
    KineticSpelluarWall = 12,
    ExtrapolPressure1ST = 16,
    ExtrapolPressure2ND = 17,
    MDPeriodic = 18,
    FDPeriodic = 19,
    BounceBack = 20,
    FreeFlux = 21,
    ZouHeVelocity = 22,
    EQNNoSlip = 23,
    EQMDiffuseRefl = 24,
    None = -1
};

struct BlockBoundary {
    int blockIndex;
    int componentID;
    std::vector<Real> givenVars;
    BoundarySurface boundarySurface;
    BoundaryScheme boundaryScheme;
    std::vector<VariableTypes> macroVarTypesatBoundary;
    VertexType boundaryType;
};


int BoundaryHaloNum();
void SetBoundaryHaloNum(const int boundaryhaloNum);
void BoundaryNormal3D(const VertexGeometryType vg, int* unitNormal);

void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         BoundaryScheme boundaryScheme,
                         const std::vector<VariableTypes>& macroVarTypes,
                         const std::vector<Real>& macroVarValues,
                         const VertexType boundaryType = VertexType::Wall);

void DefineBlockBoundary(
    int blockIndex, int componentID, BoundarySurface boundarySurface,
    const VertexType boundaryType = VertexType::VirtualBoundary);
const std::vector<BlockBoundary>& BlockBoundaries();
#ifdef OPS_3D
void TreatBlockBoundary3D(const Block& block, const int componentID,
                          const Real* givenVars,
                          const BoundaryScheme boundaryScheme,
                          const BoundarySurface boundarySurface);
void ImplementBoundary3D();
#endif

#ifdef OPS_2D
void TreatBlockBoundary(const Block& block, const int componentID,
                          const Real* givenVars,
                          const BoundaryScheme boundaryScheme,
                          const BoundarySurface boundarySurface);
void ImplementBoundary();
#endif
#endif  // BOUNDARY_H
