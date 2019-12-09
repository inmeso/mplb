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
#include "flowfield.h"
#include "model.h"
/*!
 * Discrete velocity type at a solid wall boundary
 */
enum BndryDvType {
    BndryDv_Incoming = 1,
    BndryDv_Outgoing = 2,
    BndryDv_Parallel = 3,
    BndryDv_Invalid = -1
};
/*!
 * @brief Determining discrete velocity type a solid wall boundary: 3D
 * @param vg Geometry property, e.g., corner type
 * @param discreteVelocity components of a discrete velocity
 * @return BndryDvType discrete velocity type
 */
BndryDvType FindBdyDvType3D(const VertexGeometryTypes vg,
                            const Real* discreteVelocity);
/*!
 * @brief Determining discrete velocity type a solid wall boundary
 * @param vg Geometry property, e.g., corner type
 * @param discreteVelocity components of a discrete velocity
 * @return BndryDvType discrete velocity type
 */
BndryDvType FindBdyDvType(const VertexGeometryTypes vg,
                          const Real* discreteVelocity);
#ifdef OPS_2D
// CutCell block boundary condition
/*!
 * @brief  Equilibrium diffuse reflection boundary condition
 * @param givenMacroVars  specified velocity
 * @param nodeType if the current node is set to be EDR node
 * @param geometryProperty e.g., corner types
 * @param f distribution function
 * see Meng, Gu Emerson, Peng and Zhang, https://arxiv.org/abs/1803.00390.
 */
void KerCutCellEQMDiffuseRefl(const Real* givenMacroVars,
                              const ACC<int>& nodeType,
                              const ACC<int>& geometryProperty, ACC<Real>& f,
                              const int* componentId);
void KerCutCellPeriodic(const ACC<int>& nodeType,
                        const ACC<int>& geometryProperty, ACC<Real>& f);
/*!
 * Standard bounce back boundary condition for block boundaries
 * As discussed by Meng, Gu and Emerson, Journal of Computational Science
 * 2018 (28): 476-482, its definition is incomplete which will induc
 * non-physical slip velocity at wall
 */
void KerCutCellBounceBack(const ACC<int>& nodeType,
                          const ACC<int>& geometryProperty, ACC<Real>& f);

void KerCutCellKinetic(const Real* givenMacroVars, const int* nodeType,
                       const int* geometryProperty, Real* f);
void KerCutCellCorrectedKinetic(const Real* givenMacroVars, const Real* dt,
                                const ACC<int>& nodeType,
                                const ACC<int>& geometryProperty,
                                const ACC<Real>& tau, const ACC<Real>& feq,
                                ACC<Real>& f);
void KerCutCellExtrapolPressure1ST(const Real* givenBoundaryVars,
                                   const ACC<int>& nodeType,
                                   const ACC<int>& geometryProperty,
                                   ACC<Real>& f);
void KerCutCellExtrapolPressure2ND(const Real* givenBoundaryVars,
                                   const ACC<int>& nodeType,
                                   const ACC<int>& geometryProperty,
                                   ACC<Real>& f);


/*!
 * The Zou-He boundary condition can only be valid for the D2Q9 lattice and
 * the second order equilibrium function
 */
void KerCutCellZouHeVelocity(const Real* givenMacroVars,
                             const ACC<int>& nodeType,
                             const ACC<int>& geometryProperty,
                             const ACC<Real>& macroVars, ACC<Real>& f);
void KerCutCellZeroFlux(const ACC<int>& nodeType,
                        const ACC<int>& geometryProperty, ACC<Real>& f);
/*
 * @givenMacroVars:given macroscopic variables for the boundary condition.
 */

// CutCell immersed solid boundary condition
void KerCutCellEmbeddedBoundary(const ACC<int>& nodeType,
                                const ACC<int>& geometryProperty, ACC<Real>& f);

#endif /* OPS_2D  */
#ifdef OPS_3D
// CutCell block boundary condition
/*!
 * @brief First order extrapolation pressure flow boundary:3D
 * @param givenBoundaryVars specified pressure
 * @param nodeType if the current node is set to be pressure flow boundary node
 * @param geometryProperty e.g., corner types
 * @param f distribution
 */
void KerCutCellExtrapolPressure1ST3D(const Real* givenBoundaryVars,
                                     const ACC<int>& nodeType,
                                     const ACC<int>& geometryProperty,
                                     ACC<Real>& f);
/*!
 * @brief  Equilibrium diffuse reflection boundary condition: 3D
 * @param givenMacroVars  specified velocity
 * @param nodeType if the current node is set to be EDR node
 * @param geometryProperty e.g., corner types
 * @param f distribution function
 */
void KerCutCellEQMDiffuseRefl3D(ACC<Real>& f, const ACC<int>& nodeType,
                                const ACC<int>& geometryProperty,
                                const Real* givenMacroVars,
                                const int* componentId);

/*!
 * @brief  The 3D no-slip boundary condition using EQN scheme
 * @param givenMacroVars  Values of macroscopic variables
 * @param nodeType if the current node is a Dirichlet node
 * @param f distribution function
 */
void KerCutCellNoslipEQN3D(ACC<Real>& f, const ACC<int>& nodeType,
                           const Real* givenMacroVars, const int* componentId);

void KerCutCellPeriodic3D(ACC<Real>& f, const ACC<int>& nodeType,
                          const ACC<int>& geometryProperty,
                          const int* componentId);
#endif /* OPS_3D*/

const int BoundaryHaloNum();
void SetBoundaryHaloNum(const int boundaryhaloNum);
void BoundaryNormal3D(const VertexGeometryTypes vg, int* unitNormal);
#endif  // BOUNDARY_H
