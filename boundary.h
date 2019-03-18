// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

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
 *
 * @param vg Geometry property, e.g., corner type
 * @param discreteVelocity components of a discrete velocity
 * @return BndryDvType discrete velocity type
 */
BndryDvType FindBdyDvType3D(const VertexGeometryTypes vg,
                            const Real* discreteVelocity);
/*!
 * @brief Determining discrete velocity type a solid wall boundary
 *
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
 *
 * @param givenMacroVars  specified velocity
 * @param nodeType if the current node is set to be EDR node
 * @param geometryProperty e.g., corner types
 * @param f distribution function
 */
void KerCutCellEQMDiffuseRefl(const Real* givenMacroVars, const int* nodeType,
                              const int* geometryProperty, Real* f,
                              const int* componentId);
void KerCutCellPeriodic(const int* nodeType, const int* geometryProperty,
                        Real* f);
void KerCutCellBounceBack(const int* nodeType, const int* geometryProperty,
                          Real* f);
void KerCutCellBounceBackNew(const int* nodeType, const int* geometryProperty,
                             Real* f);
void KerCutCellKinetic(const Real* givenMacroVars, const int* nodeType,
                       const int* geometryProperty, Real* f);
void KerCutCellCorrectedKinetic(const Real* givenMacroVars, const Real* dt,
                                const int* nodeType,
                                const int* geometryProperty, const Real* tau,
                                const Real* feq, Real* f);
void KerCutCellExtrapolPressure1ST(const Real* givenBoundaryVars,
                                   const int* nodeType,
                                   const int* geometryProperty, Real* f);
void KerCutCellExtrapolPressure2ND(const Real* givenBoundaryVars,
                                   const int* nodeType,
                                   const int* geometryProperty, Real* f);
void KerCutCellNonEqExtrapolPressure(const Real* givenMacroVars,
                                     const int* nodeType,
                                     const int* geometryProperty,
                                     const Real* macroVars, const Real* feq,
                                     Real* f);

/*!
 * The Zou-He boundary condition can only be valid for the D2Q9 lattice and
 * the second order equilibrium function
 */
void KerCutCellZouHeVelocity(const Real* givenMacroVars, const int* nodeType,
                             const int* geometryProperty, const Real* macroVars,
                             Real* f);
void KerCutCellZeroFlux(const int* nodeType, const int* geometryProperty,
                        Real* f);
/*
 * @givenMacroVars:given macroscopic variables for the boundary condition.
 */
void KerEquibriumVelocity(const Real* givenMacroVars, Real* f);
// CutCell immersed solid boundary condition
void KerCutCellEmbeddedBoundary(const int* nodeType, const int* geometryProperty,
                               Real* f);
/*!
 * Non-equilibrium Extraploation boundary condition by Guo et al.
 */
void KerCutCellNonEqExtrapol(const Real* givenMacroVars, const int* nodeType,
                             const int* geometryProperty, const Real* macroVars,
                             const Real* feq, Real* f);

// Boundary fitting
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
                                     const int* nodeType,
                                     const int* geometryProperty, Real* f);
/*!
 * @brief  Equilibrium diffuse reflection boundary condition: 3D
 * @param givenMacroVars  specified velocity
 * @param nodeType if the current node is set to be EDR node
 * @param geometryProperty e.g., corner types
 * @param f distribution function
 */
void KerCutCellEQMDiffuseRefl3D(const Real* givenMacroVars, const int* nodeType,
                                const int* geometryProperty, Real* f,
                                const int* componentId);

/*!
 * @brief  The 3D no-slip boundary condition using EQN scheme
 * @param givenMacroVars  Values of macroscopic variables
 * @param nodeType if the current node is a Dirichlet node
 * @param f distribution function
 */
void KerCutCellNoslipEQN3D(const Real* givenMacroVars, const int* nodeType,
                           Real* f, const int* componentId);

void KerCutCellPeriodic3D(const int* nodeType, const int* geometryProperty,
                          Real* f, const int* componentId);
#endif /* OPS_3D*/
/*!
 * For updating halo points, f must have wind direction
 * the geometry variable can use central difference
 * According the boundary condition, may need different method
 */
void KerUpdateHalo();
const int BoundaryHaloNum();
void SetBoundaryHaloNum(const int boundaryhaloNum);
void SetupBoundary();
#endif  // BOUNDARY_H
