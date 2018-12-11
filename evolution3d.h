// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief   Head file for wrap functions.
  * @author  Jianping Meng
  * @details Declare wrap functions for implementing the main evolution
  * cycle
  */

#ifndef EVOLUTION3D_H_
#define EVOLUTION3D_H_
#include "boundary.h"
#include "flowfield.h"
#include "scheme.h"
#include "type.h"

/*!
 * In this module, we provide subroutines that implement a update system over
 * the whole block and block by block if multi-block
 * The subroutines belongs to three categories:
 * 1. updating boundary with the chosen boundary conditions
 * 2. updating bulk with the chose spatial scheme while the time-marching will
 *    be implemented in this module.
 * 3. updating the macroscopic variables, equilibrium function, body force term
 *    and other parameters as required.
 * Requirement:
 * 1. The ops_par_loop will be called here
 * 2. The subroutines should not directly alter any variables defined in other
 *    modules.  Instead, the kernel functions should be used. This mechanism
 *    allows the subroutines to directly read variables such as g_f
 */

#ifdef OPS_3D
//Stream-collision scheme related
/*!
 * Overall wrap for stream-collision scheme
 */
void StreamCollision3D();
/*!
 * Ops_par_loop for the stream step
 */
void Stream3D();
/*!
 * Ops_par_loop for the collision step
 */
void Collision3D();
//finite-difference scheme with cutting-cell technique related
/*!
 * Overall wrap for a finite difference scheme
 */
void TimeMarching3D();
void ForwardEuler3D();
// Common routines
/*!
 * Calculating the residual errors for each macroscopic variables
 */
void CalcResidualError3D();
/*!
 * Calculating the total mass of the whole problem domain (all the blocks)
 */
void CalcTotalMass3D(Real* totalMass);
void InitialiseSolution3D();
/*!
 * Mainly for a steady simulation
 */
void DispResidualError3D(const int iter, const Real timePeriod);
void NormaliseF3D(Real* ratio);
void UpdateMacroVars3D();
void UpdateTau3D();
void UpdateFeqandBodyforce3D();
void UpdateHalos3D();
void ImplementBoundary3D();
void CopyDistribution3D(const ops_dat* fSrc, ops_dat* fDest);
#endif /* OPS_3D */
#endif /* EVOLUTION3D_H_ */
