// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief   Head file for wrap functions.
  * @author  Jianping Meng
  * @details Declare wrap functions for implementing the main evolution
  * cycle
  */

#ifndef EVOLUTION_H_
#define EVOLUTION_H_
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
// Routines for the stream-collision scheme.
void Stream();
void Collision();
void ImplementBoundary();
void CalcResidualError();
/*!
 * Routine for completing one full time step
 */
void StreamCollision();
void StreamCollisionSWE();
// Routines for the general finite-difference scheme.
void UpdateBoundary();
void TimeMarching();
// Common routines
void InitialiseSolution();
void UpdateMacroVars();
void UpdateTau();
void UpdateSWETau();
void UpdateFeqandBodyforce();
void UpdateSWEFeqandBodyforce();
void UpdateHalos();
void CopyDistribution(const ops_dat *fSrc, ops_dat *fDest);
void DispResidualError(const int iter, const Real timePeriod);
void NormaliseF(Real *ratio);
void TreatDomainBoundary(const Real *givenVars, int *range,
                         const VertexTypes boundaryType);
#endif /* EVOLUTION_H_ */
