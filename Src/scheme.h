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

// Declaring kernel functions which to be called by ops_par_loop
#ifdef OPS_2D
// Schemes for cutting cell mesh technique
// stream-collision scheme
/*!
 * @fn Collision step for the stream-collision scheme
 * @param dt time step
 * @param nodeType the node
 * @param f distribution function
 * @param feq equilibrium function
 * @param relaxationTime relaxation time
 * @param bodyForce force term
 * @param fStage temporary storage
 */
void KerCollide(const Real* dt, const int* nodeType, const Real* f,
                const Real* feq, const Real* relaxationTime,
                const Real* bodyForce, Real* fStage);

/*!
 * @fn KerStream
 * @brief Stream step for the stream-collision scheme
 * @param nodeType node type
 * @param geometry geometry type, e.g., if it is a corner
 * @param fStage temporary storage
 * @param f distribution function
 */
void KerStream(ACC<Real>& f, const ACC<Real>& fStage, const ACC<int>& nodeType,
               const ACC<int>& geometry);
#endif

#ifdef OPS_3D
/*!
 * @fn KerCollide3D
 * @brief Collision step for the stream-collision scheme: 3D case
 * @param dt time step
 * @param nodeType the node
 * @param f distribution function
 * @param feq equilibrium function
 * @param relaxationTime relaxation time
 * @param bodyForce force term
 * @param fStage temporary storage
 */
void KerCollide3D(const Real* dt, const int* nodeType, const Real* f,
                  const Real* feq, const Real* relaxationTime,
                  const Real* bodyForce, Real* fStage);
/*!
 * @fn KerStream3D
 * @brief Stream step for the stream-collision scheme: 3D case
 * @param nodeType node type
 * @param geometry geometry type, e.g., if it is a corner
 * @param fStage temporary storage
 * @param f distribution function
 */
void KerStream3D(ACC<Real>& f, const ACC<Real>& fStage,
                 const ACC<int>& nodeType, const ACC<int>& geometry);
#endif
#ifdef OPS_2D
// Finite difference scheme for the cutting cell mesh
/*!
 * First order upwind scheme for term c dot df/dx
 */
void KerCutCellCVTUpwind1st(const ACC<Real>& coordinateXYZ,
                            const ACC<int>& nodeType, const ACC<int>& geometry,
                            const ACC<Real>& f, ACC<Real>& fGradient)
/*!
 * Second order upwind scheme for term c dot df/dx
 */
KerCutCellCVTUpwind2nd(const ACC<Real>& coordinateXYZ, const ACC<int> nodeType,
                       const ACC<int>& geometry, const ACC<Real>& f,
                       ACC<Real>& fGradient);
/*!
 * A semi-implicit time scheme where f in the n+1 time step is used in
 * collision term
 */
void KerCutCellSemiImplicitTimeMach(
    const Real* dt, const Real* schemeCoeff, const ACC<int>& nodeType,
    const ACC<int>& geometry, const ACC<Real>& fGradient, const ACC<Real>& feq,
    const ACC<Real>& relaxationTime, const ACC<Real>& bodyForce, Real* f);
/*!
 * A general explicit time scheme which can be defined by using schemeCoeff
 */
void KerCutCellExplicitTimeMach(
    const Real* dt, const Real* schemeCoeff, const ACC<int>& nodeType,
    const ACC<int>& geometry, const ACC<Real>& fGradient, const ACC<Real>& feq,
    const ACC<Real>& relaxationTime, const ACC<Real>& bodyForce, ACC<Real>& f);
// End Finite difference scheme for the cutting cell mesh
#endif
//HiLeMMS interface see https://gitlab.com/jpmeng/hilemms
void DefineScheme(const SchemeType scheme);
/*!
 * Utility kernel function for setting distribution to a fixed value
 */
void KerSetfFixValue(const Real* value, ACC<Real>& f);
/*!
 * Utility kernel function for re-normalise distribution function by a ratio
 */
void KerNormaliseF(const Real* ratio, ACC<Real>& f);
/*!
 * Utility kernel function for copying distribution function
 */
void KerCopyf(ACC<Real>& dest, const ACC<Real>& src);
/*!
 * Utility kernel function for calculating numerator of L2 norm
 */
void KerCalcMacroVarSquareofDifference(const ACC<Real>& macroVars,
                                       const ACC<Real>& macroVarsCopy,
                                       const int* varId, double* sumSquareDiff);

/*!
 * Utility kernel function for calculating denominator of L2 norm
 * Mainly for a steady problem
 */
void KerCalcMacroVarSquare(const ACC<Real>& macroVars, const int* varId,
                           double* sumSquare);
/*!
 * Utility kernel function for calculating whole block sum of density
 */
void KerCalcSumofDensity(const ACC<Real>& macroVars, double* densitySum);
/*!
 * Utility kernel function for copying geometry and node property data
 */
void KerCopyProperty(const ACC<int>& src, ACC<int>& dest);
/*!
 * Utility kernel function for copying macroscopic variables
 */
void KerCopyMacroVars(const ACC<Real> &src, ACC<Real>& dest);
/*!
 * Utility kernel function for copying distribution with a displacement
 */
void KerCopyDispf(const ACC<Real>& src, ACC<Real>& dest, const int* disp);
/*!
 * Utility kernel function for copying coordinates
 */
void KerCopyCoordinateXYZ(const ACC<Real>& src, ACC<Real>& dest);
/*!
 * Set a scalar variable to a specific value.Mainly used in initialisation.
 */
void KerSetGeometryProperty(const int* value, ACC<int>& var);
void KerSetNodeType(const int* value, ACC<int>& var, const int* compoId);
/*!
 * Set macroscopic variable to values specified by value.
 */
void KerSetMacroVarToConst(const Real* value, ACC<Real>& macroVar);
/*!
 * Kernel function for getting the value at a grid point
 */
void KerGetPointMacroVarValue(const ACC<Real>& macroVars, Real* pointValue);
const int SchemeHaloNum();
void SetSchemeHaloNum(const int schemeHaloNum);
const SchemeType Scheme();
#endif
