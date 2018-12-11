// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

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
/*!
 * g_HaloDepth the number of halo points
 * Determined by the chosen scheme and boundary condition
 * Periodic boundary condition requires one halo point
 * but it is mainly used for transferring information
 */

// Declaring kernel functions which to be called by ops_par_loop

// Finite difference schemes for boundary fitting mesh
/*!
 * General function to calculate the gradient
 * len is the dimension of var, the dimension of grad must be 2*len for 2D case
 * the halo points must be filled beforehand according to the chosen scheme.
 * To be improved in the future
 */
void KerGradCentral2nd(const Real* var, Real* grad, const int len);
void KerGradCentral4th(const Real* var, Real* grad, const int len);
void KerGradCentral6th(const Real* var, Real* grad, const int len);
void KerCVTUpwind2nd(const Real* metrics, const Real* f, Real* xidotgrad,
                     const Real* XI);
void KerCVTUpwind4th(const Real* metrics, const Real* f, Real* xidotgrad,
                     const Real* XI);
void KerCVTUpwind6th(const Real* metrics, const Real* f, Real* xidotgrad,
                     const Real* XI);
// End finite difference schemes for boundary fitting mesh
// Schemes for cutting cell mesh technique
// stream-collision scheme
#ifdef OPS_2D
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
#endif
#ifdef OPS_2D
/*!
 * @fn KerStream
 * @brief Stream step for the stream-collision scheme
 * @param nodeType node type
 * @param geometry geometry type, e.g., if it is a corner
 * @param fStage temporary storage
 * @param f distribution function
 */
void KerStream(const int* nodeType, const int* geometry, const Real* fStage,
               Real* f);
#endif
#ifdef OPS_3D
/*!
 * @fn KerStream3D
 * @brief Stream step for the stream-collision scheme: 3D case
 * @param nodeType node type
 * @param geometry geometry type, e.g., if it is a corner
 * @param fStage temporary storage
 * @param f distribution function
 */
void KerStream3D(const int* nodeType, const int* geometry, const Real* fStage,
                 Real* f);
#endif
#ifdef OPS_2D
// Finite difference scheme for the cutting cell mesh
/*!
 * First order upwind scheme for term c dot df/dx
 */
void KerCutCellCVTUpwind1st(const Real* coordinateXYZ, const int* nodeType,
                            const int* geometry, const Real* f,
                            Real* fGradient);
/*!
 * Second order upwind scheme for term c dot df/dx
 */
void KerCutCellCVTUpwind2nd(const Real* coordinateXYZ, const int* nodeType,
                            const int* geometry, const Real* f,
                            Real* fGradient);
/*!
 * A semi-implicit time scheme where f in the n+1 time step is used in
 * collision term
 */
void KerCutCellSemiImplicitTimeMach(const Real* dt, const Real* schemeCoeff,
                                    const int* nodeType, const int* geometry,
                                    const Real* fGradient, const Real* feq,
                                    const Real* relaxationTime,
                                    const Real* bodyForce, Real* f);
/*!
 * A general explicit time scheme which can be defined by using schemeCoeff
 */
void KerCutCellExplicitTimeMach(const Real* dt, const Real* schemeCoeff,
                                const int* nodeType, const int* geometry,
                                const Real* fGradient, const Real* feq,
                                const Real* relaxationTime,
                                const Real* bodyForce, Real* f);
// End Finite difference scheme for the cutting cell mesh
#endif

/*!
 * Utility kernel function for setting distribution to a fixed value
 */
void KerSetfFixValue(const Real* value, Real* f);
/*!
 * Utility kernel function for re-normalise distribution function by a ratio
 */
void KerNormaliseF(const Real* ratio, Real* f);
/*!
 * Utility kernel function for copying distribution function
 */
void KerCopyf(const Real* src, Real* dest);
/*!
 * Utility kernel function for calculating numerator of L2 norm
 */
void KerCalcMacroVarSquareofDifference(const Real* macroVars,
                                       const Real* macroVarsCopy,
                                       const int* varId, double* sumSquareDiff);

/*!
 * Utility kernel function for calculating denominator of L2 norm
 * Mainly for a steady problem
 */
void KerCalcMacroVarSquare(const Real* macroVars, const int* varId,
                           double* sumSquare);
/*!
 * Utility kernel function for calculating whole block sum of density
 */
void KerCalcSumofDensity(const Real* macroVars, double* densitySum);
/*!
 * Utility kernel function for copying geometry and node property data
 */
void KerCopyProperty(const int* src, int* dest);
/*!
 * Utility kernel function for copying macroscopic variables
 */
void KerCopyMacroVars(const Real* src, Real* dest);
/*!
 * Utility kernel function for copying distribution with a displacement
 */
void KerCopyDispf(const Real* src, Real* dest, const int* disp);
/*!
 * Utility kernel function for copying coordinates
 */
void KerCopyCoordinateXYZ(const Real* src, Real* dest);
/*!
 * Set a scalar variable to a specific value.Mainly used in initialisation.
 */
void KerAssignProperty(const int* value, int* var);
/*!
 * Set macroscopic variable to values specified by value.
 */
void KerSetMacroVarToConst(const Real* value, Real* macroVar);
/*!
 * Kernel function for getting the value at a grid point
 */
void KerGetPointMacroVarValue(const Real* macroVars, Real* pointValue);

void SetupScheme();
/*!
 * Periodic boundary condition may need to adjust the g_HaloDepth value
 */
const int SchemeHaloNum();
void SetSchemeHaloNum(const int schemeHaloNum);
#endif
