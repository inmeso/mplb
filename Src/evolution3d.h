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
#include "hilemms.h"
#include "model.h"

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
// Stream-collision scheme related
/*!
 * Overall wrap for stream-collision scheme
 */
void StreamCollision3D(const Real time);
/*!
 * Ops_par_loop for the stream step
 */
void Stream3D();
/*!
 * Ops_par_loop for the collision step
 */
void PreDefinedCollision3D();
// finite-difference scheme with cutting-cell technique related
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
void PreDefinedInitialCondition3D();
/*!
 * Mainly for a steady simulation
 */
void DispResidualError3D(const int iter, const Real timePeriod);
void NormaliseF3D(Real* ratio);
void UpdateMacroVars3D();
void PreDefinedBodyForce3D();
void UpdateHalos3D();
void ImplementBoundary3D();
void CopyDistribution3D(ops_dat* fDest, const ops_dat* fSrc);

void TreatBlockBoundary3D(const int blockIndex, const int componentID,
                          const Real* givenVars, int* range,
                          const VertexTypes boundaryType);

void Iterate(const SizeType steps, const SizeType checkPointPeriod);
void Iterate(const Real convergenceCriteria, const SizeType checkPointPeriod);

void UpdateMacroscopicBodyForce(const Real time);
void SetInitialMacrosVars();

template <typename T>
void Iterate(void (*cycle)(T), const SizeType steps,
             const SizeType checkPointPeriod) {
    ops_printf("Starting the iteration...\n");
    for (int iter = 0; iter < steps; iter++) {
        const Real time{iter * TimeStep()};
        cycle(time);
        if ((iter % checkPointPeriod) == 0 && iter != 0) {
            UpdateMacroVars3D();
            CalcResidualError3D();
            DispResidualError3D(iter, checkPointPeriod * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
        }
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
    DestroyFlowfield();
}

template <typename T>
void Iterate(void (*cycle)(T), const Real convergenceCriteria,
             const SizeType checkPointPeriod) {
    int iter{0};
    Real residualError{1};
    do {
        const Real time{iter * TimeStep()};
        cycle(time);
        if ((iter % checkPointPeriod) == 0) {
            UpdateMacroVars3D();
            CalcResidualError3D();
            residualError = GetMaximumResidual(checkPointPeriod * TimeStep());
            DispResidualError3D(iter, checkPointPeriod * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
        }

        iter = iter + 1;
    } while (residualError >= convergenceCriteria);

    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
    DestroyFlowfield();
}

#endif /* OPS_3D */
#endif /* EVOLUTION3D_H_ */
