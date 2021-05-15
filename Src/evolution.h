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
//#include "boundary.h"
//#include "flowfield.h"
//#include "model.h"
//#include "scheme.h"
#include "type.h"
#include "field.h"
#include "flowfield.h"
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

// Stream-collision scheme related
/*!
 * Overall wrap for stream-collision scheme
 */
void StreamCollision(const Real time);

void Iterate(const SizeType steps, const SizeType checkPointPeriod,
             const SizeType start = 0);
void Iterate(const Real convergenceCriteria, const SizeType checkPointPeriod,
             const SizeType start = 0);

template <typename T>
void Iterate(void (*cycle)(T), const SizeType steps,
             const SizeType checkPointPeriod, const SizeType start = 0) {
    ops_printf("Starting the iteration...\n");
    for (SizeType iter = start; iter < start + steps; iter++) {
        const Real time{iter * TimeStep()};
        cycle(time);
        if (((iter + 1) % checkPointPeriod) == 0) {
            ops_printf("%d iterations!\n", iter + 1);
#ifdef OPS_3D
            UpdateMacroVars3D();
#endif
#ifdef OPS_2D
            UpdateMacroVars();
#endif
            WriteFlowfieldToHdf5((iter + 1));
            WriteDistributionsToHdf5((iter + 1));
            WriteNodePropertyToHdf5((iter + 1));
        }
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
}

template <typename T>
void Iterate(void (*cycle)(T), const Real convergenceCriteria,
             const SizeType checkPointPeriod, const SizeType start = 0) {
    SizeType iter{start};
    Real residualError{1};
    do {
        const Real time{iter * TimeStep()};
        cycle(time);
        iter = iter + 1;
        if ((iter % checkPointPeriod) == 0) {
#ifdef OPS_3D
            UpdateMacroVars3D();
#endif
#ifdef OPS_2D
            UpdateMacroVars();
#endif
            CalcResidualError();
            residualError = GetMaximumResidual(checkPointPeriod);
            DispResidualError(iter, checkPointPeriod);
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
        }
    } while (residualError >= convergenceCriteria);

    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
}

#endif /* EVOLUTION3D_H_ */
