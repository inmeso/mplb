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

/*!
 * @brief   Wrap functions for main evolution cycle.
 * @author  Jianping Meng
 * @details Define wrap functions for implementing the main evolution
 * cycle
 */
#include "evolution.h"
#include "type.h"
#include "scheme.h"
#include "block.h"
#include "field.h"
#include "boundary.h"
#include "flowfield.h"
#include "model.h"

/*
 * In the following routines, there are some variables are defined
 * for the convenience of the translator which may not be able to
 * understand a function parameter in the ops_par_loop call
 * Even though, a variable rather than a numerical literacy will need
 * some modifications in the Python translator.
 */
// void TreatEmbeddedBoundary3D() {
//     for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
//         int* iterRng.data() = BlockIterRng(blockIdx, IterRngBulk());
//         ops_par_loop(
//             KerCutCellImmersedBoundary3D, "KerCutCellImmersedBoundary3D",
//             g_Block[blockIdx], SpaceDim(), iterRng.data(),
//             ops_arg_dat(g_NodeType()[blockIdx], 1, LOCALSTENCIL, "int",
//             OPS_READ), ops_arg_dat(g_GeometryProperty[blockIdx], 1,
//             LOCALSTENCIL, "int",
//                         OPS_READ),
//             ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL, "double",
//             OPS_RW));
//     }
// }


void Iterate(const SizeType steps, const SizeType checkPointPeriod,
             const SizeType start) {
    const SchemeType scheme = Scheme();
    ops_printf("Starting the iteration...\n");
    switch (scheme) {
        case Scheme_StreamCollision: {
            for (SizeType iter = start; iter < start + steps; iter++) {
                const Real time{iter * TimeStep()};
                StreamCollision(time);
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
        } break;
        default:
            break;
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();

}

void Iterate(const Real convergenceCriteria, const SizeType checkPointPeriod,
             const SizeType start) {
    const SchemeType scheme = Scheme();
    ops_printf("Starting the iteration...\n");
    switch (scheme) {
        case Scheme_StreamCollision: {
            SizeType iter{start};
            Real residualError{1};
            do {
                const Real time{iter * TimeStep()};
                StreamCollision(time);
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
        } break;
        default:
            break;
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
}

void StreamCollision(const Real time) {
#if DebugLevel >= 1
    ops_printf("Calculating the macroscopic variables...\n");
#endif
#ifdef OPS_3D
    UpdateMacroVars3D();
#endif
#ifdef OPS_2D
    UpdateMacroVars();
#endif
    CopyBlockEnvelopDistribution(g_fStage(), g_f());
#if DebugLevel >= 1
    ops_printf("Calculating the mesoscopic body force term...\n");
#endif
    UpdateMacroscopicBodyForce(time);
#ifdef OPS_3D
    PreDefinedBodyForce3D();
#endif
#ifdef OPS_2D
    PreDefinedBodyForce();
#endif
#if DebugLevel >= 1
    ops_printf("Calculating the collision term...\n");
#endif
#ifdef OPS_3D
    PreDefinedCollision3D();
#endif
#ifdef OPS_2D
    PreDefinedCollision();
#endif

#if DebugLevel >= 1
    ops_printf("Updating the halos...\n");
#endif
    TransferHalos();

#if DebugLevel >= 1
    ops_printf("Streaming...\n");
#endif
#ifdef OPS_3D
    Stream3D();
#endif
#ifdef OPS_2D
    Stream();
#endif

#if DebugLevel >= 1
    ops_printf("Implementing the boundary conditions...\n");
#endif

#ifdef OPS_3D
    ImplementBoundary3D();
#endif
#ifdef OPS_2D
    ImplementBoundary();
#endif
}

