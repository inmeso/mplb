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
#include "evolution3d.h"


/*
 * In the following routines, there are some variables are defined
 * for the convenience of the translator which may not be able to
 * understand a function parameter in the ops_par_loop call
 * Even though, a variable rather than a numerical literacy will need
 * some modifications in the Python translator.
 */
#ifdef OPS_3D

void Collision3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        // ops_par_loop(KerCollide3D, "KerCollide3D", g_Block[blockIndex],
        //              SPACEDIM, iterRng,
        //              ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
        //              ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
        //                          LOCALSTENCIL, "int", OPS_READ),
        //              ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
        //                          OPS_READ),
        //              ops_arg_dat(g_feq[blockIndex], NUMXI, LOCALSTENCIL,
        //                          "double", OPS_READ),
        //              ops_arg_dat(g_Tau[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
        //                          "double", OPS_READ),
        //              ops_arg_dat(g_Bodyforce[blockIndex], NUMXI, LOCALSTENCIL,
        //                          "double", OPS_READ),
        //              ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL,
        //                          "double", OPS_WRITE));
    }
}

void Stream3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerStream3D, "KerStream3D", g_Block[blockIndex], SPACEDIM,
                     iterRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_fStage[blockIndex], NUMXI,
                                 ONEPTLATTICESTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
                                 OPS_RW));
    }
}

void UpdateMacroVars3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerCalcMacroVars3D, "KerCalcMacroVars3D",
                     g_Block[blockIndex], SPACEDIM, iterRng,
                     ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
                                 OPS_READ),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW));
    }
}

void UpdateFeqandBodyforce3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        // time is not used in the current force
        Real* timeF{0};
        ops_par_loop(KerCalcBodyForce3D, "KerCalcBodyForce3D",
                     g_Block[blockIndex], SPACEDIM, iterRng,
                     ops_arg_gbl(&timeF, 1, "double", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_Bodyforce[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_RW));
    }
}
// TODO Delete unimplemented functionalities for release
void TreatBlockBoundary3D(const int blockIndex, const int componentID,
                          const Real* givenVars, int* range,
                          const VertexTypes boundaryType) {
    switch (boundaryType) {
        case Vertex_ExtrapolPressure1ST: {
            ops_par_loop(
                KerCutCellExtrapolPressure1ST3D,
                "KerCutCellExtrapolPressure1ST3D", g_Block[blockIndex],
                SPACEDIM, range,
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                            ONEPTREGULARSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_f[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
                            "double", OPS_RW));
        } break;
        case Vertex_EQMDiffuseRefl: {
            ops_par_loop(
                KerCutCellEQMDiffuseRefl3D, "KerCutCellEQMDiffuseRefl3D",
                g_Block[blockIndex], SPACEDIM, range,
                ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case Vertex_NoslipEQN: {
            ops_par_loop(
                KerCutCellNoslipEQN3D, "KerCutCellNoslipEQN3D",
                g_Block[blockIndex], SPACEDIM, range,
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case Vertex_Periodic: {
            ops_par_loop(KerCutCellPeriodic3D, "KerCutCellPeriodic3D",
                         g_Block[blockIndex], SPACEDIM, range,
                         ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL,
                                     "double", OPS_RW),
                         ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        default:
            break;
    }
    //}
}

// void TreatEmbeddedBoundary3D() {
//     for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
//         int* iterRng = BlockIterRng(blockIdx, IterRngBulk());
//         ops_par_loop(
//             KerCutCellImmersedBoundary3D, "KerCutCellImmersedBoundary3D",
//             g_Block[blockIdx], SPACEDIM, iterRng,
//             ops_arg_dat(g_NodeType[blockIdx], 1, LOCALSTENCIL, "int",
//             OPS_READ), ops_arg_dat(g_GeometryProperty[blockIdx], 1,
//             LOCALSTENCIL, "int",
//                         OPS_READ),
//             ops_arg_dat(g_f[blockIdx], NUMXI, LOCALSTENCIL, "double",
//             OPS_RW));
//     }
// }


void InitialiseSolution3D() {
    UpdateFeqandBodyforce3D();
    CopyDistribution3D(g_feq, g_f);
}

void CopyDistribution3D(const ops_dat* fSrc, ops_dat* fDest) {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerCopyf, "KerCopyf", g_Block[blockIndex], SPACEDIM,
                     iterRng,
                     ops_arg_dat(fSrc[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_dat(fDest[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_WRITE));
    }
}

void CalcTotalMass3D(double* totalMass) {
    ops_reduction massHandle =
        ops_decl_reduction_handle(sizeof(double), "double", "massHandle");
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
        ops_par_loop(KerCalcSumofDensity, "KerCalcSumofDensity",
                     g_Block[blockIdx], SPACEDIM, iterRng,
                     ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_reduce(massHandle, 1, "double", OPS_INC));
    }
    ops_reduction_result(massHandle, (double*)totalMass);
}
// TODO Delete for release
void NormaliseF3D(Real* ratio) {
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
        ops_par_loop(
            KerNormaliseF, "KerNormaliseF", g_Block[blockIdx], SPACEDIM,
            iterRng, ops_arg_gbl(ratio, 1, "double", OPS_READ),
            ops_arg_dat(g_f[blockIdx], NUMXI, LOCALSTENCIL, "double", OPS_RW));
    }
}
void CalcResidualError3D() {
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
            int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
            ops_par_loop(KerCalcMacroVarSquareofDifference,
                         "KerCalcMacroVarSquareofDifference", g_Block[blockIdx],
                         SPACEDIM, iterRng,
                         ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_MacroVarsCopy[blockIdx], NUMMACROVAR,
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_gbl(&macroVarIdx, 1, "int", OPS_READ),
                         ops_arg_reduce(g_ResidualErrorHandle[macroVarIdx], 1,
                                        "double", OPS_INC));
        }
    }
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        ops_reduction_result(g_ResidualErrorHandle[macroVarIdx],
                             (double*)&g_ResidualError[2 * macroVarIdx]);
    }
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
        ops_par_loop(KerCopyMacroVars, "KerCopyMacroVars3D", g_Block[blockIdx],
                     SPACEDIM, iterRng,
                     ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_MacroVarsCopy[blockIdx], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW));
    }
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
            int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
            ops_par_loop(KerCalcMacroVarSquare, "KerCalcMacroVarSquare3D",
                         g_Block[blockIdx], SPACEDIM, iterRng,
                         ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_gbl(&macroVarIdx, 1, "int", OPS_READ),
                         ops_arg_reduce(g_ResidualErrorHandle[macroVarIdx], 1,
                                        "double", OPS_INC));
        }
    }
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        ops_reduction_result(g_ResidualErrorHandle[macroVarIdx],
                             (double*)&g_ResidualError[2 * macroVarIdx + 1]);
        // ops_printf("\n macro id = %i, abs res error = %.28f, rel abs error =
        // %.28f", macroVarIdx, g_ResidualError[2 * macroVarIdx]*10E21,
        // g_ResidualError[2 * macroVarIdx+1]*10E21); ops_printf("\n macro id =
        // %i, displayed error = %.28f \n",macroVarIdx, g_ResidualError[2 *
        // macroVarIdx]/g_ResidualError[2 * macroVarIdx+1]);
    }
}

void DispResidualError3D(const int iter, const Real checkPeriod) {
    ops_printf("##########Residual Error at %i time step##########\n", iter);
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        Real residualError = g_ResidualError[2 * macroVarIdx] /
                             g_ResidualError[2 * macroVarIdx + 1] /
                             (checkPeriod * TimeStep());
        ops_printf("Residual of %s = %.17g\n",
                   MacroVarName()[macroVarIdx].c_str(), residualError);
    }
}

void StreamCollision3D() {
#if DebugLevel >= 1
    ops_printf("Calculating the macroscopic variables...\n");
#endif
    UpdateMacroVars3D();
    // Real TotalMass{0};
    // CalcTotalMass(&TotalMass);
    // Real Ratio{TotalMass/TotalMeshSize()};
    // NormaliseF(&Ratio);
    // UpdateMacroVars();
    CopyDistribution3D(g_f, g_fStage);
#if DebugLevel >= 1
    ops_printf("Calculating the equilibrium function and the body force term...\n");
#endif
    UpdateFeqandBodyforce3D();

#if DebugLevel >= 1
    ops_printf("Colliding...\n");
#endif
    Collision3D();
#if DebugLevel >= 1
    ops_printf("Streaming...\n");
#endif
    Stream3D();
#if DebugLevel >= 1
    ops_printf("Updating the halos...\n");
#endif
    if (nullptr != HaloGroup()) {
        ops_halo_transfer(HaloGroup());
    }
#if DebugLevel >= 1
    ops_printf("Implementing the boundary conditions...\n");
#endif
    ImplementBoundaryConditions();
}

#endif /* OPS_3D */
