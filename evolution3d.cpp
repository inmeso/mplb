// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*!
 * @brief   Wrap functions for main evolution cycle.
 * @author  Jianping Meng
 * @details Define wrap functions for implementing the main evolution
 * cycle
 */
#include "evolution3d.h"
#include "hilemms.h"
#include "model.h"

/*
 * In the following routines, there are some variables are defined
 * for the convenience of the translator which may not be able to
 * understand a function parameter in the ops_par_loop call
 * Even though, a variable rather than a numerical literacy will need
 * some modifications in the Python translator.
 */
#ifdef OPS_3D
void UpdateTau3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerCalcTau3D, "KerCalcTau3D", g_Block[blockIndex],
                     SPACEDIM, iterRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_gbl(TauRef(), NUMCOMPONENTS, "double", OPS_READ),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_Tau[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
                                 "double", OPS_RW));
    }
}

void Collision3D() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerCollide3D, "KerCollide3D", g_Block[blockIndex],
                     SPACEDIM, iterRng,
                     ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double",
                                 OPS_READ),
                     ops_arg_dat(g_feq[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_dat(g_Tau[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_dat(g_Bodyforce[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_WRITE));
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
        ops_par_loop(KerCalcFeq3D, "KerCalcFeq3D", g_Block[blockIndex],
                     SPACEDIM, iterRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_feq[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_RW));

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

void TreatBlockBoundary3D(const int blockIndex, const int componentID,
                          const Real* givenVars, int* range,
                          const VertexTypes boundaryType) {
    // void TreatBlockBoundary3D(const Real* givenVars, int* range, const
    // VertexTypes boundaryType) for (int blockIdx = 0; blockIdx < BlockNum();
    // blockIdx++) {
    int blockIdx;
    // This way will require minimum changes in the MPLB code.
    blockIdx = blockIndex;
    switch (boundaryType) {
        case Vertex_NoneqExtrapol: {
            //                ops_par_loop(
            //                    KerCutCellNonEqExtrapol,
            //                    "KerCutCellNonEqExtrapol", g_Block[blockIdx],
            //                    SPACEDIM, range, ops_arg_gbl(givenVars,
            //                    NUMMACROVAR, "double", OPS_READ),
            //                    ops_arg_dat(g_NodeType[blockIdx], 1,
            //                    LOCALSTENCIL, "int",
            //                                OPS_READ),
            //                    ops_arg_dat(g_GeometryProperty[blockIdx], 1,
            //                    LOCALSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_MacroVars[blockIdx],
            //                    NUMMACROVAR,
            //                                ONEPTLATTICESTENCIL, "double",
            //                                OPS_READ),
            //                    ops_arg_dat(g_feq[blockIdx], NUMXI,
            //                    ONEPTLATTICESTENCIL,
            //                                "double", OPS_RW),
            //                    ops_arg_dat(g_f[blockIdx], NUMXI,
            //                    ONEPTLATTICESTENCIL,
            //                                "double", OPS_RW));
        } break;
        case Vertex_ExtrapolPressure1ST: {
            ops_par_loop(
                KerCutCellExtrapolPressure1ST3D,
                "KerCutCellExtrapolPressure1ST3D", g_Block[blockIdx], SPACEDIM,
                range, ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_dat(g_NodeType[blockIdx], NUMCOMPONENTS,
                            ONEPTREGULARSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty[blockIdx], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_f[blockIdx], NUMXI, ONEPTREGULARSTENCIL, "double",
                            OPS_RW));
        } break;
        case Vertex_ExtrapolPressure2ND: {
            //                ops_par_loop(
            //                    KerCutCellExtrapolPressure2ND,
            //                    "KerCutCellExtrapolPressure2ND",
            //                    g_Block[blockIdx], SPACEDIM, range,
            //                    ops_arg_gbl(givenVars, NUMMACROVAR, "double",
            //                    OPS_READ), ops_arg_dat(g_NodeType[blockIdx],
            //                    1, ONEPTREGULARSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_GeometryProperty[blockIdx], 1,
            //                    LOCALSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_f[blockIdx], NUMXI,
            //                    TWOPTREGULARSTENCIL,
            //                                "double", OPS_RW));
        } break;
        case Vertex_NonEqExtrapolPressure: {
            //                ops_par_loop(
            //                    KerCutCellNonEqExtrapolPressure,
            //                    "KerCutCellNonEqExtrapolPressure",
            //                    g_Block[blockIdx], SPACEDIM, range,
            //                    ops_arg_gbl(givenVars, NUMMACROVAR, "double",
            //                    OPS_READ), ops_arg_dat(g_NodeType[blockIdx],
            //                    1, LOCALSTENCIL, "int",
            //                                OPS_READ),
            //                    ops_arg_dat(g_GeometryProperty[blockIdx], 1,
            //                    LOCALSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_MacroVars[blockIdx],
            //                    NUMMACROVAR,
            //                                ONEPTLATTICESTENCIL, "double",
            //                                OPS_READ),
            //                    ops_arg_dat(g_feq[blockIdx], NUMXI,
            //                    ONEPTLATTICESTENCIL,
            //                                "double", OPS_RW),
            //                    ops_arg_dat(g_f[blockIdx], NUMXI,
            //                    ONEPTLATTICESTENCIL,
            //                                "double", OPS_RW));
        } break;
        case Vertex_ZouHeVelocity: {
            //                ops_par_loop(
            //                    KerCutCellZouHeVelocity,
            //                    "KerCutCellZouHeVelocity,", g_Block[blockIdx],
            //                    SPACEDIM, range, ops_arg_gbl(givenVars,
            //                    NUMMACROVAR, "double", OPS_READ),
            //                    ops_arg_dat(g_NodeType[blockIdx], 1,
            //                    LOCALSTENCIL, "int",
            //                                OPS_READ),
            //                    ops_arg_dat(g_GeometryProperty[blockIdx], 1,
            //                    LOCALSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_MacroVars[blockIdx],
            //                    NUMMACROVAR,
            //                                ONEPTLATTICESTENCIL, "double",
            //                                OPS_READ),
            //                    ops_arg_dat(g_f[blockIdx], NUMXI,
            //                    ONEPTLATTICESTENCIL,
            //                                "double", OPS_RW));
        } break;
        case Vertex_KineticDiffuseWall: {
            //                ops_par_loop(
            //                    KerCutCellKinetic, "KerCutCellKinetic",
            //                    g_Block[blockIdx], SPACEDIM, range,
            //                    ops_arg_gbl(givenVars, NUMMACROVAR, "double",
            //                    OPS_READ), ops_arg_dat(g_NodeType[blockIdx],
            //                    1, LOCALSTENCIL, "int",
            //                                OPS_READ),
            //                    ops_arg_dat(g_GeometryProperty[blockIdx], 1,
            //                    LOCALSTENCIL,
            //                                "int", OPS_READ),
            //                    ops_arg_dat(g_f[blockIdx], NUMXI,
            //                    LOCALSTENCIL, "double",
            //                                OPS_RW));
        } break;
        case Vertex_EQMDiffuseRefl: {
            ops_par_loop(
                KerCutCellEQMDiffuseRefl3D, "KerCutCellEQMDiffuseRefl3D",
                g_Block[blockIdx], SPACEDIM, range,
                ops_arg_dat(g_f[blockIdx], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType[blockIdx], NUMCOMPONENTS, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty[blockIdx], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case Vertex_NoslipEQN: {
            ops_par_loop(
                KerCutCellNoslipEQN3D, "KerCutCellNoslipEQN3D",
                g_Block[blockIdx], SPACEDIM, range,
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_dat(g_NodeType[blockIdx], NUMCOMPONENTS, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_dat(g_f[blockIdx], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case Vertex_FreeFlux: {
            //                ops_par_loop(KerCutCellZeroFlux,
            //                "KerCutCellZeroFlux",
            //                             g_Block[blockIdx], SPACEDIM, range,
            //                             ops_arg_dat(g_NodeType[blockIdx], 1,
            //                             LOCALSTENCIL,
            //                                         "int", OPS_READ),
            //                             ops_arg_dat(g_GeometryProperty[blockIdx],
            //                             1,
            //                                         LOCALSTENCIL, "int",
            //                                         OPS_READ),
            //                             ops_arg_dat(g_f[blockIdx], NUMXI,
            //                             LOCALSTENCIL,
            //                                         "double", OPS_RW));
        } break;
        case Vertex_Periodic: {
            ops_par_loop(KerCutCellPeriodic3D, "KerCutCellPeriodic3D",
                         g_Block[blockIdx], SPACEDIM, range,
                         ops_arg_dat(g_f[blockIdx], NUMXI, LOCALSTENCIL,
                                     "double", OPS_RW),
                         ops_arg_dat(g_NodeType[blockIdx], NUMCOMPONENTS,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty[blockIdx], 1,
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

/*
void ImplementBoundary3D() {
    // TreatEmbeddedBoundary3D();
    // Real givenInletVars[]{1.00005, 0, 0};
    int* inletRng = BlockIterRng(0, IterRngImin());
    Real givenInletVars[]{1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenInletVars, inletRng, Vertex_Periodic);

    int* outletRng = BlockIterRng(0, IterRngImax());
    Real givenOutletVars[] = {1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenOutletVars, outletRng, Vertex_Periodic );

    int* topRng = BlockIterRng(0, IterRngJmax());
    // Real givenTopWallBoundaryVars[]{1, 0, 0};
    Real givenTopWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenTopWallBoundaryVars, topRng,
                          Vertex_NoslipEQN);

    int* bottomRng = BlockIterRng(0, IterRngJmin());
    Real givenBotWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenBotWallBoundaryVars, bottomRng,
                         Vertex_NoslipEQN);

    int* backRng = BlockIterRng(0, IterRngKmin());
    Real givenBackWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenBackWallBoundaryVars, backRng,
                          Vertex_NoslipEQN);

    int* frontRng = BlockIterRng(0, IterRngKmax());
    Real givenFrontWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
    TreatBlockBoundary3D(givenFrontWallBoundaryVars, frontRng,
                         Vertex_NoslipEQN);
}
*/

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

// void ForwardEuler() {
//    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
//        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
//        ops_par_loop(KerCutCellCVTUpwind, "KerCutCellCVTUpwind2nd",
//                     g_Block[blockIndex], SPACEDIM, iterRng,
//                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
//                                 ONEPTREGULARSTENCIL, "double", OPS_READ),
//                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL,
//                     "int",
//                                 OPS_READ),
//                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
//                                 LOCALSTENCIL, "int", OPS_READ),
//                     ops_arg_dat(g_f[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
//                                 "double", OPS_READ),
//                     ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL,
//                                 "double", OPS_RW));
//        Real schemeCoeff{1};
//        ops_par_loop(KerCutCellExplicitTimeMach, "KerCutCellExplicitTimeMach",
//                     g_Block[blockIndex], SPACEDIM, iterRng,
//                     ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
//                     ops_arg_gbl(&schemeCoeff, 1, "double", OPS_READ),
//                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL,
//                     "int",
//                                 OPS_READ),
//                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
//                                 LOCALSTENCIL, "int", OPS_READ),
//                     ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL,
//                                 "double", OPS_READ),
//                     ops_arg_dat(g_feq[blockIndex], NUMXI, LOCALSTENCIL,
//                                 "double", OPS_READ),
//                     ops_arg_dat(g_Tau[blockIndex], NUMCOMPONENTS,
//                     LOCALSTENCIL,
//                                 "double", OPS_READ),
//                     ops_arg_dat(g_Bodyforce[blockIndex], NUMXI, LOCALSTENCIL,
//                                 "double", OPS_READ),
//                     ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL,
//                     "double",
//                                 OPS_RW));
//    }
//}

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
    UpdateMacroVars3D();
    // Real TotalMass{0};
    // CalcTotalMass(&TotalMass);
    // Real Ratio{TotalMass/TotalMeshSize()};
    // NormaliseF(&Ratio);
    // UpdateMacroVars();
    CopyDistribution3D(g_f, g_fStage);
    UpdateFeqandBodyforce3D();
    UpdateTau3D();
    Collision3D();
    Stream3D();
    ops_halo_transfer(HaloGroup());
    ImplementBoundaryConditions();
}

// void TimeMarching() {
//    UpdateMacroVars();
//    //CopyDistribution(g_f, g_fStage);
//    UpdateSWEFeqandBodyforce();
//    UpdateSWETau();
//    ForwardEuler();
//    //ops_halo_transfer(HaloGroups);
//    ImplementBoundary();
// }

// void TimeMarching() {
//     UpdateMacroVars();
//     //CopyDistribution(g_f, g_fStage);
//     UpdateFeqandBodyforce();
//     UpdateTau();
//     ForwardEuler();
//     //ops_halo_transfer(HaloGroups);
//     ImplementBoundary();
// }

#endif /* OPS_3D */
