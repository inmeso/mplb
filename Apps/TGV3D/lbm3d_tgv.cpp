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

/** @brief An example main source code of stimulating 3D lid-driven cavity flow
 *  @author Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "mplb.h"
#include "ops_seq_v2.h"
#include "tgv3d_kernel.inc"
RealField KineticEnergy{"KineticEnergy"};
RealField KineticEnergyDissipation{"KineticEnergyDissipation"};
RealField Enstrophy{"Enstrophy"};
ops_reduction KeHandle{
    ops_decl_reduction_handle(sizeof(Real), "double", "KeHandle")};
ops_reduction KedHandle{
    ops_decl_reduction_handle(sizeof(Real), "double", "KedHandle")};
ops_reduction EnHandle{
    ops_decl_reduction_handle(sizeof(Real), "double", "EnHandle")};
// Provide macroscopic initial conditions
void SetInitialMacrosVars() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars().at(rhoId).at(blockIdx), 1,
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());
        }
    }
}

void TGVInitialCondition() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int compoId{compo.id};
            const Real tau{compo.tauRef};
            ops_par_loop(
                KerTGV3DInit, "KerTGV3DInit", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_WRITE),
                ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex), 1,
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_gbl(&tau, 1, "double", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ));
        }
    }
    // TODO this may be better arranged.
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
#endif  // OPS_3D
}

void CalcTurburlentQuantities() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(
                KerCalcTurburlentQuantities, "KerCalcTurburlentQuantities",
                block.Get(), SpaceDim(), iterRng.data(),
                ops_arg_dat(KineticEnergy.at(blockIdx), 1, TWOPTREGULARSTENCIL,
                            "double", OPS_READ),
                ops_arg_dat(KineticEnergyDissipation.at(blockIdx), 1,
                            TWOPTREGULARSTENCIL, "double", OPS_READ),
                ops_arg_dat(Enstrophy.at(blockIdx), 1, TWOPTREGULARSTENCIL,
                            "double", OPS_READ),
                ops_arg_reduce(EnHandle, 1, "double", OPS_INC),
                ops_arg_reduce(KeHandle, 1, "double", OPS_INC),
                ops_arg_reduce(KedHandle, 1, "double", OPS_INC),
                ops_arg_dat(g_MacroVars().at(rhoId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(&Config().meshSize, 1, "double", OPS_READ));
        }
    }
}

// Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time) {}

void simulate(const Configuration& config) {
    DefineCase(config.caseName, config.spaceDim, config.transient);
    DefineBlocks(config.blockIds, config.blockNames, config.blockSize,
                 config.meshSize, config.startPos);
    DefineBlockConnection(config.fromBlockIds, config.fromBoundarySurface,
                          config.toBlockIds, config.toBoundarySurface,
                          config.blockConnectionType);
    DefineComponents(config.compoNames, config.compoIds, config.lattNames,
                     config.tauRef, config.currentTimeStep);
    DefineMacroVars(config.macroVarTypes, config.macroVarNames,
                    config.macroVarIds, config.macroCompoIds,
                    config.currentTimeStep);
    DefineCollision(config.CollisionTypes, config.CollisionCompoIds);
    DefineBodyForce(config.bodyForceTypes, config.bodyForceCompoIds);
    DefineScheme(config.schemeType);
    DefineInitialCondition(config.initialTypes, config.initialConditionCompoId);
    for (auto& bcConfig : config.blockBoundaryConfig) {
        DefineBlockBoundary(bcConfig.blockIndex, bcConfig.componentID,
                            bcConfig.boundarySurface, bcConfig.boundaryScheme,
                            bcConfig.macroVarTypesatBoundary,
                            bcConfig.givenVars, bcConfig.boundaryType);
    }

    KineticEnergy.SetDataDim(1);
    KineticEnergy.CreateFieldFromScratch(g_Block());
    KineticEnergy.SetDataHalo(2);
    KineticEnergy.CreateHalos();
    KineticEnergyDissipation.SetDataDim(1);
    KineticEnergyDissipation.CreateFieldFromScratch(g_Block());
    KineticEnergyDissipation.SetDataHalo(2);
    KineticEnergyDissipation.CreateHalos();
    Enstrophy.SetDataDim(1);
    Enstrophy.CreateFieldFromScratch(g_Block());
    Enstrophy.SetDataHalo(2);
    Enstrophy.CreateHalos();

    Partition();
    ops_diagnostic_output();
    if (config.currentTimeStep == 0) {
        SetInitialMacrosVars();
        TGVInitialCondition();
    };
    WriteFieldsToHdf5(0);
    SetTimeStep(config.meshSize / SoundSpeed());
    if (config.schemeType == Scheme_StreamCollision_Swap) {
        if (config.transient) {
            Iterate(SwapStreamCollision, config.timeStepsToRun,
                    config.checkPeriod, config.currentTimeStep);
        } else {
            Iterate(SwapStreamCollision, config.convergenceCriteria,
                    config.checkPeriod, config.currentTimeStep);
        }
    }

    if (config.schemeType == Scheme_StreamCollision) {
        if (config.transient) {
            Iterate(StreamCollision, config.timeStepsToRun, config.checkPeriod,
                    config.currentTimeStep);
        } else {
            Iterate(StreamCollision, config.convergenceCriteria,
                    config.checkPeriod, config.currentTimeStep);
        }
    }
}

int main(int argc, const char** argv) {
    // OPS initialisation where a few arguments can be passed to set
    // the simulation
    ops_init(argc, argv, 4);
    bool configFileFound{false};
    std::string configFileName;
    GetConfigFileFromCmd(configFileFound, configFileName, argc, argv);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    // start a new simulaton from a configuration file
    if (configFileFound) {
        ReadConfiguration(configFileName);
        simulate(Config());
    }
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(std::cout);
    ops_exit();
}