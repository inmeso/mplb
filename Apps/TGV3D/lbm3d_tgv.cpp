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
#include <fstream>
#include <string>
#include "mplb.h"
#include "ops_seq_v2.h"
#include "tgv3d_kernel.inc"
long long totalMeshSize;
ops_reduction KeHandle;
ops_reduction KedHandle;
ops_reduction EnHandle;
std::ofstream datCentral2nd, datCentral4th, datCentral6th, datKed;

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

void CalcTurburlentQuantities4th() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(
                KerCalcTurburlentQuantities4th,
                "KerCalcTurburlentQuantities4th", block.Get(), SpaceDim(),
                iterRng.data(), ops_arg_reduce(KeHandle, 1, "double", OPS_INC),
                ops_arg_reduce(KedHandle, 1, "double", OPS_INC),
                ops_arg_reduce(EnHandle, 1, "double", OPS_INC),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                            TWOPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                            TWOPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx), 1,
                            TWOPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(&Config().meshSize, 1, "double", OPS_READ));
        }
    }
}

void CalcTurburlentQuantities6th() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(
                KerCalcTurburlentQuantities6th,
                "KerCalcTurburlentQuantities6th", block.Get(), SpaceDim(),
                iterRng.data(), ops_arg_reduce(KeHandle, 1, "double", OPS_INC),
                ops_arg_reduce(KedHandle, 1, "double", OPS_INC),
                ops_arg_reduce(EnHandle, 1, "double", OPS_INC),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                            THREEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                            THREEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx), 1,
                            THREEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(&Config().meshSize, 1, "double", OPS_READ));
        }
    }
}

void CalcTurburlentQuantities2nd() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(
                KerCalcTurburlentQuantities2nd,
                "KerCalcTurburlentQuantities2nd", block.Get(), SpaceDim(),
                iterRng.data(), ops_arg_reduce(KeHandle, 1, "double", OPS_INC),
                ops_arg_reduce(KedHandle, 1, "double", OPS_INC),
                ops_arg_reduce(EnHandle, 1, "double", OPS_INC),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                            ONEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                            ONEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx), 1,
                            ONEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(&Config().meshSize, 1, "double", OPS_READ));
        }
    }
}


// Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time) {}

void TGVIterate(const SizeType steps, const SizeType checkPointPeriod,
                const SizeType start) {
    const SchemeType scheme = Scheme();
    switch (scheme) {
        case Scheme_StreamCollision: {
            for (SizeType iter = start; iter < start + steps; iter++) {
                const Real time{iter * TimeStep()};
                StreamCollision(time);
                if (((iter + 1) % checkPointPeriod) == 0) {
                    Real ked2nd, ked4th, ked6th;
#ifdef OPS_3D
                    UpdateMacroVars3D();
                    if ((90-TimeStep())<=(iter+1) *TimeStep()) && (iter+1) *TimeStep())<=(90+TimeStep())){
                            WriteFieldsToHdf5(iter + 1);
                        }
                    g_MacroVars().at(g_Components().at(0).uId).TransferHalos();
                    g_MacroVars().at(g_Components().at(0).vId).TransferHalos();
                    g_MacroVars().at(g_Components().at(0).wId).TransferHalos();
                    CalcTurburlentQuantities4th();
                    Real enRes, keRes, kedRes;
                    ops_reduction_result(KeHandle, &keRes);
                    ops_reduction_result(KedHandle, &kedRes);
                    ops_reduction_result(EnHandle, &enRes);
                    ked4th = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
                    if (ops_my_global_rank == MPI_ROOT) {
#endif
                        datCentral4th
                            << TimeStep() * (iter + 1) << " "
                            << keRes / totalMeshSize / 0.01 << " " << ked4th
                            << " " << enRes / totalMeshSize / 0.01 << std::endl;
#ifdef OPS_MPI
                    }
#endif
                    CalcTurburlentQuantities2nd();
                    ops_reduction_result(KeHandle, &keRes);
                    ops_reduction_result(KedHandle, &kedRes);
                    ops_reduction_result(EnHandle, &enRes);
                    ked2nd = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
                    if (ops_my_global_rank == MPI_ROOT) {
#endif
                        datCentral2nd
                            << TimeStep() * (iter + 1) << " "
                            << keRes / totalMeshSize / 0.01 << " " << ked2nd
                            << " " << enRes / totalMeshSize / 0.01 << std::endl;
#ifdef OPS_MPI
                    }
#endif

                    CalcTurburlentQuantities6th();
                    ops_reduction_result(KeHandle, &keRes);
                    ops_reduction_result(KedHandle, &kedRes);
                    ops_reduction_result(EnHandle, &enRes);
                    ked6th = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
                    if (ops_my_global_rank == MPI_ROOT) {
#endif
                        datCentral6th
                            << TimeStep() * (iter + 1) << " "
                            << keRes / totalMeshSize / 0.01 << " " << ked6th
                            << " " << enRes / totalMeshSize / 0.01 << std::endl;
#ifdef OPS_MPI
                    }
#endif

#endif
#ifdef OPS_2D
                    UpdateMacroVars();
#endif
                }
            }
        } break;
        default:
            break;
    }
    DestroyModel();
}

void simulate(const Configuration& config) {
    DefineCase(config.caseName, config.spaceDim, config.transient);
    DefineBlocks(config.blockIds, config.blockNames, config.blockSize,
                 config.meshSize, config.startPos);
#ifdef OPS_MPI
    if (ops_my_global_rank == MPI_ROOT) {
#endif
        datCentral2nd.open(CaseName() +
                           std::to_string(g_Block().at(0).Size().at(0)) +
                           "2nd.dat");
        datCentral4th.open(CaseName() +
                           std::to_string(g_Block().at(0).Size().at(0)) +
                           "4th.dat");
        datCentral6th.open(CaseName() +
                           std::to_string(g_Block().at(0).Size().at(0)) +
                           "6th.dat");
        datCentral2nd.precision(16);
        datCentral4th.precision(16);
        datCentral6th.precision(16);
        datCentral4th
            << "#Time kineticEnergy KineticEnergyDissipation Enstrophy "
            << std::endl;
        datCentral2nd
            << "#Time kineticEnergy KineticEnergyDissipation Enstrophy "
            << std::endl;
        datCentral6th
            << "#Time kineticEnergy KineticEnergyDissipation Enstrophy "
            << std::endl;
#ifdef OPS_MPI
    }
#endif

    KeHandle = ops_decl_reduction_handle(sizeof(Real), "double", "KeHandle");
    KedHandle = ops_decl_reduction_handle(sizeof(Real), "double", "KedHandle");
    EnHandle = ops_decl_reduction_handle(sizeof(Real), "double", "EnHandle");

    totalMeshSize =
        config.blockSize[0] * config.blockSize[1] * config.blockSize[2];
    DefineBlockConnection(config.fromBlockIds, config.fromBoundarySurface,
                          config.toBlockIds, config.toBoundarySurface,
                          config.blockConnectionType);
    DefineComponents(config.compoNames, config.compoIds, config.lattNames,
                     config.tauRef, config.currentTimeStep);
    DefineMacroVars(config.macroVarTypes, config.macroVarNames,
                    config.macroVarIds, config.macroCompoIds,
                    config.currentTimeStep, 3);
    g_MacroVars().at(g_Components().at(0).uId).CreateHalos();
    g_MacroVars().at(g_Components().at(0).vId).CreateHalos();
    g_MacroVars().at(g_Components().at(0).wId).CreateHalos();
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
    Partition();
    ops_diagnostic_output();
    if (config.currentTimeStep == 0) {
        SetInitialMacrosVars();
        TGVInitialCondition();
    };
    g_MacroVars().at(g_Components().at(0).uId).TransferHalos();
    g_MacroVars().at(g_Components().at(0).vId).TransferHalos();
    g_MacroVars().at(g_Components().at(0).wId).TransferHalos();
    Real ked2nd, ked4th, ked6th;
    CalcTurburlentQuantities4th();
    Real enRes, keRes, kedRes;
    ops_reduction_result(KeHandle, &keRes);
    ops_reduction_result(KedHandle, &kedRes);
    ops_reduction_result(EnHandle, &enRes);
    ked4th = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
    if (ops_my_global_rank == MPI_ROOT) {
#endif
        datCentral4th << 0.0 << " " << keRes / totalMeshSize / 0.01 << " "
                      << ked4th << " " << enRes / totalMeshSize / 0.01
                      << std::endl;
#ifdef OPS_MPI
    }
#endif
    CalcTurburlentQuantities2nd();
    ops_reduction_result(KeHandle, &keRes);
    ops_reduction_result(KedHandle, &kedRes);
    ops_reduction_result(EnHandle, &enRes);
    ked2nd = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
    if (ops_my_global_rank == MPI_ROOT) {
#endif
        datCentral2nd << 0.0 << " " << keRes / totalMeshSize / 0.01 << " "
                      << ked2nd << " " << enRes / totalMeshSize / 0.01
                      << std::endl;
#ifdef OPS_MPI
    }
#endif

    CalcTurburlentQuantities6th();
    ops_reduction_result(KeHandle, &keRes);
    ops_reduction_result(KedHandle, &kedRes);
    ops_reduction_result(EnHandle, &enRes);
    ked6th = kedRes * Config().tauRef[0] / totalMeshSize / 0.01;
#ifdef OPS_MPI
    if (ops_my_global_rank == MPI_ROOT) {
#endif
        datCentral6th << 0.0 << " " << keRes / totalMeshSize / 0.01 << " "
                      << ked6th << " " << enRes / totalMeshSize / 0.01
                      << std::endl;
#ifdef OPS_MPI
    }
#endif

    //WriteFieldsToHdf5(0);
    SetTimeStep(config.meshSize / SoundSpeed());
    if (config.schemeType == Scheme_StreamCollision_Swap) {
    }

    if (config.schemeType == Scheme_StreamCollision) {
        TGVIterate(config.timeStepsToRun, config.checkPeriod,
                   config.currentTimeStep);
    }
#ifdef OPS_MPI
    if (ops_my_global_rank == MPI_ROOT) {
#endif
        datCentral2nd.close();
        datCentral4th.close();
        datCentral6th.close();
        datKed.close();
#ifdef OPS_MPI
    }
#endif
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