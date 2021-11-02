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
#include "Poiseuille_kernel.inc"
//Provide macroscopic initial conditions
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

void UpdateMacroscopicBodyForce(const Real time) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};

            ops_par_loop(KerUpdateMacroBodyForce, "KerUpdateMacroBodyForce",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx), 3,
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());
        }
    }
}

void simulate() {

    std::string caseName{"Poiseuille"};
    SizeType spaceDim{3};
    DefineCase(caseName, spaceDim);
    std::vector<int> blockIds{0};
    std::vector<std::string> blockNames{"Cavity"};
    std::vector<int> blockSize{33, 33, 33};
    Real meshSize{1. / 32};
    std::map<int, std::vector<Real>> startPos{{0, {0.0, 0.0, 0.0}}};
    DefineBlocks(blockIds, blockNames, blockSize, meshSize, startPos);


    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d3q19"};
    std::vector<Real> tauRef{0.01};
    DefineComponents(compoNames, compoid, lattNames, tauRef);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                                 Variable_V, Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<CollisionType> collisionTypes{Collision_BGKIsothermal2nd};
    std::vector<int> collisionCompoId{0};
    DefineCollision(collisionTypes, collisionCompoId);

    std::vector<BodyForceType> bodyForceTypes{BodyForce_1st};
    std::vector<SizeType> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SchemeType scheme{Scheme_StreamCollision};
    DefineScheme(scheme);

    // Setting boundary conditions
    SizeType componentId{0};
    std::vector<VariableTypes> macroVarTypesatBoundary{Variable_U, Variable_V,
                                                       Variable_W};
    std::vector<Real> noSlipStationaryWall{0, 0, 0};

    //Top and Bottom No Slip Walls
    DefineBlockBoundary(0, componentId, BoundarySurface::Top,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Bottom,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    //Periodic Boundary Conditions
    DefineBlockBoundary(0, componentId, BoundarySurface::Front,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Back,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall);


    std::vector<InitialType> initType{Initial_BGKFeq2nd};
    std::vector<SizeType> initalCompoId{0};
    DefineInitialCondition(initType,initalCompoId);
    Partition();
    SetInitialMacrosVars();
    PreDefinedInitialCondition3D();
    SetTimeStep(meshSize/ SoundSpeed());

    const Real convergenceCriteria{1E-5};
    const SizeType checkPeriod{200};
    Iterate(StreamCollision,convergenceCriteria, checkPeriod);
}

void simulate(const Configuration & config, const SizeType timeStep=0) {
    // DefineCase(config.caseName, config.spaceDim);
    // DefineBlocks(config.blockNum, config.blockSize, config.meshSize,
    //                     config.startPos);
    // if (timeStep == 0) {
    //     DefineComponents(config.compoNames, config.compoIds, config.lattNames);
    //     DefineMacroVars(config.macroVarTypes, config.macroVarNames,
    //                     config.macroVarIds, config.macroCompoIds);
    // } else {
    //     // restart from a time step
    //     DefineComponents(config.compoNames, config.compoIds, config.lattNames,
    //                      timeStep);
    //     DefineMacroVars(config.macroVarTypes, config.macroVarNames,
    //                     config.macroVarIds, config.macroCompoIds,timeStep);
    // }

    // DefineCollision(config.CollisionTypes, config.CollisionCompoIds);
    // DefineBodyForce(config.bodyForceTypes, config.bodyForceCompoIds);
    // DefineScheme(config.schemeType);
    // DefineInitialCondition(config.initialTypes,config.initialConditionCompoId);
    // for (auto& bcConfig : config.blockBoundaryConfig) {
    //     DefineBlockBoundary(bcConfig.blockIndex, bcConfig.componentID,
    //                         bcConfig.boundarySurface, bcConfig.boundaryScheme,
    //                         bcConfig.macroVarTypesatBoundary,
    //                         bcConfig.givenVars, bcConfig.boundaryType);
    // }
    // Partition();
    // if (timeStep == 0) {
    //     SetInitialMacrosVars();
    //     PreDefinedInitialCondition3D();
    // } else{
    //     //Help function for restart a steady simulation
    //     //Mainly make the residual calculation correct at first iteration.
    //     RestartMacroVars4SteadySim();
    // }

    // SetTauRef(config.tauRef);
    // SetTimeStep(config.meshSize / SoundSpeed());
    // if (config.transient){
    //     Iterate(config.timeSteps, config.checkPeriod,timeStep);
    // } else{
    //     Iterate(config.convergenceCriteria, config.checkPeriod,timeStep);
    // }
}

int main(int argc, const char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 4);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    // start a simulation by hard-coding
    if (argc <= 1) {
        simulate();
    }
    // start a new simulaton from a configuration file
    if (argc>1 && argc <=2){
        std::string configFileName(argv[1]);
        ReadConfiguration(configFileName);
        simulate(Config());
    }
    // restart from the time step specified by argv[2]
    if (argc>2 && argc <=3){
        std::string configFileName(argv[1]);
        ReadConfiguration(configFileName);
        const SizeType timeStep{static_cast<SizeType>(std::stoi(argv[2]))};
        simulate(Config(),timeStep);
    }

    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    //Print OPS performance details to output stream
    ops_timing_output(std::cout);
    ops_exit();
}