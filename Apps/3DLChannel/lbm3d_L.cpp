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
#include "LChannel_kernel.inc"
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
// Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time) {}

void simulate() {
    std::string caseName{"3DLChannel"};
    SizeType spaceDim{3};
    DefineCase(caseName, spaceDim);
    std::vector<int> blockIds{0, 1, 2};
    std::vector<std::string> blockNames{"Top", "Middle", "Right"};
    std::vector<int> blockSize{33, 33, 33, 33, 33, 33, 33, 33, 33};
    Real meshSize{1. / 32};
    std::map<int, std::vector<Real>> startPos{
        {1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}, {0, {0.0, 1.0, 0.0}}};
    DefineBlocks(blockIds, blockNames, blockSize, meshSize, startPos);

    std::vector<int> fromBlockIds{0, 1, 1, 2};
    std::vector<int> toBlockIds{1, 0, 2, 1};
    std::vector<BoundarySurface> fromBoundarySurface{
        BoundarySurface::Bottom, BoundarySurface::Top, BoundarySurface::Right,
        BoundarySurface::Left};
    std::vector<BoundarySurface> toBoundarySurface{
        BoundarySurface::Top, BoundarySurface::Bottom, BoundarySurface::Left,
        BoundarySurface::Right};
    std::vector<VertexType> blockConnectionType{
        VertexType::VirtualBoundary, VertexType::VirtualBoundary,
        VertexType::VirtualBoundary, VertexType::VirtualBoundary};

    DefineBlockConnection(fromBlockIds, fromBoundarySurface, toBlockIds,
                          toBoundarySurface, blockConnectionType);

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

    std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
    std::vector<SizeType> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SchemeType scheme{Scheme_StreamCollision};
    DefineScheme(scheme);

    // Setting boundary conditions
    int componentId{0};
    std::vector<VariableTypes> macroVarTypesatBoundary{Variable_U, Variable_V,
                                                       Variable_W};
    std::vector<Real> noSlipStationaryWall{0, 0, 0};
    // Left noSlipStationaryWall
    DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // Right noSlipStationaryWall
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    std::vector<Real> noSlipMovingWall{0, -0.001, 0};
    DefineBlockBoundary(0, componentId, BoundarySurface::Top,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipMovingWall);
    // bottom noSlipStationaryWall
    DefineBlockBoundary(0, componentId, BoundarySurface::Bottom);
    // front noSlipStationaryWall
    DefineBlockBoundary(0, componentId, BoundarySurface::Front,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // back noSlipStationaryWall
    DefineBlockBoundary(0, componentId, BoundarySurface::Back,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    DefineBlockBoundary(1, componentId, BoundarySurface::Left,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // Right noSlipStationaryWall
    DefineBlockBoundary(1, componentId, BoundarySurface::Right);

    DefineBlockBoundary(1, componentId, BoundarySurface::Top);
    DefineBlockBoundary(1, componentId, BoundarySurface::RightTop,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // bottom noSlipStationaryWall
    DefineBlockBoundary(1, componentId, BoundarySurface::Bottom,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // front noSlipStationaryWall
    DefineBlockBoundary(1, componentId, BoundarySurface::Front,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // back noSlipStationaryWall
    DefineBlockBoundary(1, componentId, BoundarySurface::Back,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    DefineBlockBoundary(2, componentId, BoundarySurface::Left);
    // pressure boundary condition
    std::vector<VariableTypes> pressureType{Variable_Rho};
    std::vector<Real> pressureOutLet{1};
    DefineBlockBoundary(2, componentId, BoundarySurface::Right,
                        BoundaryScheme::ExtrapolPressure1ST, pressureType,
                        pressureOutLet, VertexType::OutLet);

    DefineBlockBoundary(2, componentId, BoundarySurface::Top,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // bottom noSlipStationaryWall
    DefineBlockBoundary(2, componentId, BoundarySurface::Bottom,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // front noSlipStationaryWall
    DefineBlockBoundary(2, componentId, BoundarySurface::Front,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    // back noSlipStationaryWall
    DefineBlockBoundary(2, componentId, BoundarySurface::Back,
                        BoundaryScheme::EQMDiffuseRefl, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    std::vector<InitialType> initType{Initial_BGKFeq2nd};
    std::vector<int> initalCompoId{0};
    DefineInitialCondition(initType, initalCompoId);
    Partition();
    ops_diagnostic_output();
    SetInitialMacrosVars();
    PreDefinedInitialCondition3D();
    SetTimeStep(meshSize / SoundSpeed());

    const Real convergenceCriteria{1E-8};
    const SizeType checkPeriod{1000};
    Iterate(StreamCollision, convergenceCriteria, checkPeriod);
}

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
    Partition();
    ops_diagnostic_output();
    if (config.currentTimeStep == 0) {
        SetInitialMacrosVars();
        PreDefinedInitialCondition3D();
    };
    SetTimeStep(config.meshSize / SoundSpeed());
    if (config.transient) {
        Iterate(config.timeStepsToRun, config.checkPeriod,
                config.currentTimeStep);
    } else {
        Iterate(config.convergenceCriteria, config.checkPeriod,
                config.currentTimeStep);
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
    // start a simulation by hard-coding
    if (!configFileFound) {
        simulate();
    }
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