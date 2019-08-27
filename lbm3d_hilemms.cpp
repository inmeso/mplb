/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the *following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,    *    this list of conditions and the following disclaimer.
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

/** @brief Define the main iteration
 *  @author Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "boundary.h"
#include "evolution.h"
#include "evolution3d.h"
#include "flowfield.h"
#include "hilemms.h"
#include "model.h"
#include "ops_seq.h"
#include "scheme.h"
#include "setup_comput_domain.h"
#include "type.h"

// Code_modification needed
// Currently defining OPS 3d here. We need some mechanism to generate this
// automatically.

//extern int HALODEPTH;

void simulate() {

    std::string caseName{"3D_lid_Driven_cavity"};
    int spaceDim{3};
    DefineCase(caseName, spaceDim);

    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d3q19"};
    DefineComponents(compoNames, compoid, lattNames);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<EquilibriumType> equTypes{Equilibrium_BGKIsothermal2nd};
    std::vector<int> equCompoId{0};
    DefineEquilibrium(equTypes, equCompoId);

    std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
    std::vector<int> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SchemeType scheme{Scheme_StreamCollision};
    DefineScheme(scheme);
    //Setting boundary conditions
    int blockIndex{0};
    int componentId{0};
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};

    BoundaryType boundType[6] = {
        BoundaryType_EQMDiffuseRefl, BoundaryType_EQMDiffuseRefl,
        BoundaryType_EQMDiffuseRefl, BoundaryType_EQMDiffuseRefl,
        BoundaryType_EQMDiffuseRefl, BoundaryType_EQMDiffuseRefl};

    BoundarySurface surface[6] = {BoundarySurface_Left,  BoundarySurface_Right,
                                  BoundarySurface_Top,   BoundarySurface_Bottom,
                                  BoundarySurface_Front, BoundarySurface_Back};

    std::vector<Real> inletValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[0], boundType[0],
                        MacroVarsComp, inletValMacroVarsComp);

    std::vector<Real> outletValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[1], boundType[1],
                        MacroVarsComp, outletValMacroVarsComp);

    std::vector<Real> topValMacroVarsComp{1, 0, 0, 0.01};
    DefineBlockBoundary(blockIndex, componentId, surface[2], boundType[2],
                        MacroVarsComp, topValMacroVarsComp);

    std::vector<Real> bottomValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[3], boundType[3],
                        MacroVarsComp, bottomValMacroVarsComp);

    std::vector<Real> frontValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[4], boundType[4],
                        MacroVarsComp, frontValMacroVarsComp);

    std::vector<Real> backValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[5], boundType[5],
                        MacroVarsComp, backValMacroVarsComp);
    int blockNum{1};
    std::vector<int> blockSize{33, 33, 33};
    Real meshSize{1. / 32};
    std::vector<Real> startPos{0.0, 0.0, 0.0};
    DefineProblemDomain(blockNum, blockSize, meshSize, startPos);

    //If necessary we can also precisely specify the boundary types
    //for edges and corners
    //must be after DefineProblemDomain

    // BoundaryType leftBottomBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftTopBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightBottomBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightTopBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType bottomBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType bottomFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType topBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType topFrontBcType{BoundaryType_EQMDiffuseRefl};

    // BoundaryType edgeType[12] = {
    //     leftBottomBcType,  leftTopBcType,    rightBottomBcType,
    //     rightTopBcType,    leftBackBcType,   leftFrontBcType,
    //     rightBackBcType,   rightFrontBcType, bottomBackBcType,
    //     bottomFrontBcType, topBackBcType,    topFrontBcType};

    // BoundaryType leftBottomBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftBottomFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftTopBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType leftTopFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightBottomBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightBottomFrontBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightTopBackBcType{BoundaryType_EQMDiffuseRefl};
    // BoundaryType rightTopFrontBcType{BoundaryType_EQMDiffuseRefl};

    // BoundaryType cornerType[8] = {leftBottomBackBcType,  leftBottomFrontBcType,
    //                               leftTopBackBcType,     leftTopFrontBcType,
    //                               rightBottomBackBcType, rightBottomFrontBcType,
    //                               rightTopBackBcType,    rightTopFrontBcType}

    std::vector<Real> initialMacroValues{1, 0, 0, 0};
    DefineInitialCondition(blockIndex, componentId, initialMacroValues);

    // WriteFlowfieldToHdf5(100000);
    // WriteDistributionsToHdf5(100000);
    // WriteNodePropertyToHdf5(100000);

    std::vector<Real> tauRef{0.01};
    SetTauRef(tauRef);

    SetTimeStep(meshSize / SoundSpeed());



    // currently this information is not playing major role in this
    // implementation.

    const int steps{10000};
    const int checkPeriod{1000};
    Iterate(steps, checkPeriod);

    // currently this information is not playing major role in this
    // implementation.

    // SchemeType scheme{stStreamCollision};
    // const Real convergenceCriteria{1E-5};
    // const int checkPeriod{1000};
    // Iterate(convergenceCriteria, checkPeriod);
}

int main(int argc, char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 1);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    simulate();
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(stdout);
    ops_exit();
}