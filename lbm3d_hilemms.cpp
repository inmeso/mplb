// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

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

// Code_modifcication needed
// Currently defining OPS 3d here. We need some mechanism to generate this
// automatically.

using namespace std;

extern int HALODEPTH;

BoundarySurface surface[6] = {BoundSurf_Inlet, BoundSurf_Outlet,
                              BoundSurf_Top,   BoundSurf_Bottom,
                              BoundSurf_Front, BoundSurf_Back};

BoundaryType boundType[6] = {
    BoundType_EQMDiffuseRefl, BoundType_EQMDiffuseRefl,
    BoundType_EQMDiffuseRefl, BoundType_EQMDiffuseRefl,
    BoundType_EQMDiffuseRefl, BoundType_EQMDiffuseRefl};

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

    SetupScheme();
    SetupBoundary();

    int blockNum{1};
    std::vector<int> blockSize{11, 11, 11};
    Real meshSize{0.1};
    std::vector<Real> startPos{0.0, 0.0, 0.0};
    DefineProblemDomain(blockNum, blockSize, meshSize, startPos);

    int blockIndex{0};
    SetupGeomPropAndNodeType(blockIndex, boundType);

    int compoIdInitialCond{0};
    std::vector<Real> initialMacroValues{1, 0, 0, 0};
    DefineIntialCond(blockIndex, compoIdInitialCond, initialMacroValues);
    ops_printf("%s\n", "Flowfield is Initialised now!");

    std::vector<Real> tauRef{0.001};
    SetTauRef(tauRef);

    SetTimeStep(meshSize / SoundSpeed());

    HALODEPTH = HaloPtNum();
    ops_printf("%s\n", "Starting to allocate...");
    DefineHaloTransfer3D();
    // above calls must be before the ops_partition call
    //ops_partition((char*)"LBM");
    ops_printf("%s\n", "Flowfield is setup now!");
    InitialiseSolution3D();

    blockIndex = 0;
    int componentId{0};
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};
    std::vector<Real> inletValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[0], boundType[0],
                        MacroVarsComp, inletValMacroVarsComp);

    std::vector<Real> outletValMacroVarsComp{1, 0, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[1], boundType[1],
                        MacroVarsComp, outletValMacroVarsComp);

    std::vector<Real> topValMacroVarsComp{1, 0.01, 0, 0};
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

#if 0
    // currently this information is not playin major role in this
    // implementation.
    SchemeType scheme{stStreamCollision}; 
    const int steps{10000};
    const int checkPeriod{1000};
    Iterate(scheme, steps, checkPeriod);
#endif

    // if 0
    // currently this information is not playin major role in this
    // implementation.
    SchemeType scheme{stStreamCollision};
    const Real convergenceCriteria{1E-3};
    const int checkPeriod{1000};
    Iterate(scheme, convergenceCriteria, checkPeriod);
    //#endif
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