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
#include "model.h"
#include "ops_seq.h"
#include "scheme.h"
#include "type.h"
#include "hilemms.h"
#include "setup_comput_domain.h"

// Code_modifcication needed
// Currently defining OPS 3d here. We need some mechanism to generate this automatically.
#define OPS_3D

extern int HALODEPTH;

// Face Type
int Bottom_bc_type = 1014;
int Top_bc_type = 1014;
int Left_bc_type = 1014;
int Right_bc_type = 1014;
int Front_bc_type = 1014;
int Back_bc_type = 1014;

int Face_type_bc_array[] = {Left_bc_type,   Right_bc_type,
                            Bottom_bc_type, Top_bc_type,
                            Front_bc_type,  Back_bc_type};

// Edge type
int Left_bottom_bc_type = 1014;
int Left_top_bc_type = 1014;
int Right_bottom_bc_type = 1014;
int Right_top_bc_type = 1014;
int Left_back_bc_type = 1014;
int Left_front_bc_type = 1014;
int Right_back_bc_type = 1014;
int Right_front_bc_type = 1014;
int Bottom_back_bc_type = 1014;
int Bottom_front_bc_type = 1014;
int Top_back_bc_type = 1014;
int Top_front_bc_type = 1014;

int Edge_type_bc_array[] = {Left_bottom_bc_type,  Left_top_bc_type,
                            Right_bottom_bc_type, Right_top_bc_type,
                            Left_back_bc_type,    Left_front_bc_type,
                            Right_back_bc_type,   Right_front_bc_type,
                            Bottom_back_bc_type,  Bottom_front_bc_type,
                            Top_back_bc_type,     Top_front_bc_type};

// Corner type
int Left_bottom_back_bc_type = 1014;
int Left_bottom_front_bc_type = 1014;
int Left_top_back_bc_type = 1014;
int Left_top_front_bc_type = 1014;
int Right_bottom_back_bc_type = 1014;
int Right_bottom_front_bc_type = 1014;
int Right_top_back_bc_type = 1014;
int Right_top_front_bc_type = 1014;

int Corner_type_bc_array[] = {Left_bottom_back_bc_type,  Left_bottom_front_bc_type,
                              Left_top_back_bc_type,     Left_top_front_bc_type,
                              Right_bottom_back_bc_type, Right_bottom_front_bc_type,
                              Right_top_back_bc_type,    Right_top_front_bc_type};

void simulate()
{
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

    /*
    int Num_bound_halo_pts = 1;
    int Halo_Num = 0;
    int Halo_depth = 0; 
    int Scheme_halo_points = 1;
    DefineHaloNumber(Halo_Num, Halo_depth, Scheme_halo_points,Num_bound_halo_pts);
    */

    SetupScheme();
    //cout<<"\n Scheme set up done "<<endl;

    SetupBoundary();
    //cout<<"\n Boundary set up done "<<endl;

    
    int blockNum{1};
    std::vector<int> blockSize{21,21,21};
    Real meshSize{0.05};
    std::vector<Real> startPos{0.0, 0.0, 0.0};
    DefineProblemDomain(blockNum, blockSize, meshSize, startPos);
    
    int Block_index = 0;
    ReadNodeType3D(Block_index, Face_type_bc_array, Edge_type_bc_array, Corner_type_bc_array);

    int blockIndex{0};
    int compoIdInitialCond{0};
    std::vector<Real> initialMacroValues{1,0,0,0};
    DefineIntialCond(blockIndex, compoIdInitialCond, initialMacroValues);
    ops_printf("%s\n", "Flowfield is Initialised now!");

    std::vector<Real> tauRef{0.001};
    SetTauRef(tauRef);

    SetTimeStep(meshSize/SoundSpeed());
    cout<<"DT = "<<TimeStep()<<endl;

    HALODEPTH = HaloPtNum();
    ops_printf("%s\n", "Starting to allocate...");
    DefineHaloTransfer3D();
    // above calls must be before the ops_partition call
    ops_partition((char*)"LBM");
    ops_printf("%s\n", "Flowfield is setup now!");
    InitialiseSolution3D();
    
    SchemeType scheme{stStreamCollision}; //currently this information is not playin major role in this implementation.
    const int steps{10000};
    const int checkPeriod{100};
    Iterate(scheme, steps, checkPeriod);
}


void ImplementBoundary()
{
    //cout<<"\n inside immplement boundary";
    int blockIndex{0};
    int componentId{0};
    BoundarySurface surface{BoundSurf_Inlet};
    BoundaryType type{BoundType_EQMDiffuseRefl};
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,Variable_V, Variable_W};
    std::vector<Real> inletValMacroVarsComp{1,0,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, inletValMacroVarsComp);


    surface = BoundSurf_Outlet;
    std::vector<Real> outletValMacroVarsComp{1,0,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, outletValMacroVarsComp);


    surface = BoundSurf_Top;
    std::vector<Real> topValMacroVarsComp{1,0.01,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, topValMacroVarsComp);


    surface = BoundSurf_Bottom;
    std::vector<Real> bottomValMacroVarsComp{1,0,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, bottomValMacroVarsComp);


    surface = BoundSurf_Front;
    std::vector<Real> frontValMacroVarsComp{1,0,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, frontValMacroVarsComp);


    surface = BoundSurf_Back;
    std::vector<Real> backValMacroVarsComp{1,0,0,0};
    //valuesMacroVarsComp = (Real []) {1,0,0,0};
    DefineBlockBoundary(blockIndex, componentId, surface, type, MacroVarsComp, backValMacroVarsComp);
}

int main(int argc, char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 1);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    simulate();
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n",et1-et0);
    //Print OPS performance details to output stream
    ops_timing_output(stdout);
    ops_exit();
}