# Usage of Python Script for distribution function

The script name is User_defined_function.py. It needs input from two files namely:-

1. C++ code file written by user which will be used to extract some information. This is the code which the user actually wants to run and contains the main() function.  

   In my case, I have named it to be `lbm3d_hilemms.cpp`. 

   

2. A file where user will write the equation which has to be translated, generate the code and then inserted the code at appropriate place in the MPLB code. This file is named as `Dist_fun_eqn.txt`.

Consider the following example.

```
#Please enter the CPP file with main() function here.
    CppFileName = lbm3d_hilemms.cpp;

Dist_f[CompoVeloSpaIdx(0, 18, 0,0,0)] =
    Weights[CompoVeloIdx(1, 10)] * Macro_Vars[CompoMacroSpaIdx(0, rho, 2,0,0)]*
    (
        Micro_Vel_Cx[CompoVeloIdx(0, 5)] * Macro_Vars[CompoMacroSpaIdx(0, u, 3,0,0)] +

        Micro_Vel_Cy[CompoVeloIdx(0, 5)] * Macro_Vars[CompoMacroSpaIdx(0, v, 4,0,0)] +

        Micro_Vel_Cz[CompoVeloIdx(1, 3)] * Macro_Vars[CompoMacroSpaIdx(0, w, 5,0,0)] +

        Macro_Vars[CompoMacroSpaIdx(1,T,0,0,0)] +

        { Macro_Vars[CompoMacroSpaIdx(1,T,0,0,0)] } ^ {3} +

        {Micro_Vel_Cy[CompoVeloIdx(0, 5)]}^{4} +

        ( Macro_Vars[CompoMacroSpaIdx(0, u, 0,0,0)] - Macro_Vars[CompoMacroSpaIdx(0, u, -1,0,0)] ) /
        ( Coord_X[SpaIndex(0,0,0)] - Coord_X[SpaIndex(-1,0,0)])

    );
```



This is the equation which the use will write. Some conventions used are:-

1. For Distribution function, a variable has to be preceded by Dist_ (E.g Dist_f)
2. For microscopic velocity, a variable has to be preceded by Micro_Vel_. Then a user can choose any of the three velocity magnitudes in x, y and z directions (E.g. Micro_Vel_Cx).
3. For macroscopic variables, the variable has to be preceded by Macro_Vars (e.g. `Macro_Vars[CompoMacroSpaIdx(0, u, 3,0,0)`. 
4. For referring to coordinates, a variable `X, Y or Z` has to be preceded by Coord_ (e.g. `Coord_X`).
5. For referring to Weights, use the string literal `Weights` (e.g. `Weights[CompoVeloIdx(1, 10)]`) .
6. For raising something to a power, enclose the expression in {} followed by ^{m} where m is power to which expression is to be raised (E.g. `{ Macro_Vars[CompoMacroSpaIdx(1,T,0,0,0)] } ^ {3}`).

For this example, the file `lbm3d_hilemms.cpp` was defined as follows:-

```c++
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

// Code_modification needed
// Currently defining OPS 3d here. We need some mechanism to generate this
// automatically.

extern int HALODEPTH;

void simulate() {

    std::string caseName{"3D_lid_Driven_cavity"};
    int spaceDim{3};
    DefineCase(caseName, spaceDim);

    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d3q19","d3q27"};
    DefineComponents(compoNames, compoid, lattNames);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w", "T"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0, 1};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<EquilibriumType> equTypes{Equilibrium_BGKIsothermal2nd};
    std::vector<int> equCompoId{0};
    DefineEquilibrium(equTypes, equCompoId);

    std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
    std::vector<int> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SetupScheme();
    SetupBoundary();

    //Setting boundary conditions
    int blockIndex{0};
    int componentId{0};
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};

    BoundaryType boundType[6] = {
        BoundaryType_Periodic,       BoundaryType_Periodic,
        BoundaryType_Periodic,       BoundaryType_Periodic,    
        BoundaryType_Periodic,       BoundaryType_Periodic};

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

    ops_printf("Block boundary defined!\n");
    int blockNum{1};
    std::vector<int> blockSize{64, 64, 64};
    Real meshSize{1. / 63};
    std::vector<Real> startPos{0.0, 0.0, 0.0};
    DefineProblemDomain(blockNum, blockSize, meshSize, startPos);

    std::vector<Real> initialMacroValues{1, 0, 0, 0};
    DefineInitialCondition(blockIndex, componentId, initialMacroValues);

    std::vector<Real> tauRef{0.01};
    SetTauRef(tauRef);

    SetTimeStep(meshSize / SoundSpeed());



    // currently this information is not playing major role in this
    // implementation.
    SchemeType scheme{stStreamCollision};
    const int steps{10000};
    const int checkPeriod{500};
    Iterate(scheme, steps, checkPeriod);

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
```



In the terminal, then type `python User_defined_function.py`, we will get the following output.

```c++
f[OPS_ACC_MD2(18,0,0,0)] =
    WEIGHTS[29] * macroVars[OPS_ACC_MD1(0,2,0,0)]*
    (
        XI[15] * macroVars[OPS_ACC_MD1(1,3,0,0)] +

        XI[16] * macroVars[OPS_ACC_MD1(2,4,0,0)] +

        XI[68] * macroVars[OPS_ACC_MD1(3,5,0,0)] +

        macroVars[OPS_ACC_MD1(4,0,0,0)] +

        pow( macroVars[OPS_ACC_MD1(4,0,0,0)] ,3)+

        pow(XI[16],4)+

        ( macroVars[OPS_ACC_MD1(1,0,0,0)] - macroVars[OPS_ACC_MD1(1,-1,0,0)] ) /
        ( coordinates[OPS_ACC_MD0(0,0,0,0)] - coordinates[OPS_ACC_MD0(0,-1,0,0)])

    );
```



Note:- The example provided here in this manual is correct and will work fine if you use the python script. If you decide to clone the repository from gitlab, then the file `lbm3d_hilemms.cpp` or the equation in file `Dist_fun_eqn.txt` has to be changed accordingly.