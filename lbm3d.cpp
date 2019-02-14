// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief Define the main iteration
 *  @author Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "boundary.h"
#include "evolution3d.h"
#include "flowfield.h"
#include "model.h"
#include "ops_seq.h"
#include "scheme.h"
#include "type.h"
void simulate() {
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
    SetupScheme();
    SetupBoundary();
    SetupFlowfieldfromHdf5();
    ops_printf("%s\n", "Flowfield is setupped now!");
    InitialiseSolution3D();
    const int maxIter =100000;
    const int checkPeriod = 1000;
    for (int iter = 0; iter < maxIter; iter++) {
        StreamCollision3D();//Stream-Collision scheme
        //TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % checkPeriod) == 0) {
            //#ifdef debug
            UpdateMacroVars3D();
            CalcResidualError3D();
            DispResidualError3D(iter,checkPeriod*TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
            //if ((densityResidualError + uResidualError+vResidualError) <= 1e-12) break;
            // WriteDistributionsToHdf5(iter);
            // WriteNodePropertyToHdf5(iter);
            //#endif
        }
    }
    DestroyModel();
    DestroyFlowfield();
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


