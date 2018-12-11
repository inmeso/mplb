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
#include "evolution.h"
#include "flowfield.h"
#include "model.h"
#include "ops_seq.h"
#include "scheme.h"
#include "type.h"
void simulate() {
    SetupModel();
    SetupScheme();
    SetupBoundary();
    SetupFlowfieldfromHdf5();
    InitialiseSolution();
    const int maxIter =2500;
    const int checkPeriod = 500;
    for (int iter = 1; iter <= maxIter; iter++) {
        //StreamCollision();//Stream-Collision scheme
        TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % checkPeriod) == 0) {
            //#ifdef debug
            UpdateMacroVars();
            CalcResidualError();
            DispResidualError(iter,checkPeriod*TimeStep());
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
    ops_init(argc, argv, 5);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    simulate();
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n",et1-et0);
    //Print OPS performance details to output stream
    ops_timing_output(stdout);
    ops_exit();
}

