# 1. Usage of Python Script for distribution function (Version 1.0)

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
#include "evolution.h"
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
    DefineBlocks(blockNum, blockSize, meshSize, startPos);

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





# 2. Python Script For Translating and Inserting BGK Equation (Version 2.0)

The Script file is named as `Single_eqn_multi_xi_dist_fun.py`. For an introductory example, let us assume that we have only one single component. Let us also suppose that the user is using `DnQm` lattice model. In this case, the equation for distribution function can be specified in the following 3 ways:-

1. A user might have to write `m` equations for the distribution function (This is the worst possible scenario where the value of distribution function i.e. `f` has different equation for each velocity index). 

2. There is one single equation for all `m` lattice nodes. In this case, a user will have to insert the loop of the type:-

   ```C++
   for (XiIndex = 0; XiIndex<= MaxXiIndex; XiIndex++)
   {
   	User_Defined_Eqn;
   }
   ```

3. Another way, which is to allow user to specify the range in the equation itself (See details in the box below). In this case, we are preventing user from inserting the loop and conditional statements in the code by his own. Instead a Python code will generate the appropriate code for the same.

   ```C++
   // The format to specify the didtribution function in general (without range) is:-
   Dist_f[CompoVeloSpaIdx(CompoID, XiIdx, RelPos_x, RelPos_y, RelPos_z)]
   
   //If we have to specify a range of components, we use following syntax:-
   Dist_f[CompoVeloSpaIdx({a:b}~{c|d|e},{e:f}~{g|h},RelPos_x, RelPos_y, RelPos_z)]
   
   // where
   a = Min Component ID to start from
   
   b = Max Component ID 
   
   {a:b} --> means the equation has to be used for all componenst with ID from a to b.
   
   ~ = Is used to specify a list of ID for which the equation is not valid.
   
   {c|d|e} --> Means that the equation is not valid for component ID c, d and e.
   
   // Similar logic is used for velocity ID.
   ```

   For the case of BGK equation with just one component, we use the following.

   ```C++
   // Assuming D3Q19 lattice
   Dist_f[CompoVeloSpaIdx({0:0}~{},{0:18}~{}, 0, 0, 0)]
   
   /************************************************************************************
   1. Here a=b means that there is just one single component ID.
   2. Since there is one single component, there cannot be anyhting in the list of ID's  to be exluded.
   3. We are using D3Q19 lattice, therefore the velocity range is specified as {0:18}.
   4. Since no velocity ID has a special/different equation, the list of velocity ID's to be excluded is kept empty.
   *************************************************************************************/
   ```

   

The translator Python script needs information from two files:-

1. C++ code file written by user which will be used to extract some information. This is the code which the user actually wants to run and contains the main() function.  

   In my case, I have named it to be `lbm3d_hilemms.cpp`. 

   

2. A file where user will write the equation which has to be translated, generate the code and then inserted the code at appropriate place in the MPLB code. This file is named as `Dist_fun_eqn.txt`.

Consider the following example for BGK equation.

```
#Please enter the CPP file with main() function here.
    CppFileName = lbm3d_hilemms.cpp;

Dist_f[CompoVeloSpaIdx({0:0}~{}, {0:18}~{}, 0,0,0)] =
    Weights[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(,rho, 0,0,0)]*
    
    (1.0 + {3}^{0.5} * 
    
    	   ( Micro_Vel_Cx[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, u, 0,0,0)] +
             Micro_Vel_Cy[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, v, 0,0,0)] +
             Micro_Vel_Cz[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, w, 0,0,0)]
		   ) +

      1.5 *
            ( Micro_Vel_Cx[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, u, 0,0,0)] +
              Micro_Vel_Cy[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, v, 0,0,0)] +
              Micro_Vel_Cz[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, w, 0,0,0)]
            ) *

            ( Micro_Vel_Cx[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, u, 0,0,0)] +
              Micro_Vel_Cy[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, v, 0,0,0)] +
              Micro_Vel_Cz[CompoVeloIdx(,)] * Macro_Vars[CompoMacroSpaIdx(, w, 0,0,0)]
            ) -

   0.5 *
   ( Macro_Vars[CompoMacroSpaIdx(, u, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(, u, 0,0,0)] +
     Macro_Vars[CompoMacroSpaIdx(, v, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(, v, 0,0,0)] +
     Macro_Vars[CompoMacroSpaIdx(, w, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(, w, 0,0,0)]
   ) 
   );
```

Here we have assumed that the current `Xi Index` and `Component ID` is same for all variables occurring on Right hand side of the equation as compared to the left hand side of the equation. Therefore, many of the values are just blank spaces.   



The Python Script has 4 important roles to fulfil:-

1. Code  Translation
2. User Defined Function (UDF) Definition
3. UDF Declaration
4. UDF Call insertion into code.

#### 1. Code Translation

On running the translator script by typing `Pyhton Single_eqn_multi_xi_dist_fun.py` into the terminal we get:-

```C++
Real result;
std::vector<int> CompIdsExcluded;
std::vector<int> VelIdsExcluded;
if(XiIdx >= 0 && XiIdx <= 18)
{
 result  =

 WEIGHTS[XiIdx] * macroVars[OPS_ACC_MD1(0,0,0,0)]*
 (
   1.0 +
   pow(3,0.5)* 
   ( 
     XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +
     XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +
     XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]
   ) +

   1.5 *
   (
     XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +
     XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +
     XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]                
   ) *

   (
     XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +
     XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +
     XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]
    ) -

    0.5 *
    (
      macroVars[OPS_ACC_MD1(1,0,0,0)] * macroVars[OPS_ACC_MD1(1,0,0,0)] +
      macroVars[OPS_ACC_MD1(2,0,0,0)] * macroVars[OPS_ACC_MD1(2,0,0,0)] +
      macroVars[OPS_ACC_MD1(3,0,0,0)] * macroVars[OPS_ACC_MD1(3,0,0,0)]
    ) 
    );
}
```

The vectors `CompIdsExcluded` and `VelIdsExcluded` are empty because we did not specify any ID which has to excluded. Otherwise these vectors would have been later used to check whether the ID is in the Exclusion list or not. 



#### 2. UDF Function Definition

The Pyhton code generates the following UDF definition after it finishes the translation successfully.

```C++
Real CalcUDFFeqNew (const int XiIdx, const Real* macroVars, const int polyOrder)
{
	// Above Tranlated Code.
	return result;
}
```

Result is the variable that stores the value of the distribution function and is the value returned from the function. This function definition is then inserted into `model.cpp` just before the script finds the `#include "model_kernel.h"` into the file.



#### 3. UDF Function Declaration

The python code generates the following for the UDF declaration.

```C++
Real CalcUDFFeqNew (const int XiIdx, const Real* macroVars, const int polyOrder)
```

This is inserted into the the file `model.h` just before the text `#endif`.



#### 4. UDF call insertion into HiLeMMS code

The python code looks for all the `kernel.h` files in the current folder to find the existing call `CalcBGKFeq(Original Arguements)` and replaces it with `CalcUDFFeqNew(New Arguements)`. For details, about the arguments, see the box below.

```C++
// Original Call Taken from model.h
CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)

// After replcement, it becomes
CalcUDFFeqNew(xiIndex, macroVars, polyOrder)

// In each case, we have kept the first and last arguement intact so as to keep the original logic of MPLB code intact.
```

   

# 3. BGK Translator (Version 2.1) 

A. As compared to Version 2.0, it was found that the mechanism used to allow user to specify a range using colons `:` and exclusion list using `~` symbol before `{}` was not good as it might cause some confusion to the user. Moreover, it was felt that such mechanism will complicate the task of Pyhton translator.

B. It was decided to allow a user to write his own loops  to specify a range. See the input and the output for the new translator for the BGK Equation.

####  INPUT

```C++
//Please enter the CPP file with main() function here.
//    CppFileName = lbm3d_hilemms.cpp;

for(int XiIndex =0; XiIndex<=18; XiIndex++ )
{

Dist_f[CompoVeloSpaIdx(0,XiIndex,0,0,0)] =
    
    Weights[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0,0, 0,0,0)]*
    (
    1.0 +
    
    {3}^{0.5} * 
       ( 
       Micro_Vel_Cx[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 1, 0,0,0)] +
       Micro_Vel_Cy[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 2, 0,0,0)] +
       Micro_Vel_Cz[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 3, 0,0,0)]
       ) +
    
       1.5 *
       (
       Micro_Vel_Cx[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 1, 0,0,0)] +
       Micro_Vel_Cy[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 2, 0,0,0)] +
       Micro_Vel_Cz[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 3, 0,0,0)]
       ) *
                    
       (
       Micro_Vel_Cx[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 1, 0,0,0)] +
       Micro_Vel_Cy[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 2, 0,0,0)] +
       Micro_Vel_Cz[CompoVeloIdx(0,XiIndex)] * Macro_Vars[CompoMacroSpaIdx(0, 3, 0,0,0)]
       ) -

   0.5 *
   (
   Macro_Vars[CompoMacroSpaIdx(0, 1, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(0, 1, 0,0,0)] +
   Macro_Vars[CompoMacroSpaIdx(0, 2, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(0, 2, 0,0,0)] +
   Macro_Vars[CompoMacroSpaIdx(0, 3, 0,0,0)] * Macro_Vars[CompoMacroSpaIdx(0, 3, 0,0,0)]
   ) 

   );

} 
//End of for loop. 
```



On running the following command in the terminal, we will get the following output.

#### OUTPUT

```C++
for(int XiIndex =0; XiIndex<=18; XiIndex++ )
{

f[OPS_ACC_MD2(COMPOINDEX[2 *0] + XiIndex,0,0,0)] =
        
WEIGHTS[COMPOINDEX[2 *0] + XiIndex] * 
    macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 0,0,0,0)]*
    
    (
     1.0 +
     pow(3,0.5)* 
      ( 
       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * 		                                       macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +
          
       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] *     		     						macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +
          
       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * 					        			macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]

       ) +

      1.5 *
      (
       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * 										 	macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * 				 				 		 macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * 			 						 	macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
      ) *

      (
       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * 		 					 				macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * 										macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

       XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * 										macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
       ) -

       0.5 *
       (
        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] * 								macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] * 								macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)] * 								macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
       ) 

    );

} 
//End of for loop. 
```

To Do: Task to insert the code at appropriate place by python code has to be done. Although the generated code was manually inserted to check the correctness and was found fine. 



# 4. Body Force 1st Order (Code Generation and Insertion)

Consider an example where a user gives an equation to calculate the first order body force term. Our aim is to translate and insert the equation into the MPLB code. 

As an example, consider that the user gives the following equation.

####  INPUT

```C++
Real acceleration[]{0.0001, 0, 0};

Real g_dot_c{0.0};

for(int XiIndex =0; XiIndex<=18; XiIndex++ )
{
    g_dot_c = acceleration[0] * Micro_Vel_Cx[CompoVeloIdx(0,XiIndex)] +

              acceleration[1] * Micro_Vel_Cy[CompoVeloIdx(0,XiIndex)] +

              acceleration[2] * Micro_Vel_Cz[CompoVeloIdx(0,XiIndex)];


Force_BodyForce[CompoVeloSpaIdx(0,XiIndex,0,0,0)] =  
			                Weights[CompoVeloIdx(0,XiIndex)] * 
                            Macro_Vars[CompoMacroSpaIdx(0,0, 0,0,0)] *
                            g_dot_c; 

            
}
```



After running the python script, we get the following output.

#### OUTPUT

```C++
{
 Real acceleration[]{0.0001, 0, 0};

 Real g_dot_c{0.0};

 for(int XiIndex =0; XiIndex<=18; XiIndex++ )
 {
    g_dot_c = acceleration[0] * XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * CS +

              acceleration[1] * XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * CS +

              acceleration[2] * XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * CS;


	bodyForce[OPS_ACC_MD4(COMPOINDEX[2 *0] + XiIndex,0,0,0)] =  
    						WEIGHTS[COMPOINDEX[2 *0] + XiIndex] * 
                            macroVars[OPS_ACC_MD3(VARIABLECOMPPOS[2 * 0] + 0,0,0,0)] *
                            g_dot_c; 
          
 }
}
```



