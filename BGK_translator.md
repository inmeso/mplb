# Python Script For Translating and Inserting BGK Equation

The Script file is named as `Single_eqn_multi_xi_dist_fun.py`. For an introductory example, let us assume that we have only one single component. Let us also suppose that the user is using `DnQm` lattice model. In this case, the equation for distribution function can be specified in the following 3 ways:-

1.  A user might have to write `m` equations for the distribution function (This is the worst possible scenario where the value of distribution function i.e. `f` has different equation for each velocity index). 

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

   