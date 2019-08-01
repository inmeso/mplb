# Manual for the MPLB code

## Introduction

The multi-platform lattice Boltzmann code (MPLB) is a lattice Boltzmann solver written by using the oxford parallel library for structured mesh solvers (OPS).  The code development is supported by [the UK Consortium on Mesoscale Engineering Sciences (UKCOMES)](http://www.ukcomes.org/). By utilizing the functionalities of the OPS library, the code is capable of running on heterogeneous computing platform, and supporting finite-difference  lattice Boltzmann models. We are continually developing new functionalities into the code, including high-order lattice Boltzmann models (finite-difference scheme in the physical space), coupling scheme for solid-fluid two-phase flows etc.  

## Installation

In general, the environment for developing can be setup on any of Windows, Linux or Mac OS system, if we can have a MPI library and the corresponding parallel HDF5 library. For the Windows, we suggest the combination of Windows 10 and its Linux Subsystem, which provides almost the same environment to Linux and Mac OS.

### OPS library

The code relies on the OPS library, which provides the mesh management for parallel computing and the capability of running on heterogeneous computing platform. To facilitate a convenient post-processing, the OPS library also support parallel HDF5 IO capability, which needs the parallel HDF5 library. Other dependencies include OpenMPI/MPICH and one of tools, e.g., OPENCL, CUDA if we would like to use the corresponding hardware. For the detail of the OPS library, we refer to [here](http://www.oerc.ox.ac.uk/projects/ops) where the source code and manual are provided.  

####Installation of HDF5 library

- Windows 10 + WSL 
  We assume an Ubuntu 14.04 system by default, then type
  sudo apt install lbhdf5-openmpi-dev.
  Note: Current no GPU support for WSL
- Mac OS:
  brew install hdf5 --with-mpihttps://www.cvent.com/c/express/c4ff405f-9034-4ab3-a978-c7fe170b9b7f
  *Opensuse:
  sudo zypper in \*\*\*\*
- Configuring the environment
  Using the Mac OS as an example, we need to set up a number of environmental variables. It can be either added into *.bashrc*, or into a script file to be run as "source 'the script file'".
  Example script in *.bashrc*

```bash
#setting the default compiler for openmpi
export OMPI_MPICC=clang
export OMPI_MPICXX=clang++
#setting the default compiler for OPS
export OPS_COMPILER=clang 
#setting the installation direction of MPI, OPS, HDF5, CUDA...
export MPI_INSTALL_PATH=/usr/local
export OPS_INSTALL_PATH=/Users/jpmeng/Documents/work/OPS/ops
export HDF5_INSTALL_PATH=/usr/local
export CUDA_INSTALL_PATH=/Developer/NVIDIA/CUDA-8.0
```

####Compiling the OPS library

```bash
cd $OPS_INSTALL_PATH/c
# compile the sequential version with HDF5 support
make seq
make hdf5_seq
# compile the parallel (mpi) version with HDF5 suuport
make mpi 
make hdf5_mpi
# compile the CUDA version if the tools are installed 
make cuda 
# compile the OPENCL version if the tools are installed
make opencl
# compile the mpi+cuda version
make mpi_cuda
```

### Compiling the MPLB code

The main solver can be compiled in two modes, i.e., the developing mode and the optimized mode. The developing mode is recommended to be used when debugging the code, while  optimized mode is for the production run. However, the developing mode is also fine for production running if we would like to use CPUs only and do not want to be bothered by the issues fhttps://www.cvent.com/c/express/c4ff405f-9034-4ab3-a978-c7fe170b9b7for compiling the optimized version.
####Developing mode

```bash
make lbm2d_dev_seq #(sequential)
make lbm2d_dev_mpi #(parallel)
```



#### Optimized mode

To compile the code in an optimized way, we have to fix two minor issues of the OPS library discussed as follows.

1. The current OPS python translator assumes some function arguments as literal numbers, i.e., not a variable. For instance, the "dim" parameter of the ops\_par\_loop call. To fix this, the python script FixConstantDefinition.py can be used to change all constant variables defined in the \"h\" files to its actual value in \"CPP\" files.
2. The OPS library tends to put source codes for "kernel function" into a ".h" file, which may be confusing to some extent. Therefore, as an intermediate solutiohttps://www.cvent.com/c/express/c4ff405f-9034-4ab3-a978-c7fe170b9b7fn, we put all the kernel codes into a "XXX\_kernel.h" file. For instance, all kernel functions related to the model module are put into the "model\_kernel.h" file. However, this cause the issue that the OPS translator cannot find/include the correct function declaration. Therefore, we provide the python script "FixKernelDeclaration.py".
   For convenience, we provide a bash script to automatize the process. The script will create a specific directory "opsversion" to hold all the files, and invoke the two python script.

#### Pre-processor and Solver

All the information required during the pre-processing stage along with the information required for setting the solver has to be provided in a '.cpp' file. Please see the example given below for details. 

#### Post-processor

There is a simple post-processor written in Python, which can display contour plot and vector plot in both 2D (using matplotlib) and 3D (mayavi) for checking results. The post-processor can also convert the output to the format more friendly to other available software, e.g., plain HDF5 format (readable by Matlab) and  Tecplot HDF5 format. 

These functionalities reply on a complete Python installation，which may be configured by using the [Canopy suite](https://store.enthought.com/downloads/) or the [Anaconda distribution](https://www.anaconda.com/download/). In general either Python 3 or Python 2 will work. 

# 1. 3D lid-driven cavity flow using the stream-collision scheme

For running any simulation using the HiLeMMS interface, a user will have to write his own C++ file. Since the current example is focused on 3D Lid Driven Cavity, let us call the file as `lbm3d_hilemms.cpp`. In this file, the end  user is responsible for defining all the necessary simulation parameters such as number of spatial dimensions, lattice type for the simulation, number of components, number of macroscopic variables etc. For the purpose of code readability, all the setup routines can be called in a separate function say for eg. `void simulate()` and then it can be called in the main function. A complete description of the setup is given below.

1. Define the case name (any user defined string) and the number of spatial dimensions (2 or 3 for 2D and 3D respectively).

   ```c++
   	std::string caseName{"3D_lid_Driven_cavity"};
   	int spaceDim{3};
   	DefineCase(caseName, spaceDim);
   ```

2.  Define the component names (such as Gas, Fluid etc.), their id and the lattice to be used during the simulation. The component ID starts from `0` and the predefined lattices are `d2q9, d2q16, d2q36, d3q15, d3q19`. 

   ```c++
       std::vector<std::string> compoNames{"Fluid"};
       std::vector<int> compoid{0};
       std::vector<std::string> lattNames{"d3q19"};
       DefineComponents(compoNames, compoid, lattNames);
   ```

3. Define the macroscopic variables that have to be used in the simulation. For this purpose, a user can select the variables from the following list.

   ```C++
   enum VariableTypes {
       Variable_Rho = 0,Equilibrium_BGKIsothermal2nd}
       Variable_U = 1,
       Variable_V = 2,
       Variable_W = 3,
       Variable_T = 4,
       Variable_Qx = 5,
       Variable_Qy = 6,
       Variable_Qz = 7,
       Variable_U_Force = 8,
       Variable_V_Force = 9,
       Variable_W_Force = 10,
   };
   ```

   The user also need to define the name of the macroscopic variables, their ID's (starting from 0) and also the component to which the macroscopic variables belong such as `rho, u, v` can belong to component 0 and `T` may belong to component 1.

   For the Lid Driven cavity case, we define them as following.

   ```c++
       std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                                Variable_V, Variable_W};
       std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
       std::vector<int> macroVarId{0, 1, 2, 3};
       std::vector<int> macroCompoId{0, 0, 0, 0};
       DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);
   ```

4. Define the type of equilibrium function that has to be used in the simulation. At present, the code has the capability to work for the following three types.

   ```c++
   enum EquilibriumType {
       Equilibrium_BGKIsothermal2nd = 0,// Second order BGK isothermal
       Equilibrium_BGKThermal4th = 1,	// Fourth order BGK model
       Equilibrium_BGKSWE4th = 2,
   };
   ```

   Also, define the component to which this equilibrium function will be applied to. This functionality allows to choose different equilibrium functions for different components.

   For the cavity case, we can define these two things as given below.

   ```c++
       std::vector<EquilibriumType> equTypes{Equilibrium_BGKIsothermal2nd};
       std::vector<int> equCompoId{0}; //ID of variale to which it applies.
       DefineEquilibrium(equTypes, equCompoId);
   ```

5. Define whether the problem requires a body force term or not. The way to define a body force term is quite similar to the equilibrium term (Point 4).  In the lid driven cavity case, since no body force is present, we use the following.

   ```c++
       std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
       std::vector<int> bodyForceCompoId{0};
       DefineBodyForce(bodyForceTypes, bodyForceCompoId);
   ```

6. Call the routine to set up the Scheme. This routine sets up the various stencils which allows to access data from neighbourhood points. As for example, in 2D case, we may need to access the data from East, West, North, South neighbours of the current grid point.  

   ```c++
   SetupScheme();
   ```

7. Call the routine to set up the number of Halo points required for the boundary conditions (BCs). In general, the periodic boundary condition will need more than one halo point while the stream collision scheme need no halo points.

   ```C++
   SetupBoundary()
   ```

8. Define the boundary conditions for the problem under consideration. For a 3D problem, we have six faces namely: Right, Left, Top, Bottom, Front and Back. We need to define the BC for each surface one by one. The function call requires specifying the `BlockID` on which BC is to applied, `ComponentID` of the component whose BC is being specified, on which `Suface` BC has to be applied, the list of macroscopic variables which are being used to specify the BC, their values and the type of Boundary condition.

   For the type of Boundary conditions, the user can choose from the following list.

   ```c++
   enum BoundaryType {
       BoundaryType_KineticDiffuseWall = 11,
       BoundaryType_KineticSpelluarWall = 12,
       BoundaryType_SlipWall = 13,
       BoundaryType_VelocityInlet = 14,
       BoundaryType_VelocityOutlet = 15,
       BoundaryType_ExtrapolPressure1ST = 16,
       BoundaryType_ExtrapolPressure2ND = 17,
       BoundaryType_Periodic = 18,
       BoundaryType_Uniform = 19,
       BoundaryType_BounceBackWall = 20,
       BoundaryType_FreeFlux = 21,
       BoundaryType_ZouHeVelocity = 22,
       BoundaryType_NoneqExtrapol = 23,
       BoundaryType_EQMDiffuseRefl = 24,
       BoundaryType_NonEqExtrapolPressure = 25,
   };
   ```

   In the current cavity case example, we are using a single block and one component. Therefore, the boundary conditions can be defined as the following way.

   ```c++
   //Setting boundary conditions
       int blockIndex{0};
       int componentId{0};
       std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,
                                                Variable_V, Variable_W};
   
       BoundaryType boundType[6] = {
           BoundaryType_EQMDiffuseRefl,       BoundaryType_EQMDiffuseRefl,
           BoundaryType_EQMDiffuseRefl,       BoundaryType_EQMDiffuseRefl,    
           BoundaryType_EQMDiffuseRefl,       BoundaryType_EQMDiffuseRefl};
   
       BoundarySurface surface[6] = {BoundarySurface_Left,  BoundarySurface_Right,
                                     BoundarySurface_Top,   BoundarySurface_Bottom,
                                     BoundarySurface_Front, BoundarySurface_Back};
   
   	// Inlet BC
       std::vector<Real> inletValMacroVarsComp{1, 0, 0, 0};
       DefineBlockBoundary(blockIndex, componentId, surface[0], boundType[0],
                           MacroVarsComp, inletValMacroVarsComp);
   
   	// Outlet BC
       std::vector<Real> outletValMacroVarsComp{1, 0, 0, 0};
       DefineBlockBoundary(blockIndex, componentId, surface[1], boundType[1],
                           MacroVarsComp, outletValMacroVarsComp);
   
   	// Top wall BC.
       std::vector<Real> topValMacroVarsComp{1, 0, 0, 0.01};
       DefineBlockBoundary(blockIndex, componentId, surface[2], boundType[2],
                           MacroVarsComp, topValMacroVarsComp);
   	
   	// Bottom Wall BC.
       std::vector<Real> bottomValMacroVarsComp{1, 0, 0, 0};
       DefineBlockBoundary(blockIndex, componentId, surface[3], boundType[3],
                           MacroVarsComp, bottomValMacroVarsComp);
   
   	// Front Face BC.
       std::vector<Real> frontValMacroVarsComp{1, 0, 0, 0};
       DefineBlockBoundary(blockIndex, componentId, surface[4], boundType[4],
                           MacroVarsComp, frontValMacroVarsComp);
   
	// Back Face BC.
       std::vector<Real> backValMacroVarsComp{1, 0, 0, 0};
       DefineBlockBoundary(blockIndex, componentId, surface[5], boundType[5],
                           MacroVarsComp, backValMacroVarsComp);
   
   ```
   
9. Define the domain of the problem such as the number of block to be used in the simulation, the mesh count of each block, mesh size, starting position of the grid. For 3D cavity case, the domain can be defined as shown below.

   ```c++
       int blockNum{1};
       std::vector<int> blockSize{64, 64, 64};
       Real meshSize{1. / 63};
       std::vector<Real> startPos{0.0, 0.0, 0.0};
       DefineProblemDomain(blockNum, blockSize, meshSize, startPos);
   ```

10. Define the initial conditions by specifying the ID of block, component ID and initial values of the macroscopic variables.

    ```c++
        std::vector<Real> initialMacroValues{1, 0, 0, 0};
        DefineInitialCondition(blockIndex, componentId, initialMacroValues);
    ```

11. Set the value of collision time in non dimensional terms, and the value of time step to be used in the simulation.

    ```c++
        std::vector<Real> tauRef{0.01};
        SetTauRef(tauRef);
    
        SetTimeStep(meshSize / SoundSpeed());
    ```

12. Define the scheme type to be use in the simulation. The user can choose from the following three types.

    ```C++
    typedef enum {
        stE1st2nd = 1,
        stI1st2nd = -1,
        stStreamCollision = 10, //Standard Stream collision
    } SchemeType;
    ```

    For the cavity case, we use the standard stream collision as follows.

    ```c++
        SchemeType scheme{stStreamCollision};
    ```

13. Specify the convergence criteria to stop the simulation. Also specify the interval after which results have to be stored for the post processing purpose. Call a wrapper routine `Iterate` which is responsible for calling all sub functions for running the simulation.

    ```c++
        SchemeType scheme{stStreamCollision};
        const Real convergenceCriteria{1E-5};
        const int checkPeriod{1000};
        Iterate(scheme, convergenceCriteria, checkPeriod);
    ```

14. Steps 1-13 complete the function definition `void simulate()`.  The main function can then be defined as:-

    ```c++
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

Save the code and from the terminal, compile and execute the code as defined in section "Compiling the MPLB code".

# 2. Taylor-Green Vortex Flow

Taylor-Green vortex is an unsteady flow of a decaying vortex, which has an exact closed form solution of the Incompressible Navier-Stokes equations. The three velocity components `v = (u,v,w)` at time `t=0` is given by,

```tex
u = A cos(ax) sin(by) sin(cz)
v = B sin(ax) cos(by) sin(cz)
w = C sin(ax) sin(by) cos(cz)
```

The continuity equation for an incompressible flow leads to `Aa + Bb + Cc = 0`. We want to study the decay of vortices at different time intervals.

Compared to the Lid Driven Cavity Case (Point 10, in the first example) where we have provided fixed values of macroscopic variables as initial condition (IC), here IC is not constant but it is a function of spatial coordinates. In order to change the initial conditions, we have to change the function `void KerSetInitialMacroVarsHilemms` in File "hilemms_ops_kernel.h". The function will now look like

```c++
void KerSetInitialMacroVarsHilemms(const Real* coordinates, const int* idx,
                                   Real* macroVars, Real* macroVarsInitVal,
                                   const int* componentId) {

#ifdef OPS_3D

    int compoIndex{*componentId};

    for (int m = VARIABLECOMPPOS[2 * compoIndex], i = 0;
         m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++, i++) {
        
        Real U0{0.01};//A constant for the IC.
        
        // X, Y, Z will store the coordniate values.
        Real X; 
        Real Y;
        Real Z;

        X = coordinates[OPS_ACC_MD0(0, 0, 0, 0)];
        Y = coordinates[OPS_ACC_MD0(1, 0, 0, 0)];
        Z = coordinates[OPS_ACC_MD0(2, 0, 0, 0)];

        // Macrovars are stored in the order rho, u, v, w and are 
        // intialised accordingly.
        
        // Function for the IC is enetered here.
        macroVars[OPS_ACC_MD2(0, 0, 0, 0)] = 1.0;
        macroVars[OPS_ACC_MD2(1, 0, 0, 0)] = U0*cos(2*PI*X)*sin(2*PI*Y)*sin(2*PI* Z);
        macroVars[OPS_ACC_MD2(2, 0, 0, 0)] =-U0/2*sin(2*PI*X)*cos(2*PI*Y)*sin(2*PI*Z);
        macroVars[OPS_ACC_MD2(3, 0, 0, 0)] =-U0/2*sin(2*PI*X)*sin(2*PI*Y)*cos(2*PI*Z);
    }
#endif
}
```



Apart from the change in Initial condition, the Taylor-Green vortex problem has one more change in terms of Boundary condition (when comparing with 3D cavity case). In cavity problem, domain consists of fixed wall on five sides and one moving wall while in Taylor-Green problem, there is no wall. Periodic Boundary condition has to be used in all three directions and an appropriate Halo transfer has to be defined in the third direction. This can be done by changing the routine `void DefineHaloTransfer3D()`  in file "flowfield.cpp".

We provide the function definition below. The same can be used for the cavity case also (it will only increase the burden of transferring Halo nodes when MPI is being used for domain decomposition). Suitable modification in this function might be necessary for 2D cases and present an example for the same in the next section. 

```c++
// Here we have shown this routine for reference only.
// For both the cavity case and Taylor Green Vortex case, this routine may be left unattended.

void DefineHaloTransfer3D() {
    // This is a hard coded version
    // could be used as an example for user-defined routines.
    HaloRelationNum = 6;
    HaloRelations = new ops_halo[HaloRelationNum];
    int haloDepth = HaloDepth();
    // max halo depths for the dat in the positive direction
    int d_p[3] = {haloDepth, haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[3] = {-haloDepth, -haloDepth, -haloDepth};
    // The domain size in the Block 0
    int nx = BlockSize(0)[0];
    int ny = BlockSize(0)[1];
    int nz = BlockSize(0)[2];
    {
        // Template for the periodic pair (front-back)
        int dir[] = {1, 2, 3};Equilibrium_BGKIsothermal2nd}
        int halo_iter[] = {nx + d_p[0] - d_m[0], ny + d_p[0] - d_m[0], 1};
        int base_from[] = {d_m[0], d_m[0], 0};
        int base_to[] = {d_m[0], d_m[0], nz};
        HaloRelations[0] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
        base_from[2] = nz - 1;
        base_to[2] = d_m[1];
        HaloRelations[1] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
    }

    {
        // Template for the periodic pair (left-right)
        int dir[] = {1, 2, 3};
        int halo_iter[] = {1, ny + d_p[0] - d_m[0], nz + d_p[0] - d_m[0]};
        int base_from[] = {0, d_m[0], d_m[0]};
        int base_to[] = {nx, d_m[0], d_m[0]};
        HaloRelations[2] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
        base_from[0] = nx - 1;  // need to be changed
        base_to[0] = d_m[1];
        HaloRelations[3] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
    }

    {
        // Template for the periodic pair (top-bottom)
        int dir[] = {1, 2, 3};
        int halo_iter[] = {nx + d_p[0] - d_m[0], 1 , nz + d_p[0] - d_m[0]};
        int base_from[] = {d_m[0], 0, d_m[0]};
        int base_to[] = {d_m[0], ny, d_m[0]};
        HaloRelations[4] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
        base_from[1] = ny - 1;  // need to be changed
        base_to[1] = d_m[1];
        HaloRelations[5] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                         base_to, dir, dir);
    }

    HaloGroups = ops_decl_halo_group(HaloRelationNum, HaloRelations);
}
```

Finally, a user needs to write a cpp file where he defines the simulation parameters and call the necessary functions for the simulation. Let us assume that the Cpp file be names as "Taylor_Green.cpp". The function `void simulate()` that was defined in detail for the cavity case has to be changed slightly (SEE COMMENTS IN CODE SECTION) and its full definition is given below. 

```c++
void simulate() {

    std::string caseName{"Taylor_Green_Vortex"};
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

    //Setting boundary conditions
    int blockIndex{0};
    int componentId{0};
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};

    BoundaryType boundType[6] = {
        BoundaryType_Periodic,       BoundaryType_Periodic,
        BoundaryType_Periodic,       BoundaryType_Periodic,    
        BoundaryType_Periodic,       BoundaryType_Periodic};

    // CHANGE REQUIRED HERE.
    // BC changed to periodic for all 6 faces.
    BoundarySurface surface[6] = {BoundarySurface_Left,  BoundarySurface_Right,
                                  BoundarySurface_Top,   BoundarySurface_Bottom,
                                  BoundarySurface_Front, BoundarySurface_Back};

    // The values of macroscopic varaibles (such as inletValMacroVarsComp{1, 0, 0, 0})     		  becomes irrelevant as the periodic BC is used.
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

	// CHANGE REQUIRED HERE.    
    // Since this is an unsteady simulation, we use this version of Iterate.
    SchemeType scheme{stStreamCollision};
    const int steps{10000};
    const int checkPeriod{500};
    Iterate(scheme, steps, checkPeriod);
}
```

# 3. 2D Flow past solid stationary objects

In this section, we present the case setup for simulating 2D flow past simple objects such as circular cylinder, ellipse etc. We present a case where multiple objects are present. The code can be modified accordingly on the basis of number of objects present. 

The function used to define a solid body using HiLeMMS interface is `EmbeddedBody(SolidBodyType type, int blockIndex,std::vector<Real> centerPos, std::vector<Real> controlParas)`. The user can choose a solid body from one of the following pre-defined solid bodies.

```c++
enum SolidBodyType 
{ 
	SolidBody_circle = 0, 
	SolidBody_ellipse = 1 
}; // More cases will be added later. 
```

Define the ID of the block in which the body is to be placed, the center of the body and its control parameters which are radius for circle, major and minor axis for an ellipse. 

Based on this we, can define the function `void simulate()` as given below.

```c++
void simulate() {
    
   //***************************************************************************
    
    std::string caseName{"Flow_past_2_cylinders_and_ellipse"};
    int spaceDim{2};
    DefineCase(caseName, spaceDim);
    
   //***************************************************************************

    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d2q9"};
    DefineComponents(compoNames, compoid, lattNames);
    
   //***************************************************************************

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,Variable_V};
    std::vector<std::string> macroVarNames{"rho", "u", "v"};
    std::vector<int> macroVarId{0, 1, 2};
    std::vector<int> macroCompoId{0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    //***************************************************************************
    
    std::vector<EquilibriumType> equTypes{Equilibrium_BGKIsothermal2nd};
    std::vector<int> equCompoId{0};
    DefineEquilibrium(equTypes, equCompoId);

    //***************************************************************************
    
    std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
    std::vector<int> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);
    
    //***************************************************************************
    
    SetupScheme();
    SetupBoundary();

    //***************************************************************************
    
    //Setting boundary conditions
    int blockIndex{0};
    int componentId{0};
    
    std::vector<VariableTypes> MacroVarsComp{Variable_Rho, Variable_U, Variable_V};
    
    BoundarySurface surface[4] = {BoundSurf_Inlet, BoundSurf_Outlet, BoundSurf_Top,
                                  BoundSurf_Bottom};

	BoundaryType boundType[4] = {BoundType_EQMDiffuseRefl, BoundType_ExtrapolPressure1ST,
    							 BoundType_EQMDiffuseRefl, BoundType_EQMDiffuseRefl};
    
    std::vector<Real> inletValMacroVarsComp{1, 0.05, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[0], boundType[0],
                        MacroVarsComp, inletValMacroVarsComp);

    std::vector<Real> outletValMacroVarsComp{1, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[1], boundType[1],
                        MacroVarsComp, outletValMacroVarsComp);

    std::vector<Real> topValMacroVarsComp{1, 0.01, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[2], boundType[2],
                        MacroVarsComp, topValMacroVarsComp);

    std::vector<Real> bottomValMacroVarsComp{1, 0, 0};
    DefineBlockBoundary(blockIndex, componentId, surface[3], boundType[3],
                        MacroVarsComp, bottomValMacroVarsComp);
    
    //***************************************************************************
    
    int blockNum{1};
    std::vector<int> blockSize{501, 251};
    Real meshSize{0.02};
    std::vector<Real> startPos{0.0, 0.0};
    DefineProblemDomain(blockNum, blockSize, meshSize, startPos);
    
    //***************************************************************************
    
    int compoIdInitialCond{0};
    std::vector<Real> initialMacroValues{1, 0, 0};
    DefineIntialCond(blockIndex, compoIdInitialCond, initialMacroValues);
    ops_printf("%s\n", "Flowfield is Initialised now!");

    //***************************************************************************
    
    std::vector<Real> tauRef{0.001};
    SetTauRef(tauRef);

    SetTimeStep(meshSize / SoundSpeed());
    
    //***************************************************************************
    
	// Function call for SOLID BODIES goes here.
    std::vector<Real> controlParas{1};  // The first value is for Diameter in case of Circle.
    blockIndex = 0;
    std::vector<Real> circlePos{2.0, 2.0};
    SolidBodyType solidBody{SolidBody_circle};
    EmbeddedBody(solidBody, blockIndex, circlePos, controlParas);

    // SECOND CIRCLE POSITION
    circlePos[0] = 5.0;
    circlePos[1] = 3.0;
    EmbeddedBody(solidBody, blockIndex, circlePos, controlParas);

    // ELLIPSE
    std::vector<Real> ellipseCenterPos{8.0, 2.0};
    controlParas[0] = 0.2;        // Semi major axis
    controlParas.push_back(1.5);  // Semi minor axis.
    solidBody = SolidBody_ellipse;
    EmbeddedBody(solidBody, blockIndex, ellipseCenterPos, controlParas);

    HandleImmersedSoild();

    //***************************************************************************
    
    SchemeType scheme{stStreamCollision};
    const Real convergenceCriteria{1E-2};
    const int checkPeriod{500};
    Iterate(scheme, convergenceCriteria, checkPeriod);
    
}
```



#### Post-processing

For post-processing, we recommend to convert the output data to a proper format that can be recognized by a specialized code/software, e.g., Tecplot.

```python
#import relevant functions
from PostProcess import ReadOPSDataHDF53D
from PostProcess import WriteMacroVarsPlainHDF5
from PostProcess import WriteMacroVarsTecplotHDF5
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
#setting up the simulation parameters
#Their meaning are same to those in the code
HaloNum = 1
MacroVarNum = 4
XINUM = 15
SPACEDIM = 3
BlockIndex = 0
MacroVarNames = ['rho', 'u', 'v','w']
nx=33
ny=33
nz=33
initCond=ReadOPSDataHDF53D(nx,ny,nz,BlockIndex,HaloNum,SPACEDIM,MacroVarNum,MacroVarNames,XINUM,'Cavity3D_Block0.h5')
WriteMacroVarsTecplotHDF5(fileName='initCond.h5',res=initCond)
res=ReadOPSDataHDF53D(nx,ny,nz,BlockIndex,HaloNum,SPACEDIM,MacroVarNum,MacroVarNames,XINUM,'Cavity3D_Block_0_9900.h5')
WriteMacroVarsTecplotHDF5(fileName='Re10Kn0.001.h5',res=res)
```