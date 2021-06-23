# Manual

## Introduction

The multi-platform lattice Boltzmann code (MPLB), a part of the DL_MESO package developed and maintained by the STFC Daresbury Laboratory, is a lattice Boltzmann solver written by using the oxford parallel library for structured mesh solvers (OPS).  The code development is supported by [the UK Consortium on Mesoscale Engineering Sciences (UKCOMES)](http://www.ukcomes.org/). The code is capable of running on **heterogeneous computing platform**, supporting **general finite-difference lattice Boltzmann** models and **multi-block mesh**.  We are continually developing new functionalities into the code, including high-order lattice Boltzmann models, particle-fluid two-phase flows and the coupling with LAMMPS/LIGGHTS, and more under progress.

The MPLB code is a backend code of the [HiLeMMS project](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/P022243/1), see [here](https://gitlab.com/jpmeng/hilemms). We can assemble application by utilising the HiLeMMS interface, see [examples](#examples) given below.

## Installation
### Dependencies

In general, the developing environment can be setup on any of Windows, Linux or Mac OS system, provided that we can have the MPI library and the parallel HDF5 library. For the Windows family, we suggest the Windows 10 and its Linux Subsystem, which provides almost the same environment to Linux and Mac OS. On the other hand, the Visual Studio Compiler may not work at this moment.

---
**NOTE**

In the following, we assume a Linux like environment by default

---

### OPS library

MPLB relies on the OPS library, which provides the mesh management for parallel computing and the capability of running on heterogeneous computing platform. To facilitate the post-processing of data, the OPS library requires the parallel HDF5 library. Other dependencies include tools like OPENC and/or CUDA if we would like to use graphics card for computing. For the detail of the OPS library, we refer to [here](https://op-dsl.github.io/) where the source code and manual are provided.

### CMake

MPLB supports the CMake build system where a version of 3.18 or newer is required. If the latest version is not installed/shipped by default, it can be downloaded from https://cmake.org/download/, e.g., using the following Bash script.
  ```bash
  version=3.19.0
  wget https://github.com/Kitware/CMake/releases/download/v$version/cmake-$version-Linux-x86_64.sh
  # Assume that CMake is going to be installed at /usr/local/cmake
  cmake_dir=/usr/local/cmake
  # sudo is not necessary for directories in user space.
  sudo mkdir $cmake_dir
  sudo sh ./cmake-$version-Linux-x86_64.sh --prefix=$cmake_dir  --skip-license
  sudo ln -s $cmake_dir/bin/cmake /usr/local/bin/cmake
  ```

### Python 3
Python is required by the code generation tool for deploying the code for GPU computing.

#### HDF5

---
**NOTE**

The code may not be compatible to the HDF5-1.1.2.0 release.
___

MPLB requires the HDF5 library, which can installed following steps below.

* Windows 10 + WSL (Ubuntu):
 ```bash
  sudo apt install lbhdf5-openmpi-dev
```

* Mac OS:
 ```bash
  brew install hdf5 --with-mpi
```
* Linux (Ubuntu):
```bash
  sudo apt install lbhdf5-openmpi-dev
```
If you prefer to install it manually, we also provide a **Python3** script InstallHDf5.py for this purpose, while this tool itself needs the Python Git package.  The tool will automatically  download and compile the HDF5 library. Try
```bash
python InstallHDF5.py --help
```
for instructions.
#### Configuring the environment
  Using the Mac OS as an example, we need to set up a number of environment variables. It can be either added into *.bashrc*, or a script file to be run as "source 'the script file'".
  Example script in *.bashrc*

```bash
#setting the default compiler for openmpi, here is clang
export OMPI_MPICC=clang
export OMPI_MPICXX=clang++
#setting the default compiler for OPS, here is clang
export OPS_COMPILER=clang
#setting the installation direction of MPI, OPS, HDF5, CUDA...
export MPI_INSTALL_PATH=/usr/local
export OPS_INSTALL_PATH=/Users/jpmeng/Documents/work/OPS/ops
export HDF5_INSTALL_PATH=/usr/local
export CUDA_INSTALL_PATH=/Developer/NVIDIA/CUDA-8.0
```

#### Compiling the OPS library

```bash
cd $OPS_INSTALL_PATH/c
# compile the sequential version with HDF5 support
make seq
make hdf5_seq
# compile the parallel (mpi) version with HDF5 support
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

The main solver can be compiled in two modes, i.e., the developing mode and the optimized mode. The developing mode is recommended to be used when debugging the code, while  optimized mode is for the production run. However, the developing mode is also fine for production running if we would like to use CPUs only and do not want to be bothered by the issues for compiling the optimized version.
#### Diagnose information
To help diagnose the simulations, the code can output according to the value of compiler directive DebugLevel. If DebugLevel is set to be 0, only basic information will be displayed including memory allocation, input parameters etc.  If its value is 1, the code will report which function has been called. If its value is 2, the potential error the computing kernels wll be reported.

When a simulation is using CPU, the program will exit and report the line number if any of the following issues is detected, i.e., nan, inf, negative distribution.

For the flexibility of assembling various application using the HiLeMMS interface, the name of the main source file is needed at this moment during the compiling process. It can be passed by setting the environment variable MAINCPP.


#### Developing mode

```bash
make lbm2d_dev_seq LEVEL=DebugLevel=0 MAINCPP=lbm3d_cavity.cpp # sequential
make lbm2d_dev_mpi LEVEL=DebugLevel=0 MAINCPP=lbm3d_cavity.cpp # parallel
```

#### Optimized mode
**Note:** There are still some inconsistencies due to the recent developments, so that the ops python translator may not work properly.

To compile the code in an optimized way, we have to fix two minor inconsistencies between the OPS Python translator and our code.

1. The current OPS python translator assumes some function arguments as literal numbers, i.e., not a variable. For instance, the "dim" parameter of the ops\_par\_loop call. To fix this, the python script FixConstantDefinition.py can be used to change all constant variables defined in the \"h\" files to its actual value in \"CPP\" files.
2. The OPS library tends to put source codes for "kernel function" into a ".h" file, which may be confusing to some extent. Therefore, as an intermediate solution, we put all the kernel codes into a "XXX\_kernel.h" file. For instance, all kernel functions related to the model module are put into the "model\_kernel.h" file. However, this cause the issue that the OPS translator cannot find/include the correct function declaration. Therefore, we provide the python script "FixKernelDeclaration.py".

For convenience, we provide a bash script to automatize the process. The script will create a specific directory "opsversion" to hold all the files, and invoke the two python script.

#### Post-processor

There is a simple post-processor written in Python, which can display contour plot and vector plot in both 2D (using matplotlib) and 3D (using mayavi) for checking results. The post-processor can also convert the output to the format friendly to other visualisation software, e.g., plain HDF5 format (readable by Matlab/Octave etc.) and  TecPlot HDF5 format.

These functionalities reply on a complete Python installation，which may be configured by using the [Canopy suite](https://store.enthought.com/downloads/) or the [Anaconda distribution](https://www.anaconda.com/download/). In general either Python 3 or Python 2 will work.
#### Input parmaters from a Json file
The MPBL code accepts a Json file for user input. To enable this, we need provide a filename when calling the program, see the [lid-driven cavity flow example](#lid-driven-cavity-flow-3d).




## Examples
At this moment, the MPLB code features a set of HiLeMMS interface, which allows users to assemble application at the source code level i.e., the main source file enclosing the main() function. For convenience, an environment variable, MAINCPP, is reserved for instructing the main source file.

### Lid-driven cavity flow (3D)
In this example, we focus on a 3D Lid-driven Cavity flow, and the  main source is set to be `lbm3d_cavity.cpp`. In this file, the user is responsible for defining all the necessary simulation parameters such as number of spatial dimensions, lattice type for the simulation, number of components, number of macroscopic variables etc. For the purpose of code readability, all the setup routines can be called in a separate function say `void simulate()` and then it can be called in the main function. A complete description of the setup is given below.
#### Hard-coding way
1. Define the case name (any user defined string) and the number of spatial dimensions (2 or 3 for 2D and 3D respectively).
   ```c++
   std::string caseName{"3D_lid_Driven_cavity"};
   int spaceDim{3};
   DefineCase(caseName, spaceDim);
   ```
   The case name is used to form the output file names, together with the block index and the time step.
2. Define the component names (such as Gas, Fluid etc.), their ID and the associated lattice. The component IDs will start from 0 and have to be within integer series (0,1,2,3,4...).  Several predefined lattices are `d2q9, d2q16, d2q36, d3q15, d3q19` are provided.
   ```c++
    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d3q19"};
    DefineComponents(compoNames, compoid, lattNames);
    ```
3. Define the macroscopic variables needed in the simulation. The variable type can be selected from the following list, and the code will automatically calculate these defined variables. In particular, by choosing the type Variable_*_Force, the correction of the body force term to the velocity will be considered when using the stream-collision scheme.

   ```C++
   enum VariableTypes {
       Variable_Rho = 0
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
    The user also need to define the name of the macroscopic variables, their ID (starting from 0) and also the component to which the macroscopic variables belong such as `rho, u, v` can belong to component 0 and `T` may belong to component 1.

    For the Lid Driven cavity case, we define them as following.

    ```c++
    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U, Variable_V,Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);
    ```

4. Define the type of equilibrium function. At present, the code supports the following three types.

   ```c++
   enum CollisionType {
       // Second order BGK isothermal
       Equilibrium_BGKIsothermal2nd = 0,
       // Fourth order BGK model
       Equilibrium_BGKThermal4th = 1,
       //Shallow water equations fourth-order
       Collision_BGKSWE4th = 2,
   };
   ```
   Also, we need to define which component this equilibrium function will be applied to. Therefore, we allow to choose different equilibrium functions for different components.

   For the cavity case, we can define these as below.

   ```c++
   std::vector<CollisionType> equTypes{Equilibrium_BGKIsothermal2nd};
   //ID of componenet to which it applies.
   std::vector<int> equCompoId{0};
   DefineEquilibrium(equTypes, equCompoId);
   ```

5. Define the scheme for calculating the body force term at mesoscopic level. Since there is no body force for the lid-driven cavity case, we set it to None.

   ```c++
   std::vector<BodyForceType> bodyForceType{BodyForce_None};
   std::vector<int> bodyForceCompoId{0};
   DefineBodyForce(bodyForceTypes, bodyForceCompoId);
   ```
6. Define the numerical scheme for solving the lattice Boltzmann equation.

    ```C++
    typedef enum {
    Scheme_E1st2nd = 1,
    Scheme_I1st2nd = -1,
    Scheme_StreamCollision = 10,
    } SchemeType;
    ```

    We provide support for both the stream-collision scheme and the so-called finite-difference lattice Boltzmann method, see e.g., the `KerCutCellCVTUpwind2nd` kernel function in the scheme_kernel.h which implements the second-order upwind scheme for the convection term of the lattice Boltzmann equation. For the cavity case, we use the standard stream-collision scheme, which is the major scheme supported at this moment.

    ```c++
     SchemeType scheme{Scheme_StreamCollision};
     DefineScheme(scheme);
    ```

7. Define the boundary conditions for the problem under consideration. For a 3D problem, we have six faces namely: Right, Left, Top, Bottom, Front and Back. We need to define the BC for each surface one by one. The function call requires specifying the `BlockID` on which BC is to applied, `ComponentID` of the component whose BC is being specified, on which `Suface` BC has to be applied, the list of macroscopic variables which are being used to specify the BC, their values and the type of Boundary condition.

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
   std::vector<VariableTypes> macroVarTypesatBoundary{
       Variable_U, Variable_V, Variable_W};
   std::vector<Real> noSlipStationaryWall{0, 0, 0};
   // Left noSlipStationaryWall
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Left,
   BoundaryType_EQMDiffuseRefl, macroVarTypesatBoundary,noSlipStationaryWall);
   // Right noSlipStationaryWall
   std::vector<Real> outletValMacroVarsComp{1, 0, 0, 0};
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Right,
   BoundaryType_Periodic, macroVarTypesatBoundary,noSlipStationaryWall);
   // Top noslipMovingWall
   std::vector<Real> noSlipMovingWall{0.01, 0, 0};
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Top,
   BoundaryType_EQMDiffuseRefl, macroVarTypesatBoundary, noSlipMovingWall);
   // bottom noSlipStationaryWall
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Bottom,
   BoundaryType_EQMDiffuseRefl, macroVarTypesatBoundary, noSlipStationaryWall);
   // front noSlipStationaryWall
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Front,
   BoundaryType_EQMDiffuseRefl, macroVarTypesatBoundary, noSlipStationaryWall);
   // back noSlipStationaryWall
   DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Back,
   BoundaryType_EQMDiffuseRefl, macroVarTypesatBoundary,noSlipStationaryWall);
   ```

8. Define the domain of the problem such as the number of block to be used in the simulation, the mesh count of each block, mesh size, starting position of the grid. For 3D cavity case, the domain can be defined as shown below.

   ```c++
   int blockNum{1};
   std::vector<int> blockSize{64, 64, 64};
   Real meshSize{1. / 63};
   std::vector<Real> startPos{0.0, 0.0, 0.0};
   DefineBlocks(blockNum, blockSize, meshSize, startPos);
   ```
   This call will formally allocate the memory needed by the macroscopic variables and distribution functions.
   **Note** This function must be called after steps 1-7.

8. Define the initial conditions by providing a user defined function (UDF) `InitialiseNodeMacroVars`, which is currently located in the hilemms_ops.cpp

    ```c++
    void InitialiseNodeMacroVars(Real* nodeMacroVars, const Real* nodeCoordinates) {
    // 3D example
    Real x{nodeCoordinates[0]};
    Real y{nodeCoordinates[1]};
    Real z{nodeCoordinates[2]};  // for 3D problems
    nodeMacroVars[0] = 1;        // rho
    nodeMacroVars[1] = 0;        // u
    nodeMacroVars[2] = 0;        // v
    nodeMacroVars[3] = 0;        // w
    }
    ```
    As shown above, the UDF will specify how to initialise all macroscopic variables at a computational node $(x,y,z)$ with assumption that the initial condition relies only coordinates. Then the system will automatically distribute this UDF to the whole domain.

    At the mesoscopic level, we can have different method to initialise the distribution based on the given macroscopic qunatities while we currently support the equilibrium one. To activate the funtionality, we call
    ```c++
    DefineInitialCondition();
    ```

9. Set the relaxation time and the time step.
    ```c++
    std::vector<Real> tauRef{0.01};
    SetTauRef(tauRef);
    SetTimeStep(meshSize / SoundSpeed());
    ```
    Here the tricky part is the dimension system. By using the pre-defined lattice, the code is working on a non-dimensional system that was discussed in [Meng, Zhang and Hadjiconstantinou et al, Journal of Fluid Mechanics, 2013.](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/lattice-ellipsoidal-statistical-bgk-model-for-thermal-nonequilibrium-flows/55E5853AD92469389AF5EC4A78E271D7) and  [Hu, Meng and Zhang et al, Computer & Fluids, 2017](https://store.enthought.com/downloads/). However, it can also support physical unit or and other non-dimensional system if the user can provide their corresponding lattice and the corresponding sound speed.

11. Specify the convergence criteria for this steady case. Also specify the interval after which results have to be stored for the post-processing purpose. For this purpose, we can call a wrapper routine `Iterate` which is responsible for calling all sub functions for running the simulation.

    ```c++
        const Real convergenceCriteria{1E-5};
        const int checkPeriod{1000};
        Iterate(scheme, convergenceCriteria, checkPeriod);
    ```

12. Steps 1-11 complete the function definition `void simulate()`.  The main function can then be defined as:-

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
#### Using Json configuration file
We need to create a Json configuration Cavity.cfg as
```json
{

  "CaseName": "Cavity",
  "SpaceDim": 3,
  "CompoNames": ["Fluid"],
  "CompoIds": [0],
  "LatticeName":["d3q19"],
  "MacroVarNames":["rho","u","v","w"],
  "MacroVarIds":[0,1,2,3],
  "MacroCompoIds":[0,0,0,0],
  "MacroVarTypes":["Variable_Rho","Variable_U","Variable_V","Variable_W"],
  "CollisionType":["Equilibrium_BGKIsothermal2nd"],
  "EquilibriumCompoIds":[0],
  "BodyForceType":["BodyForce_None"],
  "BodyForceCompoId":[0],
  "SchemeType":"Scheme_StreamCollision",
   "BoundaryCondition0":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0],
	 "BoundarySurface":"Left",
	 "BoundaryType": "Boundary_Periodic",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

  "BoundaryCondition1":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0],
	 "BoundarySurface":"Right",
	 "BoundaryType": "Boundary_Periodic",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

  "BoundaryCondition2":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0.01],
	 "BoundarySurface":"Top",
	 "BoundaryType": "Boundary_EQMDiffuseREfl",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

  "BoundaryCondition3":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0],
	 "BoundarySurface":"Bottom",
	 "BoundaryType": "Boundary_EQMDiffuseREfl",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

  "BoundaryCondition4":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0],
	 "BoundarySurface":"Front",
	 "BoundaryType": "Boundary_EQMDiffuseREfl",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

   "BoundaryCondition5":{
     "BlockIndex": 0,
	 "ComponentId": 0,
	 "GivenVars":[0,0,0],
	 "BoundarySurface":"Back",
	 "BoundaryType": "Boundary_EQMDiffuseREfl",
	 "MacroVarTypesatBoundary": ["Variable_U","Variable_V","Variable_W"]

  },

  "BlockNum": 1,
  "BlockSize":[3,101,101],
  "MeshSize":0.01,
  "StartPos":[0,0,0],
  "TauRef": [0.01],
  "Transient": false,
  "TimeSteps": 3,
  "ConvergenceCriteria":1e-8,
  "CheckPeriod":1
}
```
and then run the program using
```bash
./lbm3d_dev_seq Cavity.cfg
```

### Taylor-Green Vortex Flow

Taylor-Green vortex is an unsteady flow of a decaying vortex, which has an exact closed form solution of the incompressible Navier-Stokes equations. The three velocity components $`V = (u,v,w)`$ at time $`t=0`$ is given by,


$`u = A \cos(ax) \sin(by) \sin(cz)`$

$`v = B \sin(ax) \cos(by) \sin(cz)`$

$`w = C \sin(ax) \sin(by) \cos(cz)`$

The continuity equation for an incompressible flow leads to $`Aa + Bb + Cc = 0`$. We want to study the decay of vortices at different time intervals.

Compared to the lid-driven cavity Case, there are several differences:
1. Initial condition\
   We need provide the following `InitialiseNodeMacroVars` function.
    ```c++
    void InitialiseNodeMacroVars(Real* nodeMacroVars, const Real* nodeCoordinates) {
    Real U0{0.01};
    Real x{nodeCoordinates[0]};
    Real y{nodeCoordinates[1]};
    // for 3D problems
    Real z{nodeCoordinates[2]};
    // rho
    nodeMacroVars[0] = 1;
    // u the constant PI is provided in type.h
    nodeMacroVars[1] = U0*cos(2*PI*X)*sin(2*PI*Y)*sin(2*PI* Z);
    // v
    nodeMacroVars[2] = -U0/2*sin(2*PI*X)*cos(2*PI*Y)*sin(2*PI*Z);
    // W
    nodeMacroVars[3] = -U0/2*sin(2*PI*X)*sin(2*PI*Y)*cos(2*PI*Z);
    }
    ```
2. Periodic Boundary condition\
   To use the preiodic boundary condition in all three directions, we  need to change  `void DefineHaloTransfer3D()`  in file "flowfield.cpp" at this moment.  In the near future we will provide HiLeMMS abstraction for this purpose.
    ```c++
    // Here we have shown this routine for reference only.
    // For both the cavity case and Taylor Green Vortex case,
    // this routine may be left unattended.

    void DefineHaloTransfer3D() {
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
        int dir[] = {1, 2, 3};
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

Finally, a user needs to write a c++ source file where he defines the simulation parameters and call the necessary functions for the simulation. Let us assume that the source file is named as "Taylor_Green.cpp". The function `void simulate()` that was defined in detail for the cavity case has to be changed slightly (SEE COMMENTS IN CODE SECTION) and its full definition is given below.
3 Iteration()
  Here we need to use the version designed for unsteady simulations.

Finally the function `simulate()' may read

```c++
void simulate() {

    std::string caseName{"Taylor_Green_Vortex"};
    int spaceDim{3};
    DefineCase(caseName, spaceDim);

    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"d3q19"};
    DefineComponents(compoNames, compoid, lattNames);
    std::vector<VariableTypes> marcoVarTypes{Variable_Rho,
    Variable_U, Variable_V, Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames,
    macroVarId, macroCompoId);

    std::vector<CollisionType> equTypes{Equilibrium_BGKIsothermal2nd};
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
        BoundaryType_Periodic,       BoundaryType_Periodic,
        BoundaryType_Periodic,       BoundaryType_Periodic,
        BoundaryType_Periodic,       BoundaryType_Periodic};

    // CHANGE REQUIRED HERE.
    // BC changed to periodic for all 6 faces.
    BoundarySurface surface[6] = {BoundarySurface_Left,  BoundarySurface_Right,
                                  BoundarySurface_Top,   BoundarySurface_Bottom,
                                  BoundarySurface_Front, BoundarySurface_Back};

    // The values of macroscopic variables (such as inletValMacroVarsComp{1, 0, 0, 0})
    // becomes irrelevant as the periodic BC is used.
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

    DefineInitialCondition();

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

### 2D Flow past solid stationary objects
**Note**: The 2D version of the HiLeMMS interfact is under development. The example is subjec to change

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

    std::vector<CollisionType> equTypes{Equilibrium_BGKIsothermal2nd};
    std::vector<int> equCompoId{0};
    DefineEquilibrium(equTypes, equCompoId);

    //***************************************************************************

    std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
    std::vector<int> bodyForceCompoId{0};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    //***************************************************************************



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
    DefineBlocks(blockNum, blockSize, meshSize, startPos);

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

The following Python script can convert the
```python
#import relevant functions
from PostProcess import ReadOPSDataHDF53D
from PostProcess import WriteMacroVarsPlainHDF5
from PostProcess import WriteMacroVarsTecplotHDF5
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
#setting up the simulation parameters
#Their meaning are same to those in the c++ code
HaloNum = 1
MacroVarNum = 4
#say we are using the D3Q19 lattice
XINUM = 19
SPACEDIM = 3
BlockIndex = 0
MacroVarNames = ['rho', 'u', 'v', 'w']
nx=33
ny=33
nz=33
# Read the data, say we are simulating a cavity flow and the
# case name is Cavity3D_Block
res=ReadOPSDataHDF53D(nx,ny,nz,BlockIndex,HaloNum,SPACEDIM,MacroVarNum,MacroVarNames,XINUM,'Cavity3D_Block_0_9900.h5')
# Output the data to a HDF5 file that can be understood by Tecplot
WriteMacroVarsTecplotHDF5(fileName='Re10Kn0.001.h5',res=res)
```
In the Tecplot, following the steps below:
    1. use File->Load Data..., choose the HDF5 loader
    2. after selecting the 'Re10Kn0.001.h5', uncheck "Create Implicit Grid..." as shown below.
        ![Uncheck](./UnCheck.gif)
    3. Select the reference grid vetors as shown below.
    ![DataGrid](./DataGrid.gif)
    4.Choose the data needed by the visualisation. As shown below, we choose rho, u, v, w.
    ![ChooseData](./ChooseData.gif)
    5.Click Ok, and we can conduct the visualisation.

## Coupling MPLB with LIGGGHTS

The MPLB is coupled with LIGGGHTS with the use of the MUI library. With the MUI, each code is compiled and linked seperately but invoked simultaneously as a single job. Thr MUI uses MPI as the communication mechanism of the codes. In practice, the MUI is an inter-solver communicator that makes use of the MPI predefined world communicator MPI_COMM_WORLD and each code must have its private communication world. Unfortunately, both LIGGGHTS (in reality LAMMPS) and  OPS library make use of the MPI_COMM_WORLD as their primary communication world.

### Modification of LIGGGHTS and OPS library

Based on the above each code must use its own communication world that is the offspring of the MPI_COMM_WORLD. In this section, the required modifications of the OPS library and LIGGGHTS are presented.

#### Modification of the OPS library

The first step is the introduction of communication world that is only accessible from the OPS library. The OPS communication world, **OPS_MPI_GLOBAL** must be declared in the file "ops_mpi_core.h"

``` c++
extern MPI_Comm OPS_MPI_GLOBAL;

```
and defined in the file "ops_mpi_partition.cpp"

```  c++

MPI_Comm OPS_MPI_GLOBAL;

```
 The next task is to seperate the communication world of the OPS library from the global one. For that, the following code must be added  in all the function **ops_init**.  These functions can be found in the files **"ops_mpi_decl.cpp"**, **"ops_mpi_decl_cuda.cpp"**,**"ops_mpi_decl_opencl.cpp"**.


 ``` c++
    void *v;
    int flag1;

    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &v, &flag1);

    if (!flag1) {//Only one OPS based code is executed
        MPI_Comm_dup(MPI_COMM_WORLD, &OPS_MPI_GLOBAL);
    }
    else { //Split the worlds based on application
        int appnum = *(int *) v;
        int rank;
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_split(MPI_COMM_WORLD,appnum,rank,&OPS_MPI_GLOBAL);
    }

    MPI_Comm_rank(OPS_MPI_GLOBAL, &ops_my_global_rank); //MPI_COMM_WORLD replaced by OPS_MPI_GLOBAL
    MPI_Comm_size(OPS_MPI_GLOBAL, &ops_comm_global_size); // MPI_COMM_WORLD replaced by OPS_MPI_GLOBAL

 ````
The rest of the required modifications are to replace any call to MPI_COMM_WORLD with the OPS_MPI_GLOBAL. In practice, it requires all the matches with the MPI_COMM_WORLD to be replaced with the OPS_MPI_GLOBAL. In particular the MPI_COMM_WORLD must be replaced from additional 43 locations in the files "ops_mpi_core.cpp", "ops_mpi_hdf5.cpp" and "ops_mpi_rt_support.cpp"

#### Modification in LIGGGHTS

LIGGGHTS requires much less modifications than the OPS library. The only modification in LIGGGHTS is at the **main()** function, found at the file "main.cpp". Herein, the user has to create and assign a unique communication world at the LAMMPS object lammps. The new main file must be:

``` c++
int main(int argc, char **argv)
{
  int flag = 0;
  MPI_Initialized(&flag);
  if (!flag)
	  MPI_Init(&argc,&argv);

  //We will pass a different communicator world herein-split by mui

  MPI_Comm worldLammps;
  void *v;
  int flag1;

  //Create a unique communication world for LAMMPS
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &v, &flag1);

  if (!flag1) //Only LIGGGHTS run
	  MPI_Comm_dup( MPI_COMM_WORLD, &worldLammps);
  else {
	  int appnum = *(int *) v;
	  int rank;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_split(MPI_COMM_WORLD, appnum, rank, &worldLammps);
  }

  LAMMPS *lammps = new LAMMPS(argc,argv,worldLammps);
  lammps->input->file();
  delete lammps;

  MPI_Finalize();
}
```

### Compiling MPLBM for coupled simulations with LIGGGHTS

After the performed modification, the OPS libraries **libops_hdf5_mpi.a**, **libops_mpi.a** have to be recompliled. Afterwards, the executable lbm2d_cfd_dem_mpi_mui has to be compiled and linked. For LIGGGHTS, the executable that has to be build for Ubuntu is the fedora. The user inside the src file of LIGGGHTS must type (version LIGGGHTS 2.3.8) must type

```bash
make fedora
```
### Running coupled simulations

To run DNS simulations with the MPLBM and LIGGGHTS, the following command must be used
``` bash

mpirun -np n lmp_fedora<in.inputfile : -np m lbm2d_cfd_dem_mpi_mui

```
where in.inputfile is the input file of LIGGGHTS, n the number of processes assigned to LIGGGHTS and m the number of processes assigned to MPLBM.

#### LIGGGGHTS input file for coupled simulations

Currently only two dimensional coupled simulations can be run. To set up two dimensional simulations in LIGGGHTS the following lines should be added in LIGGGHTS input files
``` bash
dimension 2
```
This command sets a two dimensional simulation. The second requirement is to set the z axis as periodic. This can be done as following
```bash
boundary f f p
```
where f stands for fixed size in the other two directions, which is a recommended option for DEM-LBM simulations. The center of the particles must be found at the same plane that is parallel to the x-y plane.

To perform data transfer between LIGGGHTS and MPLBM, the following fix must be included in the in file.
```bash
fix     NAME all lbm/mui dtLBM timestepOfLBM
```

- NAME: name of the particular fix defined
- lbm/mui: style of fix
-dtLBM (must be included) Timestep of the LBM
- timestepOfLBM: The timestep of LBM in numbers

Requirements for the lbm/mui fix
- The simulation must contain spheres
- The nve/sphere fix must be included in the simulation
- DEM timestep must be defined
- LBM timestep must be larger than zero

For finalizing the 2d simulations, the forces and moments at the z-direction must be set to zero. The fix enforce2d sets the forces and moments in the z-direction to zero. The syntax for the freeze fix is
``` bash
fix			nameofFix	group enforce2d
```
where nameofFix is user defined, the group of particles that the fix will operate, in our case group must be set to all.

#### Part II: MPLBM input files for running coupled simulations

For two dimensional DEM-LBM coupled simulations, input parameters are parsed from the file "input_params.txt". This file must include the following
- **Name of the conducted simulation** A hdf5 file of similar name must be also created
- **Number of blocks**: The number of blocks
- **Number of nodes in the x direction**
- **Number of nodes in the y direction**
- **x coordinate of the last node in the x-direction**
- **Fluid/Structure Interaction flag** 1: EQN scheme is on 2: Classical LBM code
- **Fluid viscocity**
- **Inlet velocity**
- **Number of Iterations for simulation relating to moving particles**: It does not used in DEM/LBM coupling through MUI. The DEM code sets the number of cycles
- **Number of increments** that will pass before exporting the flow field
- **Number of Gauss integration** that are used in the calculation of the solid fraction in the EQN scheme
- **Number of iterations** of the fluid code in the initialization phase of DEM/LBM simulations. Particles are stationary in this stage
- **Body force in the x-direction**
- **Body force in the y-direction**

An "input_params.txt" will look like this
``` bash
LD_couette
1
241
121
20.0
1
0.001
0.01
100000
10000
5
1
0.0
-0.315
```
Before mergining with HiLemms, to run DEM-LBM coupled simulations, it is required to create an h5 file. An h5 file can be created by running the grid generator geom2d_dev_seq.

#### Part III: Generating a h5 file

To generate an h5 file for two dimensional coupled simulations, the grid generator geom2d_dev_seq was developed. The grid generatator acceepts as inputs the name of the block, the grid size, the number of segments and the position of the upper right corner of the rectangle (The lower left corner point is by default set at (0,0)), the boundary conditions and the initial conditions.

#### Part IV: Implementation of boundary conditions

The use of the right boundary condition at the current 2d version requires a modification to the code and the generation of a new h5 file. The boundary conditions in the h5 file correspond to a given code. The most commonly boundary conditions that used in DEM/LBM simulations are
- **Equilibrium Diffuse Reflection BC** : Enforces a user defined velocity at the fluid boundary. **Code** 1014, **TypeNum**  Vertex_EQMDiffuseRefl
- **1st oder Density Extrapolation** : Used in the outlet to set the density to 1. **Code** 1006, **TypeNum** Vertex_ExtrapolPressure1ST
- **Far field BC**: Mimicks far field behavior **Code** 1016, **TypeNum** Vertex_FarField

The implementation of a given boundary condition, requires from the user  to modify the function ImplementBoundary (in evolution.cpp) and modify the BC in the required wall. The user can alter prescribed velocity and the type of BC. The  function looks like this
``` c++
    int* inletRng = BlockIterRng(0, g_BlockIterRngImin);
    //Input params: {Density, velocity in the x-direction, velocity in y direction, Temperature}
    Real givenInletVars[]{1, ux, 0 ,1}; // Input Parameters of inlet "wall"
    TreatDomainBoundary(givenInletVars,inletRng,Vertex_EQMDiffuseRefl);
    int* outletRng = BlockIterRng(0, g_BlockIterRngImax);
    Real givenOutletVars[]{1, ux, 0, 1};// Input Parameters of outlet "Wall"
    TreatDomainBoundary(givenOutletVars,outletRng,Vertex_ExtrapolPressure1ST);
    int* topRng = BlockIterRng(0, g_BlockIterRngJmax);
    //Real givenTopWallBoundaryVars[]{1, 0, 0};
    Real givenTopWallBoundaryVars[]{1,0.0, 0.0, 1};  // Input Parameters for top wall
    TreatDomainBoundary(givenTopWallBoundaryVars, topRng,Vertex_EQMDiffuseRefl);
    int* bottomRng = BlockIterRng(0, g_BlockIterRngJmin);
    Real givenBotWallBoundaryVars[]{1, 0.0, 0.0, 1};  // Input Parameters
    TreatDomainBoundary(givenBotWallBoundaryVars,bottomRng,Vertex_EQMDiffuseRefl);
```
