# Manual for the MPLB code

## Introduction
The multi-platform lattice Boltzmann code (MPLB) is a lattice Boltzmann solver written by using the oxford parallel library for structured mesh solvers (OPS).  The code development is supported by [the UK Consortium on Mesoscale Engineering Sciences (UKCOMES)](http://www.ukcomes.org/). By utilizing the functionalities of the OPS library, the code is capable of running on heterogeneous computing platform, and supporting finite-difference  lattice Boltzmann models. We are continually developing new functionalities into the code, including high-order lattice Boltzmann models (finite-difference scheme in the physical space), coupling scheme for solid-fluid two-phase flows etc.  

## Installation

In general, the environment for developing can be setup on any of Windows, Linux or Mac OS system, if we can have a MPI library and the corresponding parallel HDF5 library. For the Windows, we suggest the combination of Windows 10 and its Linux Subsystem, which provides almost the same environment to Linux and Mac OS.

### OPS library

The code relies on the OPS library, which provides the mesh management for parallel computing and the capability of running on heterogeneous computing platform. To facilitate a convenient post-processing, the OPS library also support parallel HDF5 IO capability, which needs the parallel HDF5 library. Other dependencies include OpenMPI/MPICH and one of tools, e.g., OPENCL, CUDA if we would like to use the corresponding hardware. For the detail of the OPS library, we refer to [here](http://www.oerc.ox.ac.uk/projects/ops) where the source code and manual are provided.  

####Installation of HDF5 library

* Windows 10 + WSL 
  We assume an Ubuntu 14.04 system by default, then type
  sudo apt install lbhdf5-openmpi-dev.
  Note: Current no GPU support for WSL
* Mac OS:
  brew install hdf5 --with-mpi
  *Opensuse:
  sudo zypper in \*\*\*\*
* Configuring the environment
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
The main solver can be compiled in two modes, i.e., the developing mode and the optimized mode. The developing mode is recommended to be used when debugging the code, while the optimized mode for production run. However, the developing mode is also fine for production running if we only would like to use CPUS and do not want to be bothered by the issues for compiling the optimized version.
####Developing mode
make lbm2d\_dev\_seq (sequential)
make lbm2d\_dev\_mpi (parallel)
#### Optimized mode
To compile the code in an optimized way, we have to fix two minor issues of the OPS library discussed as follows.
1. The current OPS python translator assumes some function arguments as literal numbers, i.e., not a variable. For instance, the "dim" parameter of the ops\_par\_loop call. To fix this, the python script FixConstantDefinition.py can be used to change all constant variables defined in the \"h\" files to its actual value in \"CPP\" files.
2. The OPS library tends to put source codes for "kernel function" into a ".h" file, which may be confusing to some extent. Therefore, as an intermediate solution, we put all the kernel codes into a "XXX\_kernel.h" file. For instance, all kernel functions related to the model module are put into the "model\_kernel.h" file. However, this cause the issue that the OPS translator cannot find/include the correct function declaration. Therefore, we provide the python script "FixKernelDeclaration.py".
  For convenience, we provide a bash script to automatize the process. The script will create a specific directory "opsversion" to hold all the files, and invoke the two python script.
#### Pre-processor
There is a simple pre-processor to help the preparation of the simulation.   To compile it, we run

```bash
make setupdomain
```

#### Post-processor
There is a simple post-processor written in Python, which can display contour plot and vector plot in both 2D (using matplotlib) and 3D (mayavi) for checking results. The post-processor can also convert the output to the format more friendly to other available software, e.g., plain HDF5 format (readable by Matlab) and  Tecplot HDF5 format. 

These functionalities reply on a complete Python installation，which may be configured by using the [Canopy suite](https://store.enthought.com/downloads/) or the [Anaconda distribution](https://www.anaconda.com/download/). In general either Python 3 or Python 2 will work. 

# A lid-driven cavity flow using stream-collision scheme
Currently the prepration of a simulation includes two major steps for both the solver and the preprocessor. 
* Setting the solover
1. Setting the model for the solver
  Since the MPLB code allow user to customize the component, the macroscopic variables, etc, users will need to define these in codes. Also, it will be the user responsibilities to keep the consistency for the whole simulation, which is of importance  in a number of routines (e.g., calculating the equilibrium function). In general, we keep the order of  0 - rho, 1- u  2 - v, and 3 - T for 2D problems, and  0 - rho, 1- u  2 - v, 3 - w, and 4 - T  for the BGK-type collision term. 
  Currently, we need to do the customization by hard-coding.  In model.cpp, we will first need to modify the routine SetupModel(), as shown below for the lid-driven cavity flow

  ```c++
  void SetupModel() {
      NUMCOMPONENTS = 1; //setting the total number of component 
      FEQORDER = 2; // using the second-order equilibrium function
      THERMALPROBLEM = 0; // 
      COMPOINDEX = new int[2 * NUMCOMPONENTS];
      COMPOINDEX[0] = 0;
      COMPOINDEX[1] = 14;
      SetupD3Q15Latt(); // setting lattice 
      //setting macroscopic variables, e.g., name, type and number  
      SetupMacroVars(); 
      DefineModelConstants();
  }
  ```
  In the routine SetupMacroVars(), we customize the macroscopic variables as below. 
  ```c++
  void SetupMacroVars() {
      // rho,u,v,w,T macroscopic variables must be stored in a specific order
      NUMMACROVAR = 4;
      VARIABLETYPE = new int[NUMMACROVAR];
      VARIABLETYPE[0] = (int)Variable_Rho;
      VARIABLETYPE[1] = (int)Variable_U;
      VARIABLETYPE[2] = (int)Variable_V;
      VARIABLETYPE[3] = (int)Variable_W;   
      VARIABLECOMPINDEX = new int[NUMMACROVAR];
      VARIABLECOMPINDEX[0] = 0;
      VARIABLECOMPINDEX[1] = 0;
      VARIABLECOMPINDEX[2] = 0;
      VARIABLECOMPINDEX[3] = 0;    
      MACROVARNAME.reserve(NUMMACROVAR);
      MACROVARNAME.push_back("rho");
      MACROVARNAME.push_back("u");
      MACROVARNAME.push_back("v");
      MACROVARNAME.push_back("w");
  }
  ```

2. Setting the scheme
  We mainly need to set the halo number needed by the chosen scheme. For a single-block application, the  current implementation of both stream-collision scheme and the finite difference scheme don't require halos. In multi-block applications, we will need halos and will need to manage these halos.

3. Setting the boundary for the solver
  Again, we need to set the halo number needed by the chosen boundary conditions. As we choose to implement the boundary condition on grid, most of boundary conditions will not rely on halos. The exception is the periodic boundary condition which will require halos. 
  We also need to manually set the boundary in the following routines in the evolutions3d.cpp

  ```c++
  void ImplementBoundary3D() {
      //TreatEmbeddedBoundary3D();   
      int* inletRng = BlockIterRng(0, g_BlockIterRngImin);    
      Real givenInletVars[]{1, 0, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenInletVars, inletRng, Vertex_EQMDiffuseRefl);
  
      int* outletRng = BlockIterRng(0, g_BlockIterRngImax);
      Real givenOutletVars[] = {1, 0, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenOutletVars, outletRng,
                           Vertex_EQMDiffuseRefl);
      
      int* topRng = BlockIterRng(0, g_BlockIterRngJmax);
      // Real givenTopWallBoundaryVars[]{1, 0, 0};
      Real givenTopWallBoundaryVars[]{1, 0.01, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenTopWallBoundaryVars, topRng,
                           Vertex_EQMDiffuseRefl);
  
      int* bottomRng = BlockIterRng(0, g_BlockIterRngJmin);
      Real givenBotWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenBotWallBoundaryVars, bottomRng,
                           Vertex_EQMDiffuseRefl);
      
      int* backRng = BlockIterRng(0, g_BlockIterRngKmin);
      Real givenBackWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenBackWallBoundaryVars, backRng,
                           Vertex_EQMDiffuseRefl);
  
      int* frontRng = BlockIterRng(0, g_BlockIterRngKmax);
      Real givenFrontWallBoundaryVars[]{1, 0, 0, 0};  // Input Parameters
      TreatBlockBoundary3D(givenFrontWallBoundaryVars, frontRng,
                           Vertex_EQMDiffuseRefl);
  }
  ```

4. Setting the flow field for the solver

  Currently, we need to set a few parameters including the domain size, spatial dimension etc. 

  ```c++
  void SetupFlowfieldfromHdf5() {
      CASENAME = "Cavity3D";  // Input parameter
      SPACEDIM = 3;
      g_BlockNum = 1;  // Input parameter
      g_BlockSize = new int[g_BlockNum * SPACEDIM];
      g_BlockSize[0] = 33;  // Input parameters
      g_BlockSize[1] = 33;  // Input parameters
      g_BlockSize[2] = 33;  // Input parameters
      KN = new Real[ComponentNum()];
      KN[0] = 0.001;  // Input parameters
      // All above parameters should be written down by the python script
      Real minDx{1./32};  // Input parameters at this moment
      Real minDy{1./32};  // Input parameters at this moment   
      g_dt = minDx / SoundSpeed();  // stream-collision
      g_HaloDepth = HaloPtNum();
      DefineVariablesFromHDF5();
      // DefineHaloTransfer();
      // above calls must be before the ops_partition call
      ops_partition((char*)"LBM");
  }
  ```

  In the simulation, we tend to use non-dimensional system so that the relaxation time in the non-dimensional system is the Knudsen number.Nevertheless, one can also use the common lattice system. 5. Setting the running parameters and scheme

  In lbm3d.cpp, the main() function, we set up the maximum iteration number and the checkPeriod. To use the standard stream-collision scheme, we use the StreamCollision3D() function defined in evolution3d.cpp. 

  ```c++
  const int maxIter =100000;
  const int checkPeriod = 10000;
  for (int iter = 0; iter < maxIter; iter++) {
  	StreamCollision3D();//Stream-Collision scheme
      //TimeMarching();//Finite difference scheme + cutting cell
      if ((iter % checkPeriod) == 0) {           
           UpdateMacroVars3D();
           CalcResidualError3D();
           DispResidualError3D(iter,checkPeriod*g_dt);
           WriteFlowfieldToHdf5(iter);
           WriteDistributionsToHdf5(iter);
           WriteNodePropertyToHdf5(iter);   
       }
   }
  ```

  **After finishing all these steps, we shall compile the solver.**
* Setting the preprocessor
1. Setting  the initial condition.
  We can define the initial condition in the  go to seup\_comput\_domain_kernel.h and modify the kerSetInitialMacroVars() routine. In this routine, the variable *coordinates* store the spatial coordinates by the order 0 \-x, 1 \- y,2 \- z. The variable *idx* store the index of present node by the order 0 \-x, 1 \- y,2 \- z.. The variable *macroVars* corresponds to the macroscopic variables (a pointer pointing to an array). 
  if we define the initial condition by inserting new codes in seup\_comput\_domain_kernel.h, we need to recompile the preprocessor. 

  ```c++
  void KerSetInitialMacroVars(const Real* coordinates, const int* idx,
                              Real* macroVars) {
  #ifdef OPS_3D
      macroVars[OPS_ACC_MD2(0, 0, 0, 0)] = 1;         // rho
      macroVars[OPS_ACC_MD2(1, 0, 0, 0)] = 0;         // u
      macroVars[OPS_ACC_MD2(2, 0, 0, 0)] = 0;         // v
      macroVars[OPS_ACC_MD2(3, 0, 0, 0)] = 0;         // w
  #endif
  }
  ```

  **After we define the macroscopic initial  condition, we need to recompile the the preprocessor.** 

2. Setting the domain and boundary condition
  We run  "setup_comput_domain" to set the domain and boundary condition. It will ask a number of questions as shown below. 
  Please input CASENAME:Cavity3D
  Please input the space dimensional:3 // 2 for 2D problem
  Please input total number of blocks:1 // currently only 1 block is supported
  Please input the dimension of each block (grid point number):
  Block 0:
  At the coordinate 0: 33 // 0 - x 32 cells with 33 grid points
  At the coordinate 1:33 // 1 - y
  At the coordinate 2:33 // 2 - z
  Now constructing each block/domain...
    Constructing coordinates for Block 0
     Input total number of segments at Coordinate 0:
  ​     At the coordinate 0:1
  ​     Input the total number of cells and the position of the final grid point at each segment:
  ​    Total cell number:32
  ​    Position of the final grid:1 ​    Input total number of segments at Coordinate 1:
  ​    At the coordinate 1:1
  ​    Input the total number of cells and the position of the final grid point at each segment:
  ​    Total cell number:32
  ​    Position of the final grid:1
  ​    Input total number of segments at Coordinate 2:
  ​    At the coordinate 2:1
  ​    Input the total number of cells and the position of the final grid point at each segment:
  ​    Total cell number:32
  ​    Position of the final grid:1

  Please input boundary at left:1014 
  Please input boundary at right:1014
  Please input boundary at bottom:1014
  Please input boundary at top:1014
  Please input boundary at front:1014
  Please input boundary at back:1014
  Please input boundary at leftBottom:1014
  Please input boundary at LeftTop:1014
  Please input boundary at rightBottom:1014
  Please input boundary at rightTop:1014
  Please input boundary at leftBack:1014
  Please input boundary at leftFront:1014
  Please input boundary at rightBack:1014
  Please input boundary at rightFront:1014
  Please input boundary at bottomBack:1014
  Please input boundary at bottomFront:1014
  Please input boundary at topBack:1014
  Please input boundary at topFront:1014
  Please input boundary at leftBottomBack:1014
  Please input boundary at LeftBottomFront:1014
  Please input boundary at LeftTopBack:1014
  Please input boundary at LeftTopFront:1014
  Please input boundary at RightBottomBack:1014
  Please input boundary at RightBottomFront:1014
  Please input boundary at RightTopBack:1014
  Please input boundary at RightTopFront:1014

  The segments are intended to support non-uniform mesh in conjunction with the finite-difference LBM (see the illustration below).  That is, we can have different cell size in different segment. For the standard stream-collision, only uniform cartesian mesh are supported. 
* Post-processing
  For post-processing, we recommend to convert the output data to a propert format that can be recognized by a specialized code/software, e.g., Tecplot.

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
