/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
*/

/** @brief The main source file for direct DEM-LBM sims with the PSM
 *  @author Chrysovalantis Tsigginos
 */

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "boundary.h"
#include "evolution3d.h"
#include "flowfield.h"
#include "model.h"
#include "scheme.h"
#include "type.h"
#include "configuration.h"
#include "psm.h"
#include "dem_handle3d.h"
#include "iterate_dem.h"
#include "periodic.h"
#include "read_input.h"




void SetInitialMacrosVars() {
	Real initialVars[SPACEDIM + 1];

	initialVars[0] = 1.0;
	initialVars[1] = 0.0;
	initialVars[2] = 0.0;
	initialVars[3] = 0.0;
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
        ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                     g_Block[blockIdx], SPACEDIM, iterRng,
                     ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW),
                     ops_arg_dat(g_CoordinateXYZ[blockIdx], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_gbl(initialVars, SPACEDIM+1, "double", OPS_READ));
    }
}

void UpdateMacroscopicBodyForce(const Real time) {}


void simulate(configData input, SizeType restartStep = 0) {

	//ops_printf("I entered in the simulate\n");
	DefineCase(input.casename, input.spacedim); //set transient to true

	//Define geometry
	DefineBlocks(input.blockNum, input.gridBlock, input.dx, input.gridStart); //Define coordinates, nodeType and Geometry property.

	//Define lattice model
	std::vector<std::string> compoNames{"Fluid"};
	std::vector<SizeType> compoid{0};
	std::vector<std::string> lattNames{"d3q19"};
	DefineComponents(compoNames, compoid, lattNames, input.restartStep); //Define g_f and g_fStage

	ops_printf("Distribution functions defined at timestep = %d\n", restartStep);


	//Define macroscopic variables
	std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
	                                             Variable_V, Variable_W};
	std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
	std::vector<SizeType> macroVarId{0, 1, 2, 3};
	std::vector<SizeType> macroCompoId{0, 0, 0, 0};
	DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId, restartStep); //Define macroscopic variables and error handling

	//Define PSM variables
	DefinePSMVariables(restartStep, input.gamma);

	//Define bodyForces
	std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
	std::vector<SizeType> bodyForceCompoId{0};
	DefineBodyForce(bodyForceTypes, bodyForceCompoId);

	//Define collision model
	std::vector<CollisionType> collisionTypes{Collision_BGKIsothermal2nd};
	std::vector<SizeType> collisionCompoId{0};
	DefineCollision(collisionTypes, collisionCompoId);

	//Define bodyForces
	SetupBodyForces();

	ops_printf("MP-LBM: BodyForces have been set\n");
	//Define collision scheme
	SchemeType scheme{Scheme_StreamCollision};
	DefineScheme(scheme); //Defines halodepth

	// Setting boundary conditions
	SizeType blockIndex{0};
	SizeType componentId{0};

//Boundary conditions
	BoundaryScheme boundaryType[6] = {BoundaryScheme::EQMDiffuseRefPois, BoundaryScheme::FreeFlux, BoundaryScheme::EQMDiffuseRefl, BoundaryScheme::EQMDiffuseRefl,
									  BoundaryScheme::Periodic,BoundaryScheme::Periodic,};



	DefinePeriodicBoundariesRestart(input.boundaries.BoundaryType, input.boundaries.DirBc);


	//Define VertexType for different boundaries

	VertexType vertexType[6];

	if (periodicFlag == 0) {
		for (int iDim = 0; iDim < 6; iDim++) {
			vertexType[iDim] = VertexType::Wall;
		}
	}
	else if (periodicFlag == 1) {
		for (int iDim = 0; iDim < 3; iDim++) {
			int direction = input.boundaries.DirBc[iDim];
			if (periodic[direction]==1) {
				vertexType[2*iDim] = VertexType::Periodic;
				vertexType[2*iDim + 1] = VertexType::Periodic;
			}
			else {
				vertexType[2*iDim] = VertexType::Wall;
				vertexType[2*iDim + 1] = VertexType::Wall;
			}
		}
	}

	//Define walls

	//Input Velocities-Please note that periodic boundaries should come first

	std::vector<VariableTypes> macroVarTypesatBoundary{Variable_U, Variable_V,
	                                                       Variable_W};

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[0],
						input.boundaries.BoundaryType[0], macroVarTypesatBoundary,
						input.boundaries.velocities[0], vertexType[0]);

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[1],
						input.boundaries.BoundaryType[1], macroVarTypesatBoundary,
						input.boundaries.velocities[1], vertexType[1]);

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[2],
						input.boundaries.BoundaryType[2], macroVarTypesatBoundary,
						input.boundaries.velocities[2], vertexType[2]);

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[3],
						input.boundaries.BoundaryType[3], macroVarTypesatBoundary,
						input.boundaries.velocities[3], vertexType[3]);

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[4],
						input.boundaries.BoundaryType[4], macroVarTypesatBoundary,
						input.boundaries.velocities[4], vertexType[4]);

	DefineBlockBoundary(blockIndex, componentId, input.boundaries.BoundarySurfaces[5],
						input.boundaries.BoundaryType[5], macroVarTypesatBoundary,
						input.boundaries.velocities[5], vertexType[5]);


	//setup initial condition

	std::vector<InitialType> initType{Initial_BGKFeq2nd};
	std::vector<SizeType> initalCompoId{0};
	DefineInitialCondition(initType,initalCompoId);

	ops_printf("Ready to partition the domain\n");
	Partition();



	if (restartStep == 0) {
		SetInitialMacrosVars();
		PreDefinedInitialCondition3D();
	}
	else
		 RestartMacroVars4SteadySim();


	std::vector<Real> tauRef{input.tauRef};

	SetTauRef(tauRef);

	SetTimeStep(input.dx / SoundSpeed());

	SetDxLBM(input.dx);
	SetupParticles();



	if (periodicFlag == 1) {
		InitializePeriodic();
		ops_printf("MP-LBM: Initialize periodic variables\n");
	}


	//SetupParticles work correct !!!
	if (input.muiFlag==0)
		muiInterfaceFlag = 0;
	else
		muiInterfaceFlag = 1;


	if (muiInterfaceFlag == 1) {
		DefineMUIInterface(input.Rmax);
		ops_printf("Setting up communication\n");
		SetupMUICommunication(10000000);  //To be corrected in the next stage
	}

	ops_printf("I setup the communication\n");
	//To add here additional fumctions for running the code

	int savingFlag;
	if (input.storageOption.compare(0,input.storageOption.size(), "noSaving")==0) {
		ops_printf("No data are saved during the simulation\n");
		savingFlag = 0;
	}
	else {
		ops_printf("Data are saved during the simulation\n");
		savingFlag = 1;
	}


	SizeType maxStep;

	if (restartStep == 0) {

		const Real convergenceRate{10000};
		const int checkPeriod{1};
		int maxIters = 0;

		ops_printf("Before entering to InitializeDEM-LBM: ConvergenceRate: %e maxIters: %d, checkPeriod: %d\n",
				convergenceRate, maxIters, checkPeriod);

		InitializeDEMLBM(convergenceRate, maxIters, checkPeriod, maxStep, savingFlag);

		if (muiInterfaceFlag==0)
			maxStep = restartStep + input.iters;

	}
	else {
		//To be added
		InitializeDEMLBMRestart(maxStep, restartStep, input.iters, savingFlag);
	}

	ops_printf("Restart timeStep = %d\n", restartStep);

	IterateDEMLBM(maxStep, input.exportSteps, restartStep, savingFlag);

	if (muiInterfaceFlag == 1)
		DestroyMUIInterface();


}

void simulate() {

	//Define simulation space and simulation name
	std::string caseName{"Classic_case"};
	SizeType spaceDim{3};
	DefineCase(caseName, spaceDim); //set transient to true

	//Define geometry
	SizeType blockNum{1};
	//std::vector<SizeType> blockSize{75, 75, 75};
	//Real meshSize{20. / 74};
	std::vector<SizeType> blockSize{641, 41, 41};
	Real Lx = 16.0;

	Real meshSize{Lx / (static_cast<Real>(blockSize[0])-1.0)};
	std::vector<Real> startPos{0.0, 0.0, 0.0};
	DefineBlocks(blockNum, blockSize, meshSize, startPos); //Define coordinates, nodeType and Geometry property.

	//Define lattice model
	std::vector<std::string> compoNames{"Fluid"};
	std::vector<SizeType> compoid{0};
	std::vector<std::string> lattNames{"d3q19"};
	DefineComponents(compoNames, compoid, lattNames); //Define g_f and g_fStage

	//Define macroscopic variables
	std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
	                                             Variable_V, Variable_W};
	std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
	std::vector<SizeType> macroVarId{0, 1, 2, 3};
	std::vector<SizeType> macroCompoId{0, 0, 0, 0};
	DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId); //Define macroscopic variables and error handling

	//Define PSM variables
	DefinePSMVariables();

	//Define bodyForces
	std::vector<BodyForceType> bodyForceTypes{BodyForce_None};
	std::vector<SizeType> bodyForceCompoId{0};
	DefineBodyForce(bodyForceTypes, bodyForceCompoId);

	//Define collision model
	std::vector<CollisionType> collisionTypes{Collision_BGKIsothermal2nd};
	std::vector<SizeType> collisionCompoId{0};
	DefineCollision(collisionTypes, collisionCompoId);

	//Define bodyForces
	SetupBodyForces();

	ops_printf("MP-LBM: BodyForces have been set\n");
	//Define collision scheme
	SchemeType scheme{Scheme_StreamCollision};
	DefineScheme(scheme); //Defines halodepth

	// Setting boundary conditions
	SizeType blockIndex{0};
	SizeType componentId{0};





// Original versopm
//	BoundaryScheme boundaryType[6] = {BoundaryScheme::EQMDiffuseRefl,BoundaryScheme::EQMDiffuseRefl, BoundaryScheme::EQMDiffuseRefl,
//  BoundaryScheme::EQMDiffuseRefl, BoundaryScheme::Periodic,BoundaryScheme::Periodic};


	BoundaryScheme boundaryType[6] = {BoundaryScheme::EQMDiffuseRefPois, BoundaryScheme::FreeFlux, BoundaryScheme::EQMDiffuseRefl, BoundaryScheme::EQMDiffuseRefl,
									  BoundaryScheme::Periodic,BoundaryScheme::Periodic,};

	DefinePeriodicBoundaries(boundaryType);



	//Define VertexType for different boundaries

	VertexType vertexType[6];

	SetupVertexType(vertexType);
	if (periodicFlag == 0) {
		for (int iDim = 0; iDim < 6; iDim++) {
			vertexType[iDim] = VertexType::Wall;
		}
	}
	else if (periodicFlag == 1) {
		for (int iDim = 0; iDim < 3; iDim++) {
			if (periodic[iDim]==1) {
				vertexType[2*iDim] = VertexType::Periodic;
				vertexType[2*iDim + 1] = VertexType::Periodic;
			}
			else {
				vertexType[2*iDim] = VertexType::Wall;
				vertexType[2*iDim + 1] = VertexType::Wall;
			}
		}
	}

	//Define walls

	//Input Velocities-Please note that periodic boundaries should come first

	std::vector<VariableTypes> macroVarTypesatBoundary{Variable_U, Variable_V,
	                                                       Variable_W};
	std::vector<Real> noSlipStationaryWall{0, 0, 0};
	std::vector<Real> noSlipMovingWall{0.005, 0, 0};
	std::vector<Real> noSlipMovingWall1{-0.0, 0, 0};


	// Back Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Back,
							boundaryType[4], macroVarTypesatBoundary,
							noSlipMovingWall, vertexType[4]);

	//Front Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Front,
							boundaryType[5], macroVarTypesatBoundary,
							noSlipMovingWall1, vertexType[5]);





	// Left Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Left,
	                   	boundaryType[0], macroVarTypesatBoundary,
						noSlipMovingWall, vertexType[0]);
	// Right Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Right,
	                    boundaryType[1], macroVarTypesatBoundary,
						noSlipMovingWall1, vertexType[1]);



	// Bottom Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Bottom,
	                    boundaryType[2], macroVarTypesatBoundary,
						noSlipStationaryWall, vertexType[2]);

	//Top Wall
	DefineBlockBoundary(blockIndex, componentId, BoundarySurface_Top,
						boundaryType[3], macroVarTypesatBoundary,
						noSlipStationaryWall, vertexType[3]);














	//setup initial condition

	std::vector<InitialType> initType{Initial_BGKFeq2nd};
	std::vector<SizeType> initalCompoId{0};
	DefineInitialCondition(initType,initalCompoId);

	ops_printf("Ready to partition the domain\n");
	Partition();
	SetInitialMacrosVars();
	PreDefinedInitialCondition3D();
	std::vector<Real> tauRef{0.001};
	SetTauRef(tauRef);
	SetTimeStep(meshSize / SoundSpeed());

	SetDxLBM(meshSize);
	SetupParticles();



	if (periodicFlag == 1) {
		InitializePeriodic();
		ops_printf("MP-LBM: Initialize periodic variables\n");
	}


	//SetupParticles work correct !!!
	muiInterfaceFlag = 1;
	double Rmax = 0.3;
	if (muiInterfaceFlag == 1) {
		DefineMUIInterface(Rmax);
		ops_printf("Setting up communication\n");
		SetupMUICommunication(10000000);
	}

	ops_printf("I setup the communication\n");
	//To add here additional fumctions for running the code

	const Real convergenceRate{10000};
	const int checkPeriod{1};
	SizeType maxStep;
	int maxIters = 0;

	ops_printf("Before entering to InitializeDEM-LBM: ConvergenceRate: %e maxIters: %d, checkPeriod: %d\n",
				convergenceRate, maxIters, checkPeriod);

	InitializeDEMLBM(convergenceRate, maxIters, checkPeriod, maxStep);

	IterateDEMLBM(maxStep, 10000);

	if (muiInterfaceFlag == 1)
		DestroyMUIInterface();



}


int main(const int argc,const  char** argv) {


	//OPS initialization
	ops_init(argc, argv, 1);


	double ct0, ct1, et0, et1;

	ops_timers(&ct0, &et0);


	if (argc <= 1) { //hard-code simulation
		simulate();
	}
	else {
		//This part will later be modified to JSON file input. At the moment will be a txt reading
//		if (argc > 3) {
//			ops_printf("ERROR: Number of arguments less than two\n");
//			exit(EXIT_FAILURE);
//		}



		//Compare strings

		if (strcmp(argv[1],"restart")==0) {
			configData readData;
			ops_printf("I entered here\n");
			if (argc < 3) {
				ops_printf("ERROR: Restart Timestep not defined\n");
				exit(EXIT_FAILURE);
			}

			long int timeStep = atoi(argv[2]);
			ReadInputFile(readData, timeStep); //Later to pass to Json
			simulate(readData, readData.restartStep);

		}
		else if (strcmp(argv[1],"readFile")==0) {
			configData readData;
			ReadInputFile(readData); //Later to pass to Json
			simulate(readData);

		}
		else {
			ops_printf("ERROR: This option is not supported. Supported options: restart and readFile\n");
		}




	}


/*	if (argc>1 && argc <=2) { //start a simulation from input file
		 std::string configFileName(argv[1]);
		 ReadConfiguration(configFileName);
		 simulate(Config());
	}

	if (argc> 2 && argc <=3) {
		 std::string configFileName(argv[1]);
		 ReadConfiguration(configFileName); //NEED TO DISCUSS THIS
		 const SizeType timeStep{std::stoi(argv[2])};
		 simulate(Config(),timeStep);
		 restartFlag = 1;
	}

*/
	 ops_timers(&ct1, &et1);
	 ops_printf("\nTotal Wall time %lf\n", et1 - et0);
	    //Print OPS performance details to output stream
	 ops_timing_output(stdout);
	 ops_exit();



}

