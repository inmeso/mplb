
/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]s
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

/*! @brief   Define functions and variables for handling particles
 * @author  C. Tsigginos
 *
 * @details: Functions for handling particle data, function for iterating functions
 * 			 in DEM-LBM coupled simulations
 */

#include "dem_handle.h"
#include "flowfield_host_device.h"
#include "boundary.h"
#include "mui_wrapper.h"
#include <string>
#include <limits>
#include "dem_data.h"
#include <vector>

void SetDemLbMParams(InteractionData* data, bool flag, bool muiOn, Real convergeRate,
		SizeType checkperiod, SizeType timeStep, SizeType checkPeriodSteady,
		SizeType maximumIterations, std::string particleType) {

	data->restartFlag = flag;
	if (data->restartFlag)
		data->nStart = timeStep;
	else
		data->nStart = 0;

	data->muiFlag = muiOn;
	data->convergenceRate = convergeRate;
	data->checkPeriod = checkperiod;
	data->checkPeriodStS = checkPeriodSteady;
	data->maxIters = maximumIterations;
	//Finding Particle type
	std::string s1{"spherical"};
	std::string s2{"quadratic"};
	std::string s3{"mesh"};
	if (particleType.compare(s1) == 0)
		data->particleShape = spherical;
	else if (particleType.compare(s2) == 0)
		data->particleShape = quadratic;
	else if (particleType.compare(s3) == 0)
		data->particleShape = mesh;
	else {
		ops_printf("MPLB: This type of particle is not supported\n");
		ops_printf("MPLB: Reversing to spherical particles\n");
		data->particleShape = spherical;
	}
}



/*void DefineInteractionModel(std::vector<FSIType> FluidParticleInteractionType,
		int* porosModel, std::vector<int> fsiCompoId, Real* forceUser,
		double gamma, SizeType timeStep) {

	int noComp = FluidParticleInteractionType.size();
	int noFSICompo = fsiCompoId.size();
	FsiBase* model;
	int idCompo;

	std::map<int, Component> components = g_Components();
	int numComponents = ComponentNum();

	if (noComp != numComponents) {
		ops_printf("Error: Number of FSI components inconsistent with actual components\n");
		assert(noComp != numComponents);
	}

	if (noFSICompo != numComponents) {
		ops_printf("Error: Number of fsiCompoId differs from actual number of components\n");
		assert(noFSICompo != numComponents);
	}

	for (int iComp = 0; iComp < noComp; iComp++) {

		switch (FluidParticleInteractionType[iComp]) {
			case Model_None:
				model = new FsiBase( c.at(iComp), SpaceDim(), forceUser, false,  porosModel[0], gamma);
				break;
			case Model_PSM:
				model = new Psm(c.at(iComp), SpaceDim(),true , porosModel[0], gamma, porosModel[1], porosModel[2],
								porosModel[3]);
				break;
			case Model_Prati:
				model = new Prati(c.at(iComp), SpaceDim(), true, porosModel[0], gamma, porosModel[1], porosModel[2],
						porosModel[3]);
				break;
			default:
				ops_printf("The chosen model is not supported\n");
				exit(EXIT_FAILURE);
		}
		idCompo = components.at(iComp).id;
		fluidPartInteractionModels.emplace(idCompo, model); //TODO DEfine it

	}

}
*/

void SetupBlockParticles(InteractionData data, SizeType maxStep) {

	//Defined if block is owned
	DefineBlockOwnership();

	//SetupGlobalBox:
	DefineLocalBoxBound();
	DefineGlobalBlockBox();

	if (data.muiFlag)
		SetMuiDomains( maxStep);

}
void UserDefineInputOutputParams(std::vector<std::string>& inputDemParams,
		std::vector<std::string>& outPutDemParams) {

}
void CreateMuiInterface(InteractionData data, std::vector<FluidParticleModel> models, Real skin) {


	std::vector<std::string> inputDemParams;
	std::vector<std::string> outputDemParams;
	std::vector<std::string> particleParams;

	if (data.particleShape== quadratic) {
		particleParams.push_back("A");
		particleParams.push_back("B");
		particleParams.push_back("C");
		particleParams.push_back("a1");
		particleParams.push_back("b1");
		particleParams.push_back("c1");

		//quaternions: Particle rotations
		particleParams.push_back("q1");
		particleParams.push_back("q2");
		particleParams.push_back("q3");
	}

	UserDefineInputOutputParams(inputDemParams, outputDemParams);

	DefineMuiInterface(&BlockParticleList,
						inputDemParams,  particleParams, outputDemParams,skin);



}

void SetupParticleBoxes(InteractionData data) {


	DefineBlockOwnership();



	DefineLocalBoxBound();

	DefineGlobalBlockBox();

	if (data.muiFlag)
		SetMuiDomains(100000000);


}

void DefineGlobalBlockBox() {
	int spaceDim1 = SpaceDim();
	Real xMinMaxTmp[spaceDim1];
	Real xb[2 * spaceDim1], xMin[spaceDim1], xMax[spaceDim1], xTemp[2 * spaceDim1];
	for (auto& idBlock : BlockParticleList) {
		BlockParticles& particles = idBlock.second;
		const int blockIdx{idBlock.first};
		bool owned = particles.OwnedStatus();

		if (owned) {
			particles.GetLocalBound(xb);

			for (int iDim = 0; iDim < spaceDim1; iDim++) {
				xMin[iDim] = xb[2 * iDim];
				xMax[iDim] = xb[2 * iDim + 1];
			}
		}
		else {
			for (int iDim = 0; iDim < spaceDim1; iDim++) {
				xMin[iDim] = std::numeric_limits<Real>::max();
				xMax[iDim] = -1.0 * std::numeric_limits<Real>::max();
			}
		}


#ifdef OPS_MPI

		if (ops_num_procs()>1) {

			MPI_Allreduce(xMin, xMinMaxTmp, spaceDim1, MPI_DOUBLE, MPI_MIN, OPS_MPI_GLOBAL);

			for (int iDim = 0; iDim < spaceDim1; iDim++) {
				xTemp[2 * iDim] = xMinMaxTmp[iDim];
			}

			MPI_Allreduce(xMax, xMinMaxTmp, spaceDim1, MPI_DOUBLE, MPI_MAX, OPS_MPI_GLOBAL);
			for (int iDim = 0; iDim < spaceDim1; iDim++) {
				xTemp[2 * iDim + 1] = xMinMaxTmp[iDim];
			}
		}
		else {
			for (int iDim = 0; iDim < 2 * spaceDim1; iDim++)
				xTemp[iDim] = xb[iDim];
		}

#else
		for (int iDim = 0; iDim < 2 * spaceDim1; iDim++)
			xTemp[iDim] = xb[iDim];
#endif

#ifdef CPU
#if DebugLevel >= 2
		ops_printf("----------------------------------------------------------------------\n");
		printf(" Block %d: ", particles.GetBlock().ID());
		for (int iDir = 0; iDir < spaceDim1; iDir++)
			printf(" [%f %f] ", xTemp[2*iDir], xTemp[2 * iDir + 1]);
		printf("\n");

#endif
#endif
		particles.SetGlobalBound(xTemp);
	}
}


void StreamCollisionFSI3D(int flag) {

#if DebugLevel>1
	ops_printf("Particle mapping\n");
#endif

	ParticleMapping(flag);

#if DebugLevel >= 1
	ops_printf("Calculating the macroscopic variables...\n");
#endif

	UpdateMacroVars3D();


	CopyBlockEnvelopDistribution3D(g_fStage(), g_f());

	PostVelocityFSIFunctions();

/*
 * TODO Shift to bodyForce term global calculation
 */
/*#if DebugLevel >= 1
	ops_printf("Calculating the mesoscopic body force term...\n");
#endif
	    UpdateMacroscopicBodyForce(time);
	    PreDefinedBodyForce3D();
	#if DebugLevel >= 1
	    ops_printf("Calculating the collision term...\n");
	#endif
*/

#if DebugLevel >=1
	ops_printf("PreCollision FSI functions\n");
#endif
	PreCollisionFSIFunctions();


#if DebugLevel >=1
	ops_printf("Collisions in FSI model\n");
#endif
	FluidParticleCollisions();


	//PreDefinedCollision3D();

#if DebugLevel >= 1
	ops_printf("Updating the halos...\n");
#endif

	TransferHalos();

#if DebugLevel >= 1
	 ops_printf("Streaming...\n");
#endif

	 Stream3D();

#if DebugLevel >=  1
	ops_printf("PostStreaming FSI functions\n");
#endif
	 PostStreamingFSIFunctions();


#if DebugLevel >= 1
	 ops_printf("Implementing the boundary conditions...\n");
#endif
	    // TODO This function shall be inside evolution3D.cpp
	    // TODO The data structure for BCs shall be inside boundary module
	 ImplementBoundary3D();

#if DebugLevel >=1
	 ops_printf("Calculate Drag Force\n");
#endif
	 CalculateParticleMomentum();


}

void IterateFSI(InteractionData data, int savingFlag) {

	Real dtD;
	int flag{1};
	if (data.Nf <= 1)
		dtD = TimeStep();
	else
		dtD = data.dtDEM;

	SizeType initStep = data.nStart;
	SizeType currentStep = initStep;
	for (SizeType step = 0; step < data.nSteps; step++) {
		currentStep++;

		if (step == 0) {
			if (data.restartFlag)
				flag = 0;
		}
		ops_printf("Current step %d\n", currentStep);
		if ((currentStep % data.Npl)==0) {
			if (data.muiFlag) {
				ExtractParticleData(currentStep);
			}
			InitializeDragForce();

			for (SizeType jfl = 0; jfl < data.Nf; jfl++) {
				StreamCollisionFSI3D(flag);
			}
		}
		flag = 0;
		CalculateDragForce(dtD, currentStep);

		if (data.muiFlag) {
			SendParticleData(currentStep);
		}

		if (savingFlag==1) {
			if ((currentStep % data.checkPeriod) == 0) {

				UpdateFPIVelocities3D();

                CalcResidualError3D();
                DispResidualError3D(currentStep, data.checkPeriod);

                WriteFlowfieldToHdf5(currentStep);
                WriteDistributionsToHdf5(currentStep);
                WriteNodePropertyToHdf5(currentStep);

                //TODO ADD FSI distributions
                WriteFPIDataToHdf5(currentStep);
			}
		}
	}

	if (savingFlag == 1) {
		UpdateFPIVelocities3D();
		WriteFlowfieldToHdf5(currentStep);
		WriteDistributionsToHdf5(currentStep);
		WriteNodePropertyToHdf5(currentStep);

		//TODO ADD FSI distributions
		WriteFPIDataToHdf5(currentStep);
	}



}

void IterateFSI(Real convergenceRate,const SizeType checkPointPeriod,const SizeType maxIters) {

	Real dtD = TimeStep();


	InitializeDragForce();

	ParticleMapping(1);

	Real residualError{1};
	SizeType iter = 0;
	do {
		InitializeDragForce();

		StreamCollisionFSI3D(0);
		iter +=1;
        if ((iter % checkPointPeriod) == 0) {
        	 UpdateFPIVelocities3D();
             CalcResidualError3D();
             residualError = GetMaximumResidual(checkPointPeriod);
             DispResidualError3D(iter, checkPointPeriod);
             printf("Iter %d: Error = %12.9e Target = %12.9e\n", iter,residualError, convergenceRate);
        }

        CalculateDragForce(TimeStep(), 0);
        if (iter > maxIters)
        	break;



	} while (residualError >= convergenceRate);

	UpdateFPIVelocities3D();
	WriteFlowfieldToHdf5(0);
	WriteDistributionsToHdf5(0);
	WriteNodePropertyToHdf5(0);
	WriteFPIDataToHdf5(0);

}

void SetupRestartSimulation(SizeType timestep, SizeType endStep) {

	UpdateRegions(timestep, endStep);

	UpdateOldParticleLocation();

	UpdateParticleMappingDragForceRestart(timestep);
	//TODO Update FSI Schemes to be ready for the simulations
}

void SetupDEMLBM(InteractionData& data) {

	int myRank, timeSetup;
	Real alpha;
	int flags[2];

#ifdef OPS_MPI
	myRank = ops_get_proc();
#endif


	if (data.muiFlag) {
		//TODO Pass mui data to upgrade data
		 ExtractParticleData(data.nStart);
		 ExtractSimulationData(data.nStart, data.nStart, data.nSteps, alpha, flags,
				 data.particleShape);
		 UpdateRegions(data.nStart, data.nSteps + data.nStart);
	}
	else {
		switch(data.particleShape) {
			case spherical:
				ReadParticleDataSpherical();
				break;
			default:
				ops_printf("Only spherical particles are currently supported\n");
		}
		alpha = 1.0;
		data.nSteps = 100;

	}


	if (alpha < 1.0) {
		data.Npl = static_cast<SizeType>(1.0 / alpha);
		data.Nf = 1;
	}
	else {
		data.Npl = 1;
		data.Nf = static_cast<SizeType>(alpha);
	}

	data.dtDEM = alpha * TimeStep();

#ifdef CPU
#if DebugLevel >= 2
	ops_printf("Rank %d: Nf = %d Np = %d, Starting point: %d Number of steps: %d dtDEM = %f\n",
			data.Nf, data.Npl, data.nStart, data.nSteps, data.dtDEM);
#endif
#endif



	ops_printf("Initializing LBM-DEM coupled simulation\n");
	if (!data.restartFlag)
		IterateFSI(data.convergenceRate, data.checkPeriodStS, data.maxIters);
	else
		SetupRestartSimulation(data.nStart, data.nSteps + data.nStart);

	if (data.muiFlag)
		SendParticleData(data.nStart);




}

//Read particle Data for coupled simulations
void ReadParticleDataSpherical() {

	int myRank;
	int Nparticles, Nsize;

	myRank = 0;
#ifdef OPS_MPI
	myRank = ops_get_proc();
#endif

	FILE *cfilex;
	cfilex = fopen("input_particles.txt","r");
	if (cfilex == NULL) {
		ops_printf("ERROR: File cannot be read\n");
		exit(EXIT_FAILURE);
	}
	bool flag;


	if (myRank == 0) {
		fscanf(cfilex, "%d\n", &Nparticles);
	}

#ifdef OPS_MPI
	MPI_Bcast(&Nparticles, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
#endif

	if (Nparticles < 1)
		Nsize = 1;
	else
		Nsize = Nparticles;

	Real xTmp[Nsize], yTmp[Nsize], zTmp[Nsize], radTmp[Nsize];
	Real uTmp[Nsize], vTmp[Nsize], wTmp[Nsize];
	Real oxTmp[Nsize], oyTmp[Nsize], ozTmp[Nsize];

	if (Nparticles > 0) {

		if (myRank == 0)
			for (int iPar = 0; iPar < Nparticles; iPar++) {
				fscanf(cfilex,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						&xTmp[iPar], &yTmp[iPar], &zTmp[iPar], &uTmp[iPar],
						&vTmp[iPar], &wTmp[iPar], &oxTmp[iPar], &oyTmp[iPar],
						&ozTmp[iPar], &radTmp[iPar]);

		}

#ifdef OPS_MPI
		MPI_Bcast(xTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(yTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(zTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(uTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(vTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(wTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(oxTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(oyTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(ozTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(radTmp, Nparticles, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
#endif

	}

	for (int iDir = 0; iDir < Nparticles; iDir++)
		printf("Rank %d: Particle %d of %d: [%f %f %f] Radius = %f\n", ops_get_proc(), iDir, Nparticles,
				xTmp[iDir], xTmp[iDir], zTmp[iDir], radTmp[iDir]);

	//Assign particles to blocks
	AssignParticlesToBlocksSpheres(Nparticles, xTmp, yTmp,
			zTmp, radTmp, uTmp,vTmp,  wTmp,  oxTmp,oyTmp, ozTmp);
	fclose(cfilex);
}
