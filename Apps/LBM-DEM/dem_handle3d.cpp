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

/*!
 * @brief   Functions for particle handling in the LBM code
 * @author  Chrysovalantis Tsigginos
 * @details General functions for direct DEM-LBM simulations
 * 			User must define functions for fluid-particle interaction scheme
 */
#include "dem_handle3d.h"


#include "box_handling.h"
#include "mui_interface.h"

using namespace std;
class muiInterface;



Real dtDEM;
Real cutoff_distance;
int flag_vel;
int Nfl, Npl;
Real alpha;

int muiInterfaceFlag;

muiInterface* interface;

int restartFlag = 0;
int restartTimeStep;

void SetupDEMLBMcoupling() {

	Real skin = 0.001;
	cutoff_distance = skin;
	flag_vel = 0;
}

Real DemTimeStep() {

	return dtDEM;
}

void SetDemTimestep(Real timestep) {

	dtDEM = timestep;
}

void AllocateDomainVariables() {

	int blockNum = BlockNum();
	StenList = new int[2 * SPACEDIM];

	muiInterfaceFlag = 0;

    allocateBlockVariables(blockNum);

}

void AllocateParticles() {
	int blockNum1 = BlockNum();

	particleAllocation(blockNum1);


	AllocateFluidParticleInteractionVeriables();
}

void SetupParticles() {

	AllocateDomainVariables();
	AllocateParticles();
	SetupDEMLBMcoupling();
	FindingBoxBound();

}




//Auxiliary general functions for particle mapping into fluid domain

int CheckDistance() {

	int flag;
	int localFlag = 0;
	Real distance;
	Real cutOffSq = cutoff_distance * cutoff_distance;
	Real dx, dy, dz, dx1;
	dx1 = dxLBM();

	//printf("Rank %d: CutOffSq = %f and number of particles %d\n",ops_get_proc(), cutOffSq, Nparticles[0]);
	localFlag = 1;
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		if (Nparticles[blockIndex]==0)
			localFlag = 0;

		for (int particleIndex = 0; particleIndex < Nparticles[blockIndex]; particleIndex++) {
			dx = xp[blockIndex][particleIndex] - xp_old[blockIndex][particleIndex];
			dy = yp[blockIndex][particleIndex] - yp_old[blockIndex][particleIndex];
			dz = zp[blockIndex][particleIndex] - zp_old[blockIndex][particleIndex];

			distance = dx * dx + dy * dy + dz * dz;
			if (distance > cutOffSq * dx1 * dx1) {
				localFlag = 1;
				break;
			}
			else
				localFlag = 0;
		//	printf("Rank %d: Particle %d Old [%f %f %f] new [%f %f %f]\n", ops_get_proc(), particleIndex, xp_old[blockIndex][particleIndex],yp_old[blockIndex][particleIndex], zp_old[blockIndex][particleIndex],
		//																   xp[blockIndex][particleIndex], yp[blockIndex][particleIndex], zp[blockIndex][particleIndex]);
		 }

	}
	flag = localFlag;
#ifdef OPS_MPI
	MPI_Allreduce(&localFlag, &flag, 1, MPI_INT, MPI_MAX, OPS_MPI_GLOBAL);
#endif

	return flag;
}


void UpdateOldParticleLocation() {

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int iPar = 0; iPar < Nparticles[blockIndex]; iPar++) {
			xp_old[blockIndex][iPar] = xp[blockIndex][iPar];
			yp_old[blockIndex][iPar] = yp[blockIndex][iPar];
			zp_old[blockIndex][iPar] = zp[blockIndex][iPar];
		}

	}
}

void InitializeParticleOldLocation() {

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int iPar  = 0; iPar < Nparticles[blockIndex]; iPar++) {
			xp_old[blockIndex][iPar] = 0.0;
			yp_old[blockIndex][iPar] = 0.0;

			zp_old[blockIndex][iPar] = 0.0;

		}
	}
}

//Particle stencil
void ParticleMapping() {

	int ix, iy, iz;
	int iminl, imaxl, jminl, jmaxl, kminl, kmaxl;
	int ap;
	double dx = dxLBM();
	double x_m[BlockNum()], y_m[BlockNum()], z_m[BlockNum()];

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		x_m[blockIndex] = x_bounds[2 * SPACEDIM * blockIndex];
		y_m[blockIndex] = x_bounds[2 * SPACEDIM * blockIndex + 2];
		z_m[blockIndex] = x_bounds[2 * SPACEDIM * blockIndex + 4];


	}

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		//Finding particle stencil
		int nlocal = Nparticles[blockIndex] + Nperiodic[blockIndex];
		for (int iPar = 0; iPar < nlocal; iPar++) {
			ap = (int) ceil(Radius[blockIndex][iPar] / dx);

			ix = (int) floor((xp[blockIndex][iPar] - x_m[blockIndex]) / dx);
			iy = (int) floor((yp[blockIndex][iPar] - y_m[blockIndex]) / dx);

			iz = (int) floor((zp[blockIndex][iPar] - z_m[blockIndex]) / dx);

			iminl = - ap;
			imaxl = ap + 1;
			jminl = - ap;
			jmaxl = ap + 1;


			kminl = -ap;
			kmaxl = ap + 1;

			//Shift to global system
			iminl += ix;
			imaxl += ix;
			jminl += iy;
			jmaxl += iy;

			kminl += iz;
			kmaxl += iz;


			//Verify that particle stencil do not exceed the sub-domain.

			if (iminl < Nf[2 * SPACEDIM * blockIndex])
				iminl = Nf[2 * SPACEDIM * blockIndex];

			if (jminl < Nf[2 * SPACEDIM * blockIndex +2])
				jminl = Nf[2 * SPACEDIM * blockIndex +2];

			if (imaxl >= Nf[2 * SPACEDIM * blockIndex + 1])
				imaxl = Nf[2 * SPACEDIM * blockIndex + 1];

			if (jmaxl >= Nf[2 * SPACEDIM * blockIndex + 3])
				jmaxl = Nf[2 * SPACEDIM * blockIndex + 3];



			if (kminl < Nf[2 * SPACEDIM * blockIndex + 4])
				kminl = Nf[2 * SPACEDIM * blockIndex + 4];

			if (kmaxl >= Nf[2 * SPACEDIM * blockIndex + 5])
				kmaxl =  Nf[2 * SPACEDIM * blockIndex + 5];

			ParticleList[blockIndex][2 * SPACEDIM * iPar] = iminl;
			ParticleList[blockIndex][2 * SPACEDIM * iPar + 1] = imaxl;
			ParticleList[blockIndex][2 * SPACEDIM * iPar + 2] = jminl;
			ParticleList[blockIndex][2 * SPACEDIM * iPar + 3] = jmaxl;

			ParticleList[blockIndex][2 * SPACEDIM * iPar + 4] = kminl;
			ParticleList[blockIndex][2 * SPACEDIM * iPar + 5] = kmaxl;

/*			ops_printf("Particle %d stencil: [%d %d %d %d %d %d]\n",iPar, ParticleList[blockIndex][2 * SPACEDIM * iPar + 0],
					    ParticleList[blockIndex][2 * SPACEDIM * iPar + 1], ParticleList[blockIndex][2 * SPACEDIM * iPar + 2],
					    ParticleList[blockIndex][2 * SPACEDIM * iPar + 3], ParticleList[blockIndex][2 * SPACEDIM * iPar + 4],
						ParticleList[blockIndex][2 * SPACEDIM * iPar + 5]);
*/
		}
	}
}


void CalculateDragForce(Real dt, SizeType iter) {

	int nMax;
	nMax = 6;


	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int iPar = 0; iPar < Nparticles[blockIndex]; iPar++) {
			for (int index = 0 ; index < nMax; index++)
				FDrag[blockIndex][nMax * iPar + index] /= dt;
		}

	}

//	if (iter % 10 == 0)
//		for (int iPar = 0 ; iPar < Nparticles[0]; iPar++) {
//			printf("Rank %d at iter %d: Particle %d [%f %f %f]: Fd = [%14.9e  %14.9e %14.9e] and Md = [%14.9e %14.9e %14.9e]\n", ops_get_proc(),iter, iPar, xp[0][iPar], yp[0][iPar], 
//zp[0][iPar],
//															    FDrag[0][nMax * iPar], FDrag[0][nMax * iPar + 1], FDrag[0][nMax * iPar + 
//2], FDrag[0][nMax * iPar + 3], FDrag[0][nMax * iPar + 4], FDrag[0][nMax * iPar + 5]);

//	}

}

void InitializeDragForce() {
	int nMax;
	if (SPACEDIM == 2)
		nMax = 3;
	else if (SPACEDIM == 3)
		nMax = 6;

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int iPar = 0; iPar < Nparticles[blockIndex]; iPar++)
			for (int index  =0; index < nMax; index++)
				FDrag[blockIndex][nMax * iPar + index] = 0;

	}
}

void ReadParticleDetails() {


	FILE *cfilex;
	cfilex = fopen("input_particles.txt","r");
	if (cfilex == NULL) {
		ops_printf("ERROR: File cannot be read\n");
		exit(EXIT_FAILURE);
	}
	bool flag;
	int Np;

	if (ops_get_proc()==0)
		fscanf(cfilex,"%d\n",&Np);


	MPI_Bcast(&Np, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
	int Np1;
	if (Np < 1) {
		Np1 = 1;
	}
	else
		Np1 = Np;

	Real xTmp[Np1], yTmp[Np1], zTmp[Np1], upTmp[Np1], vpTmp[Np1], wpTmp[Np1], omegaXTmp[Np1], omegaYTmp[Np1], omegaZTmp[Np1], radTmp[Np1];

	if (ops_get_proc()==0) {
		for (int iPar = 0; iPar < Np; iPar++)
		fscanf(cfilex,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &xTmp[iPar], &yTmp[iPar], &zTmp[iPar], &upTmp[iPar], &vpTmp[iPar], &wpTmp[iPar],
				   &omegaXTmp[iPar], &omegaYTmp[iPar], &omegaZTmp[iPar], &radTmp[iPar]);
	}

	if (Np > 0)  {
		MPI_Bcast(xTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(yTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(zTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(upTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(vpTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(wpTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(omegaXTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(omegaYTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(omegaZTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
		MPI_Bcast(radTmp, Np, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
	}



	for (int iPar = 0; iPar < Np; iPar++) {


		for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {


			flag = SphereParallepidedIntersection(xTmp[iPar], yTmp[iPar], zTmp[iPar], radTmp[iPar], xBoundLocal[2 * SPACEDIM * blockIndex], xBoundLocal[2 * SPACEDIM * blockIndex + 1],
                    							  xBoundLocal[2 * SPACEDIM * blockIndex+2], xBoundLocal[2 * SPACEDIM * blockIndex + 3], xBoundLocal[2 * SPACEDIM * blockIndex + 4],
												  xBoundLocal[2 * SPACEDIM * blockIndex + 5]);



			if (flag == true) {

				Nparticles[blockIndex] += 1;
				if (Nmax[blockIndex] <= Nparticles[blockIndex]) {
					ReAllocateMemory(1000, blockIndex);
				}
				int Np1 = Nparticles[blockIndex] - 1;
				xp[blockIndex][Np1] = xTmp[iPar];
				yp[blockIndex][Np1] = yTmp[iPar];
				Radius[blockIndex][Np1] = radTmp[iPar];

				up[blockIndex][Np1] = upTmp[iPar];
				vp[blockIndex][Np1] = vpTmp[iPar];
				omegaZ[blockIndex][Np1] = omegaZTmp[iPar];

				zp[blockIndex][Np1] = zTmp[iPar];

				wp[blockIndex][Np1] = wpTmp[iPar];
				omegaX[blockIndex][Np1] = omegaXTmp[iPar];
				omegaY[blockIndex][Np1] = omegaYTmp[iPar];

//				printf("Rank %d: Particle  of radius %f inserted at location [%f %f %f] with initial vel [%f %f %f] and om=[%f %f %f] \n",
//						ops_get_proc(), Radius[blockIndex][Np1], xp[blockIndex][Np1], yp[blockIndex][Np1], zp[blockIndex][Np1],
//						up[blockIndex][Np1], vp[blockIndex][Np1], wp[blockIndex][Np1],
//						omegaX[blockIndex][Np1], omegaY[blockIndex][Np1],omegaZ[blockIndex][Np1]);
//

			}
		}



	}

}

void 	UpdateParticlMappingDragForce(SizeType timeSetup) {

	double dt = TimeStep();

	//Initialize the Drag force

	InitializeDragForce();

	//map the particles
	if (periodicFlag == 1) {
		PeriodicExchange();
		InitializePeriodicDragForce();
	}

	ParticleMapping();

    EstimateDragForce2ndOrder3D();
    if (periodicFlag == 1)
    	ReverseComm();

    CalculateDragForce(dt, timeSetup);




}

void InitializeDEMLBMRestart(SizeType &maxStep, SizeType timeSetup,SizeType RunStep, int savingFlag) {

	int myRank;;
	int tempPerFlags[SPACEDIM];
	Real xper[SPACEDIM], xcutOff[SPACEDIM];
	long int maxStep1;

#ifdef OPS_MPI
	myRank = ops_get_proc();
#endif

	if (muiInterfaceFlag == 1) {
		long int timeSetup1 = timeSetup;
		interface->extractData(timeSetup1, maxStep1, alpha, tempPerFlags);

		if (tempPerFlags[0] != periodic[0]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in x-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (tempPerFlags[1] != periodic[1]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in y-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (tempPerFlags[2] != periodic[2]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in z-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (periodicFlag==1) {
			interface->extractDataPeriodic(xper, xcutOff);
			SetPeriodicSize3D(xper[0], xper[1], xper[2]);
			SetCutOff3D(xcutOff[0], xcutOff[1], xcutOff[2]);
		}

		ExtractParticleData(timeSetup);

		maxStep = maxStep1;
	}
	else {
		alpha  = 1;
		maxStep = RunStep + timeSetup;
		ops_printf("Maxstep %d\n", maxStep);
		ReadParticleDetails();

		if (periodicFlag == 1) {
			Real xper1, yper1, zper1;

			if (periodic[0]==1)
				 xper1 = BlockSize(0)[0] + dxLBM();
			else
				xper1 = BlockSize(0)[0];

			if (periodic[1]==1)
				yper1 = BlockSize(0)[1] + dxLBM();
			else
				yper1 = BlockSize(0)[1];

			if (periodic[2]==1)
				 zper1 = BlockSize(0)[2] + dxLBM();
			else
				zper1 = BlockSize(0)[3];

			SetPeriodicSize3D(xper1, yper1, zper1);

			Real xCutOff = dxLBM();
			Real yCutOFf = dxLBM();
			Real zCutOff = dxLBM();

			SetCutOff3D(xCutOff, yCutOFf, zCutOff);

			ops_printf("I initialize distancres\n");
		}


	}

	UpdateOldParticleLocation();

	if (alpha < 0.000001) {
		ops_printf("alpha is equal to %f\n",  alpha);
		exit(EXIT_FAILURE);
	}


	if (alpha <= 1) {
		Npl = static_cast<int>(1.0/alpha);
		Nfl = 1;
	}
	else {
		Npl = 1;
		Nfl = static_cast<int>(alpha);
	}

	SetDemTimestep(alpha * TimeStep());


	if (periodicFlag == 1)
		PeriodicPartition();


	if (timeSetup == 0 ) {
		Real convergenceRate = 1e-7;
		int maxIters = 454;
		IterateDEMLBMSS(convergenceRate, 10000, maxIters,savingFlag); //fix the values

		if (muiInterfaceFlag == 1)
			SendParticleData(timeSetup);
	}
	else {
		UpdateParticlMappingDragForce(timeSetup);

		if (muiInterfaceFlag ==1)
			SendParticleData(timeSetup);
	}




}

void InitializeDEMLBM(Real convergenceRate, int maxIters, int checkPeriod,SizeType &maxStep,int savingFlag) {

	int myRank, timeSetup;
	int tempPerFlags[SPACEDIM];
	Real xper[SPACEDIM], xcutOff[SPACEDIM];
#ifdef OPS_MPI
	myRank = ops_get_proc();
#endif

	if (muiInterfaceFlag == 1) {

		long int timeSetup = 0;
		long int maxStep1;
		interface->extractData(timeSetup, maxStep1, alpha, tempPerFlags);

		maxStep = maxStep1;

		if (tempPerFlags[0] != periodic[0]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in x-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (tempPerFlags[1] != periodic[1]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in y-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (tempPerFlags[2] != periodic[2]) {
			ops_printf("ERROR MP-LBM: Inconsistent B.C in z-direction between in DEM-LBM coupled simulations\n");
			exit(EXIT_FAILURE);
		}

		if (periodicFlag==1) {
			interface->extractDataPeriodic(xper, xcutOff);
			SetPeriodicSize3D(xper[0], xper[1], xper[2]);
			SetCutOff3D(xcutOff[0], xcutOff[1], xcutOff[2]);
		}


		ExtractParticleData(0);


//		ops_printf("Rank  %d: %d particles extracted\n",ops_get_proc(), Nparticles[0]);
//		for (int iPar = 0; iPar < Nparticles[0]; iPar++) {
//			printf("Rank %d: %d [%f %f %f] of radius %f\n", ops_get_proc(), iPar, xp[0][iPar], yp[0][iPar], zp[0][iPar],Radius[0][iPar]);
//		}

		if (restartFlag == 1) {
			UpdateOldParticleLocation();
		}
		else
			InitializeParticleOldLocation();

		ops_printf("I Initialze MUI particles\n");

	}
	else {
		alpha  = 1;
		maxStep = 10000;
		ReadParticleDetails();
		InitializeParticleOldLocation();



		if (periodicFlag == 1) {

			Real xper1 = BlockSize(0)[0] + dxLBM();
			Real yper1 = BlockSize(0)[1] + dxLBM();
			Real zper1 = BlockSize(0)[2] + dxLBM();

			SetPeriodicSize3D(xper1, yper1, zper1);

			Real xCutOff = dxLBM();
			Real yCutOFf = dxLBM();
			Real zCutOff = dxLBM();

			SetCutOff3D(xCutOff, yCutOFf, zCutOff);

			ops_printf("I initialize distancres\n");
		}
	}

	if (alpha < 0.000001) {
		ops_printf("alpha is equal to %f\n",  alpha);
		exit(EXIT_FAILURE);
	}


	if (alpha <= 1) {
		Npl = static_cast<int>(1.0/alpha);
		Nfl = 1;
	}
	else {
		Npl = 1;
		Nfl = static_cast<int>(alpha);
	}

	SetDemTimestep(alpha * TimeStep());


	if (periodicFlag == 1)
		PeriodicPartition();




	ops_printf("Particle set. Ready to iterate\n");
	ops_printf("restartFlag = %d\n", restartFlag);

	ops_printf("convergence Rate=% f checkPeriod = %d, maxIters = %d\n", convergenceRate, checkPeriod, maxIters);
	IterateDEMLBMSS(convergenceRate, checkPeriod, maxIters, savingFlag); //fix the values

	if (muiInterfaceFlag == 1)
		SendParticleData(0);


}
//MUI FUNCTIONS

void DestroyParticleParams() {

	free2d(xp);
	free2d(yp);
	free2d(up);
	free2d(vp);
	free2d(omegaZ);
	free2d(Radius);



	free2d(xp_old);
	free2d(yp_old);


	free2d(zp);
	free2d(wp);
	free2d(omegaX);
	free2d(omegaY);
	free2d(zp_old);


	free2d(FDrag);

	free(Nparticles);
	free(Nmax);
	free(Nperiodic);

	free2d(ParticleList);
	delete[] StenList;



	deleteBlockVariables();


}

void DefineMUIInterface(Real Rmax) {

	ops_printf("MP_LBM: Ready to define uniface\n");

	//new mui::uniface2d("mpi://LBM/ifs");

	interface = new muiInterface(Rmax);


	//interface = new mui::uniface3d("mpi://LBM/ifs");
	ops_printf("MP-LBM: Defined uniface3d since %dD\n", SPACEDIM);



	ops_printf("MP-LBM: Defined uniface %d\n", SPACEDIM);
}

void SetupMUICommunication(int maxIter) {

	interface->setDomains(maxIter);
}

void DestroyMUIInterface() {


	delete interface;



}

//#ifdef OPS_3D
void ExtractParticles3D(int time) {

	vector<Real> xpTemp;
	vector<Real> ypTemp;
	vector<Real> zpTemp;

	vector<Real> radius;
	vector<Real> upTemp;
	vector<Real> vpTemp;
	vector<Real> wpTemp;
	vector<Real> oxTemp;
	vector<Real> oyTemp;
	vector<Real> ozTemp;


	//ops_printf("Time to extract: %d\n",time);
	bool flag;
	int nsd = 2 * SPACEDIM;

	//for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
//#ifdef OPS_MPI
//		sub_block_list sb = OPS_sub_block_list[blockIndex];
//		if (!sb->owned) continue;
//#endif


	interface->updateParticles(time, xpTemp, ypTemp, zpTemp, radius, upTemp, vpTemp, wpTemp, oxTemp, oyTemp, ozTemp);

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		Nparticles[blockIndex] = 0;
	}

	for (int iVec = 0; iVec < xpTemp.size(); ++iVec) {

		Real xtemp = xpTemp[iVec];
		Real ytemp = ypTemp[iVec];
		Real ztemp = zpTemp[iVec];
		Real radTmp = radius[iVec];



		for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
			sub_block_list sb = OPS_sub_block_list[blockIndex];
			if (!sb->owned) continue;
#endif
			flag = SphereParallepidedIntersection(xtemp,ytemp,ztemp, radTmp,
												  xBoundLocal[nsd * blockIndex], xBoundLocal[nsd * blockIndex + 1],
												  xBoundLocal[nsd * blockIndex + 2], xBoundLocal[nsd * blockIndex + 3],
												  xBoundLocal[nsd * blockIndex + 4], xBoundLocal[nsd * blockIndex + 5]);

			if (flag == true) {
				Nparticles[blockIndex] += 1;
				if (Nmax[blockIndex] <= Nparticles[blockIndex])
					ReAllocateMemory(1000, blockIndex);


				int Np = Nparticles[blockIndex] - 1;
				xp[blockIndex][Np] = xtemp;
				yp[blockIndex][Np] = ytemp;
				zp[blockIndex][Np] = ztemp;

				Radius[blockIndex][Np] = radTmp;
				up[blockIndex][Np] = upTemp[iVec];
				vp[blockIndex][Np] = vpTemp[iVec];
				wp[blockIndex][Np] = wpTemp[iVec];

				omegaX[blockIndex][Np] = oxTemp[iVec];
				omegaY[blockIndex][Np] = oyTemp[iVec];
				omegaZ[blockIndex][Np] = ozTemp[iVec];
			}
		}

	}

	//ops_printf("I extracted %d particles\n", Nparticles[0]);

}



void ExtractParticleData(SizeType time) {


		int time1 = time;

		//ops_printf("Entering Extract Particles 3D\n");
		ExtractParticles3D( time1);


}

void SendParticleData3D(int t) {

 interface->sendParticles(t);

}


void SendParticleData(SizeType timestep) {


	int t = timestep;

	SendParticleData3D(t);

}

void ForgetData(int t) {



}

//Auxiliary functions for particle mapping


bool SphereParallepidedIntersection(Real cx, Real cy, Real cz, Real radius, Real xbMin,
									Real  xbMax, Real ybMin, Real ybMax, Real zbMin, Real zbMax) {

	Real testX = cx;
	Real testY = cy;
	Real testZ = cz;

	if (cx < xbMin)
		testX = xbMin;
	else if (cx > xbMax)
		testX = xbMax;

	if (cy < ybMin)
		testY = ybMin;
	else if (cy > ybMax)
		testY = ybMax;

	if (cz < zbMin)
		testZ = zbMin;
	else if (cz > zbMax)
		testZ = zbMax;

	Real distX = cx - testX;
	Real distY = cy - testY;
	Real distZ = cz - testZ;

	Real distance = distX * distX + distY * distY + distZ * distZ;

	if (distance <= radius * radius)
		return true;

	return false;

}
