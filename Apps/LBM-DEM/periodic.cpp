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
 * @brief   Functions for the implementation of periodic conditions at DEM-LBM simulations
 * @author  Chrysovalantis Tsigginos
 * @details Wrap and Kernel functions for the implementation of periodic B.C
 */


#include "periodic.h"
#include <stdlib.h>
#include <stdio.h>

#define BIG 1.0e8


int periodicFlag = 0;
int periodic[3];

Real* cutDEM;
Real* xper; //Size of periodic domain in all directions

int** procSend; //Proc that each process sends to
int** procRecv;//Proce that each process receives form
int** pbcFlag;// Flag for sending data through the periodic BC

Real** slablo; //lower bound of slab that will be send in each swap
Real** slabhi; //upper bound of slab that will be send in each swap

int*** sendList; //List of particles that will be transmitted to the periodic rank
int** maxSend; //maximum number of particles that each rank send  at each swap
int** sendNum; //Number of particles send at each swap
int** recvNum; //Number of particles received at each swap
int** firstRecv; //First particle that received
int** reverseSendNum;
int** reverseRecvNum;

Real sizeForward; //size of forward data send per particle
Real sizeReverse; //size of reverse data send per particle

int maxSending;
int maxReceiving;

Real* buffSend; //buffer for sending
Real* buffRecv; //buffer for receiving

ops_halo* g_HalosX; //Halos in the x-directions
ops_halo_group g_HaloGroupsX;

ops_halo* g_HalosY; //Halos in the y-direction
ops_halo_group g_HaloGroupsY;

ops_halo* g_HalosZ;
ops_halo_group g_HaloGroupsZ;


void  DefinePeriodicHaloTransfer() {

	int haloDepth, dp[SPACEDIM], dm[SPACEDIM];
	int baseTo[SPACEDIM], baseFrom[SPACEDIM], haloSize[SPACEDIM], dir[SPACEDIM];
	int nx, ny, nz;
	int haloCrt;
	int nHalosX = 0;
	int nHalosY = 0;
	int nHalosZ = 0;

	if (periodic[0] ==  1) {
		nHalosX = 2;
		g_HalosX = new ops_halo[2];
	}

	if (periodic[1] == 1) {
		nHalosY = 2;
		g_HalosY  = new ops_halo[2];
	}

#ifdef OPS_3D
	if (periodic[2] == 1) {
		nHalosZ  = 2;
		g_HalosZ = new ops_halo[2];
	}
#endif

	if (BlockNum()>1) {
		ops_printf("MP-LBM: ERROR: Periodic scheme requires 1 block \nEXITING\n");
		exit(EXIT_FAILURE);
	}

	haloDepth = HaloPtNum();
	ops_printf("HaloDepth = %d\n", haloDepth);
	nx = BlockSize(0)[0];
	dir[0] = 1;
	ny = BlockSize(0)[1];
	dir[1] = 2;

	for (int idim = 0; idim < SPACEDIM; idim++) {
		dp[idim] = haloDepth;
		dm[idim] = - haloDepth;
	}

#ifdef OPS_3D
	nz = BlockSize(0)[2];
	dir[2] = 3;
#endif

	for (int idim = 0; idim < SPACEDIM; idim++) {
		if (periodic[idim] == 1) {

			if(idim == 0) {
				baseFrom[0] = 0;
				baseFrom[1] = 0;

				baseTo[0] = nx;
				baseTo[1] = 0;

				haloSize[0] = 1;
				haloSize[1] = ny;

#ifdef OPS_3D
				baseTo[2] = 0;
				baseFrom[2] = 0;
				haloSize[2] = nz;
#endif

				g_HalosX[0] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);

				//Setup the 2nd halo
				baseFrom[0] = nx-1;
				baseTo[0] = dm[0];

				g_HalosX[1] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);

				g_HaloGroupsX = ops_decl_halo_group(2, g_HalosX);
#if DebugLevel >= 1
				ops_printf("Defined x-periodic halos\n");
#endif
			}
			else if (idim == 1) {

				baseFrom[0] = dm[0];
				baseFrom[1] = 0;

				baseTo[0] = dm[0];
				baseTo[1] = ny;

				haloSize[0] = nx+ dp[0] -dm[0];
				haloSize[1] = 1;

#ifdef OPS_3D
				baseTo[2] = 0;
				baseFrom[2] = 0;
				haloSize[2] = nz;
#endif

				g_HalosY[0] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);
				baseFrom[1] = ny+dm[1];
				baseTo[1] = dm[1];

				g_HalosY[1] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);

				g_HaloGroupsY = ops_decl_halo_group(2, g_HalosY);
			}
#ifdef OPS_3D
			else if (idim ==2) {
				baseFrom[0] = dm[0];
				baseFrom[1] = dm[1];
				baseFrom[2] = 0;

				baseTo[0] = dm[0];
				baseTo[1] = dm[1];
				baseTo[2] = nz;

				haloSize[0] = nx + dp[0] - dm[0];
				haloSize[1] = ny + dp[1] - dm[1];
				haloSize[2] = 1;

				g_HalosZ[0] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);

				baseFrom[2] = nz+dm[1];
				baseTo[2] = dm[1];

				g_HalosZ[1] = ops_decl_halo(g_fStage[0], g_fStage[0], haloSize, baseFrom, baseTo, dir, dir);
				g_HaloGroupsZ = ops_decl_halo_group(2, g_HalosZ);
			}
#endif

		}

	}

}

void DestroyPeriodic() {

	free(cutDEM);

	free2d(procSend);
	free2d(procRecv);
	free2d(pbcFlag);

	free2d(maxSend);
	free2d(slabhi);
	free2d(slablo);

	free2d(sendNum);
	free2d(recvNum);
	free2d(firstRecv);
	free2d(reverseSendNum);
	free2d(reverseRecvNum);

	free(buffRecv);
	free(buffSend);
	free(xper);
}

void PeriodicHaloTransfer() {

	if (periodic[0]==1)
		ops_halo_transfer(g_HaloGroupsX);

	if (periodic[1]==1)
		ops_halo_transfer(g_HaloGroupsY);

#ifdef OPS_3D
	if (periodic[2]==1)
		ops_halo_transfer(g_HaloGroupsZ);
#endif

}



void ImplementBoundary3DPeriodic() {

	for (auto boundary : BlockBoundaries()) {

		int *range{BoundarySurfaceRange(boundary.blockIndex,
                boundary.boundarySurface)};

		if (boundary.boundaryScheme != BoundaryScheme::Periodic) {
			 TreatBlockBoundary3D(boundary.blockIndex, boundary.componentID,
			                             boundary.givenVars.data(), range,
			                             boundary.boundaryScheme, boundary.boundarySurface);

		}
	}
}

void InitializePeriodic() {

	int blockNum = BlockNum();

	cutDEM = (Real *) malloc(SPACEDIM * sizeof(Real));
	int Nswap;

	procSend = (int **) malloc(blockNum * sizeof(int *));
	procRecv = (int **) malloc(blockNum * sizeof(int *));
	pbcFlag  = (int **) malloc(blockNum * sizeof(int *));

	slablo = (Real **) malloc(blockNum * sizeof(Real *));
	slabhi = (Real **) malloc(blockNum * sizeof(Real *));

	sendList = (int ***) malloc(blockNum * sizeof(int **));
	maxSend = (int **) malloc(blockNum * sizeof(int *));

	sendNum = (int **) malloc(blockNum * sizeof(int *));
	recvNum = (int **) malloc(blockNum * sizeof(int *));
	firstRecv = (int **) malloc(blockNum * sizeof(int *));

	reverseSendNum = (int **) malloc(blockNum * sizeof(int *));
	reverseRecvNum = (int **) malloc(blockNum * sizeof(int *));


	for (int blockIndex = 0; blockIndex < blockNum; blockIndex++) {

		allocateMemory(procSend[blockIndex], 2 * SPACEDIM, 1);
		allocateMemory(procRecv[blockIndex], 2 * SPACEDIM, 1);

		Nswap = 2 * SPACEDIM;

		allocateMemory(maxSend[blockIndex], Nswap, 1);
		allocateMemory(pbcFlag[blockIndex], Nswap, SPACEDIM);
		allocateMemory(slabhi[blockIndex], Nswap , SPACEDIM);
		allocateMemory(slablo[blockIndex], Nswap, SPACEDIM);

		allocateMemory(sendNum[blockIndex], Nswap, 1);
		allocateMemory(recvNum[blockIndex], Nswap, 1);
		allocateMemory(firstRecv[blockIndex], Nswap, 1);
		allocateMemory(reverseSendNum[blockIndex], Nswap, 1);
		allocateMemory(reverseRecvNum[blockIndex], Nswap, 1);

		sendList[blockIndex] = (int **) malloc (Nswap * sizeof(int *));
		for (int iswap = 0; iswap < Nswap; iswap++) {
			allocateMemory(sendList[blockIndex][iswap],100, 1);
			maxSend[blockIndex][iswap] = 100;
		}
	}

	if (SPACEDIM == 2) {
		sizeForward = 6;
		sizeReverse = 3;
	}
	else if (SPACEDIM == 3) {
		sizeForward = 10;
		sizeReverse = 6;
	}
	else {
			ops_printf("MP-LBM ERROR: Inconsistent dimension\n");
			exit(EXIT_FAILURE);
	}

	maxSending = sizeForward * maxSend[0][0];
	allocateMemory(buffSend, maxSending, 1);
	maxReceiving = sizeForward * maxSend[0][0];
	allocateMemory(buffRecv, maxReceiving, 1);

	allocateMemory(xper, SPACEDIM, 1);

}

void DefinePeriodicBoundariesRestart(std::vector<BoundaryScheme>& boundaryType,std::vector<int> &dir) {

	int direction;

	for (int iDir = 0; iDir < SPACEDIM; iDir++) {

		direction = dir[iDir];
		if ((boundaryType[2*iDir] == BoundaryScheme::Periodic) && (boundaryType[2 * iDir + 1]== BoundaryScheme::Periodic))
			periodic[direction] = 1;
		else if ((boundaryType[2 * iDir] == BoundaryScheme::Periodic) && (boundaryType[2 * iDir + 1] != BoundaryScheme::Periodic)) {
			ops_printf("MP-LBM3D: Inconsistent %d periodic boundary condition\n", direction);
			exit(EXIT_FAILURE);
		}
		else if ((boundaryType[2 * iDir] != BoundaryScheme::Periodic) &&(boundaryType[2 * iDir + 1] == BoundaryScheme::Periodic)) {
			ops_printf("MP-LBM3D: Inconsistent %d periodic boundary condition\n", direction);
			exit(EXIT_FAILURE);
		}
		else
			periodic[direction] = 0;

	}

	periodicFlag = 0;

	if ((periodic[0]==1) || (periodic[1] == 1) || (periodic[2] == 1))
		periodicFlag = 1;


	ops_printf("Periodic Flag is %d\n", periodicFlag);
	ops_printf("x-Per: %d y-Per: %d z-Per: %d\n", periodic[0], periodic[1], periodic[2]);

	if (periodicFlag==1)
		ops_printf("x-Per: %d y-Per: %d z-Per: %d\n", periodic[0], periodic[1], periodic[2]);

	//Define periodic Partition
	if (periodicFlag == 1)
		DefinePeriodicHaloTransfer();

	ops_printf("PeriodicFlag = %d\n",periodicFlag);

}

void DefinePeriodicBoundaries(BoundaryScheme* boundaryType) {

	//Check for periodic Boundary conditions
	if ((boundaryType[0] == BoundaryScheme::Periodic) && (boundaryType[1] == BoundaryScheme::Periodic))
		periodic[0] = 1;
	else if ((boundaryType[0] == BoundaryScheme::Periodic) && (boundaryType[1] != BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent x-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else if ((boundaryType[0] != BoundaryScheme::Periodic) &&(boundaryType[1] == BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent x-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else
		periodic[0] = 0;



	if ((boundaryType[2] == BoundaryScheme::Periodic) && (boundaryType[3] == BoundaryScheme::Periodic))
		periodic[1] = 1;
	else if ((boundaryType[2] != BoundaryScheme::Periodic) && (boundaryType[3] == BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent y-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else if ((boundaryType[2] == BoundaryScheme::Periodic) && (boundaryType[3] != BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent y-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else
		periodic[1] = 0;

	if ((boundaryType[4] == BoundaryScheme::Periodic) && (boundaryType[5] == BoundaryScheme::Periodic))
		periodic[2] = 1;
	else if ((boundaryType[4] != BoundaryScheme::Periodic) && (boundaryType[5] == BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent z-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else if ((boundaryType[4] == BoundaryScheme::Periodic) && (boundaryType[5] != BoundaryScheme::Periodic)) {
		ops_printf("MP-LBM3D: Inconsistent z-periodic boundary condition\n");
		exit(EXIT_FAILURE);
	}
	else
		periodic[2] = 0;

	periodicFlag = 0;

	if ((periodic[0]==1) || (periodic[1] == 1) || (periodic[2] == 1))
		periodicFlag = 1;


	ops_printf("Periodic Flag is %d\n", periodicFlag);
	ops_printf("x-Per: %d y-Per: %d z-Per: %d\n", periodic[0], periodic[1], periodic[2]);

	if (periodicFlag==1)
		ops_printf("x-Per: %d y-Per: %d z-Per: %d\n", periodic[0], periodic[1], periodic[2]);

	//Define periodic Partition
	if (periodicFlag == 1)
		DefinePeriodicHaloTransfer();

	ops_printf("PeriodicFlag = %d\n",periodicFlag);
}


void InitializePeriodicDragForce() {

	int nlocal;
	int sizeFd;
	if (SPACEDIM == 2)
		sizeFd = SPACEDIM + 1;
	if (SPACEDIM == 3)
		sizeFd = 2 * SPACEDIM + 1;

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		nlocal = Nparticles[blockIndex] + Nperiodic[blockIndex];

		for (int iPar = Nparticles[blockIndex]; iPar < nlocal; iPar++) {
			for (int iDir = 0; iDir < sizeFd; iDir++)
				FDrag[blockIndex][sizeFd * iPar + iDir] = 0.0;

		}
	}

}


#ifdef OPS_2D
void SetPeriodicSize2D(Real x1,Real y1) {
	xper[0] = x1;
	xper[1] = y1;
}

void SetCutOff2D(Real xcut, Real ycut) {

	cutDEM[0] = xcut;
	cutDEM[1] = ycut;
}
#endif


#ifdef OPS_3D
void SetPeriodicSize3D(Real x1, Real y1, Real z1) {

	xper[0] = x1;
	xper[1] = y1;
	xper[2] = z1;
}

void SetCutOff3D(Real xcut, Real ycut, Real zcut) {

	cutDEM[0] = xcut;
	cutDEM[1] = ycut;
	cutDEM[2] = zcut;
}
#endif

void PeriodicPartition() {

	int periods[SPACEDIM], coords[SPACEDIM];
	int coordPeriodicLeft[SPACEDIM], coordPeriodicRight[SPACEDIM];
	int size1[SPACEDIM];
	int rankPeriodic;
	Real perNeigh[SPACEDIM];
	int periodicNeighborRight[BlockNum()][SPACEDIM];
	int periodicNeighborLeft[BlockNum()][SPACEDIM];
	Real dx = dxLBM();

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {

#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#else
		ops_printf("OPS-LBM-MUI: Error Require MPI version.\n");
		exit(EXIT_FAILURE);
#endif

#if DebugLevel >= 2
		if (SPACEDIM==2)
			ops_printf("Periodic x: %d Periodic y: %d\n",periodic[0], periodic[1]);
		else if (SPACEDIM == 3)
			ops_printf("Periodic x: %d Periodic y: %d Periodic z: %d\n", periodic[0], periodic[1], periodic[2]);
#endif
		MPI_Cart_get(sb->comm, SPACEDIM, size1, periods, coords);

		if (periodic[0] == 1) { //x-periodic
			if (size1[0] == 1) { //return same process
				MPI_Cart_rank(sb->comm, coords, &rankPeriodic);
				periodicNeighborRight[blockIndex][0] = rankPeriodic;
				periodicNeighborLeft[blockIndex][0] = rankPeriodic;
			}
			else {
				coordPeriodicLeft[1] = coords[1];
				coordPeriodicRight[1] = coords[1];
#ifdef OPS_3D
				coordPeriodicLeft[2] = coords[2];
				coordPeriodicRight[2] = coords[2];
#endif
				if (coords[0] == 0) { //Rank on the left boundary
					coordPeriodicLeft[0] = size1[0]-1;
					MPI_Cart_rank(sb->comm, coordPeriodicLeft, &rankPeriodic);
					periodicNeighborLeft[blockIndex][0] = rankPeriodic;
					periodicNeighborRight[blockIndex][0] = -1;
				}
				else if (coords[0] == size1[0] - 1) {
					coordPeriodicRight[0] = 0;
					MPI_Cart_rank(sb->comm, coordPeriodicRight, &rankPeriodic);
					periodicNeighborLeft[blockIndex][0] = -1;
					periodicNeighborRight[blockIndex][0] = rankPeriodic;
				}
				else {
					periodicNeighborLeft[blockIndex][0] = -1;
					periodicNeighborRight[blockIndex][0] = -1;
				}
			}

		}
		else {
			periodicNeighborLeft[blockIndex][0] = -1;
			periodicNeighborRight[blockIndex][0] =-1;
		}

		if (periodic[1] == 1) {  //y-periodic B.C
			if (size1[1] == 1) {//Return same process on MPI_Comm sb->comm.
				MPI_Cart_rank(sb->comm, coords, &rankPeriodic);
				periodicNeighborRight[blockIndex][1] = rankPeriodic;
				periodicNeighborLeft[blockIndex][1] = rankPeriodic;
			}
			else {
				coordPeriodicLeft[0] = coords[0];
				coordPeriodicRight[0] = coords[0];
#ifdef OPS_3D

				coordPeriodicLeft[2] = coords[2];
				coordPeriodicRight[2] = coords[2];
#endif
				if (coords[1] == 0) { //Rank on left boundary
					coordPeriodicLeft[1] = size1[1] - 1;
					MPI_Cart_rank(sb->comm, coordPeriodicLeft, &rankPeriodic);
					periodicNeighborLeft[blockIndex][1] = rankPeriodic;
					periodicNeighborRight[blockIndex][1] = -1;
				}
				else if (coords[1] == size1[1] - 1) {
					coordPeriodicRight[1] = 0;
					MPI_Cart_rank(sb->comm, coordPeriodicRight, &rankPeriodic);
					periodicNeighborRight[blockIndex][1] = rankPeriodic;
					periodicNeighborLeft[blockIndex][1] = -1;

				}
				else {
					periodicNeighborLeft[blockIndex][1] = -1;
					periodicNeighborRight[blockIndex][1] = -1;
				}

			}

		}
		else {
			periodicNeighborLeft[blockIndex][1] = -1;
			periodicNeighborRight[blockIndex][1] = -1;
		}
#ifdef OPS_3D
		if (periodic[2] == 1) { //z-periodic B.C
			if (size1[2] == 1) {
				MPI_Cart_rank(sb->comm, coords, &rankPeriodic);
				periodicNeighborRight[blockIndex][2] = rankPeriodic;
				periodicNeighborLeft[blockIndex][2] = rankPeriodic;
			}
			else {
				coordPeriodicRight[0] = coords[0];
				coordPeriodicLeft[0] = coords[0];

				coordPeriodicRight[1] = coords[1];
				coordPeriodicLeft[1] = coords[1];

				if (coords[2] == 0) { //Rank on bottom boundary
					coordPeriodicLeft[2] = size1[2] - 1;
					MPI_Cart_rank(sb->comm, coordPeriodicLeft, &rankPeriodic);
					periodicNeighborLeft[blockIndex][2] = rankPeriodic;
					periodicNeighborRight[blockIndex][2] = -1;
				}
				else if (coords[2] == size1[2]-1) {
					coordPeriodicRight[2] = 0;
					MPI_Cart_rank(sb->comm, coordPeriodicRight, &rankPeriodic);
					periodicNeighborLeft[blockIndex][2] = -1;
					periodicNeighborRight[blockIndex][2] = rankPeriodic;

				}
				else { //no periodic process
					periodicNeighborRight[blockIndex][2] = -1;
					periodicNeighborLeft[blockIndex][2] = -1;
				}

			}
		}
		else {
			periodicNeighborRight[blockIndex][2] = -1;
			periodicNeighborLeft[blockIndex][2] = -1;
		}
#endif

		int globalRank = ops_get_proc();
		int localRank;
		MPI_Comm_rank(sb->comm, &localRank);

	/*
	 *		Identifying the periodic swaps
	 *      Two swaps per direction
	 *      1st swap: Receives from left sends to right
	 *      2nd swap: Receives from right sends to left
	 *
	 */

		int nswap = 0;
		for (int idim = 0; idim < SPACEDIM; idim++) {
			for (int iswap = 0; iswap < 2; iswap++) {

			//Verify size of slablo and slabhi
				pbcFlag[blockIndex][nswap * SPACEDIM] = 0;
				pbcFlag[blockIndex][nswap * SPACEDIM + 1] = 0;

				slablo[blockIndex][nswap * SPACEDIM] = xBoundLocal[2 * SPACEDIM * blockIndex];// + 0.5 * dx;
				slablo[blockIndex][nswap * SPACEDIM + 1] = xBoundLocal[2 * SPACEDIM * blockIndex + 2];// + 0.5 * dx;

				slabhi[blockIndex][nswap * SPACEDIM] = xBoundLocal[2 * SPACEDIM * blockIndex + 1];// - 0.5 * dx;
				slabhi[blockIndex][nswap * SPACEDIM + 1] = xBoundLocal[2 * SPACEDIM * blockIndex + 3];// - 0.5 * dx;

#ifdef OPS_3D
				pbcFlag[blockIndex][nswap * SPACEDIM + 2] = 0;

				slablo[blockIndex][nswap * SPACEDIM + 2] = xBoundLocal[2 * SPACEDIM * blockIndex + 4];// + 0.5 * dx;
				slabhi[blockIndex][nswap * SPACEDIM + 2] = xBoundLocal[2 * SPACEDIM * blockIndex + 5];// - 0.5 * dx;
#endif

				procSend[blockIndex][nswap] = -1;
				procRecv[blockIndex][nswap] = -1;


				if (periodic[idim] == 1) {
					if ( (iswap % 2) == 0) {
						procSend[blockIndex][nswap] = periodicNeighborRight[blockIndex][idim];
						procRecv[blockIndex][nswap] = periodicNeighborLeft[blockIndex][idim];
						if (coords[idim] == size1[idim]-1) {//Only this one will send
							pbcFlag[blockIndex][nswap * SPACEDIM + idim] = -1;
							slablo[blockIndex][nswap * SPACEDIM + idim] = slabhi[blockIndex][nswap * SPACEDIM + idim] - cutDEM[idim];
							slabhi[blockIndex][nswap * SPACEDIM + idim] = BIG;
						}

					}
					else {
						procSend[blockIndex][nswap] = periodicNeighborLeft[blockIndex][idim];
						procRecv[blockIndex][nswap] = periodicNeighborRight[blockIndex][idim];
						if (coords[idim] == 0) {
							pbcFlag[blockIndex][nswap * SPACEDIM + idim] = 1;
							slabhi[blockIndex][nswap * SPACEDIM + idim] = slablo[blockIndex][nswap * SPACEDIM + idim] + cutDEM[idim];
							slablo[blockIndex][nswap * SPACEDIM + idim] = -BIG;
						}
					}

				}
				nswap++;
			}

		}

#if DebugLevel >= 2
	printf("\n\n-----------------------------------------------------------------------------------------------------------------\n");
		printf("Local rank %d: Swap 0: Recv %d Send %d, Swap 1: Recv %d Send %d, Swap 2: Recv %d Send %d, Swap 3: Recv %d Send %d Swap 4 Recv %d Send %d, Swap 5 %d %d\n",
				localRank,procRecv[blockIndex][0], procSend[blockIndex][0], procRecv[blockIndex][1], procSend[blockIndex][1],
				           procRecv[blockIndex][2], procSend[blockIndex][2], procRecv[blockIndex][3], procSend[blockIndex][3],
						   procRecv[blockIndex][4], procSend[blockIndex][4], procRecv[blockIndex][5], procSend[blockIndex][5]);

		printf("\n\n Local rank %d: Swap 0: pbcFlag [%d %d %d], Swap 1: pbcFlag [%d %d %d], Swap 2: pbcFlag [%d %d %d], Swap 3: pbcFlag [%d %d %d], Swap 4: pbcFlag [%d %d %d], Swap 5: pbcFlag[%d %d %d]\n",
				localRank, pbcFlag[blockIndex][0 * SPACEDIM], pbcFlag[blockIndex][0 * SPACEDIM + 1], pbcFlag[blockIndex][0 * SPACEDIM + 2],
				           pbcFlag[blockIndex][1 * SPACEDIM], pbcFlag[blockIndex][1 * SPACEDIM + 1], pbcFlag[blockIndex][1 * SPACEDIM + 2],
						   pbcFlag[blockIndex][2 * SPACEDIM], pbcFlag[blockIndex][2 * SPACEDIM + 1], pbcFlag[blockIndex][2 * SPACEDIM + 2],
						   pbcFlag[blockIndex][3 * SPACEDIM], pbcFlag[blockIndex][3 * SPACEDIM + 1], pbcFlag[blockIndex][3 * SPACEDIM + 2],
						   pbcFlag[blockIndex][4 * SPACEDIM], pbcFlag[blockIndex][4 * SPACEDIM + 1], pbcFlag[blockIndex][4 * SPACEDIM + 2],
						   pbcFlag[blockIndex][5 * SPACEDIM], pbcFlag[blockIndex][5 * SPACEDIM + 1], pbcFlag[blockIndex][5 * SPACEDIM + 2]);
		printf("\n\nLocal Rank %d:", localRank);
		for (int iswap = 0; iswap < nswap; iswap++) {
			if (procSend[blockIndex][iswap] > -1) {
				printf("Swap %d: boxSlab [%f %f]x[%f %f] x[%f %f] ",iswap, slablo[blockIndex][iswap * SPACEDIM], slabhi[blockIndex][iswap * SPACEDIM],
																   slablo[blockIndex][iswap * SPACEDIM + 1], slabhi[blockIndex][iswap * SPACEDIM + 1],
																   slablo[blockIndex][iswap * SPACEDIM + 2], slabhi[blockIndex][iswap * SPACEDIM + 2]);
			}
		}

		printf("\n");
		printf("---------------------------------------------------------------------------------------------------------------------\n");


/*
 * NOTES: Function periodicNeighborLeft, periodicNeighborRight
 */
#endif
	}

}

bool decide(Real xp, Real xlo, Real xhi) {

	if ((xlo < xp)&&(xp <= xhi))
		return true;

	return false;
}

void PeriodicExchange() {

	int idim;
	int nfirst, nlast;
	int nsend, nrecv1, nrecv2, nrecv;
	int meLocal;
	Real lo, hi;
	Real xl;
	bool flag;
	Real *buf;
	MPI_Status status;
	MPI_Request request;
	ops_printf("Entering PEriodic Exchange\n");

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		MPI_Comm_rank(sb->comm, &meLocal);
		Nperiodic[blockIndex] = 0;
		idim = -1;
		nlast = 0;

		for (int iswap = 0; iswap < 2 * SPACEDIM; iswap++) {
			//First and second swap transfer only local particles
			//The rest of the swaps include transfers and periodic also.
			if (iswap % 2 ==0) {
				idim++;
				nfirst = 0;
				nlast = Nparticles[blockIndex] + Nperiodic[blockIndex];
#if DebugLevel >= 2
				printf("Rank %d: Swap %d: nlast = %d\n", ops_get_proc(), iswap, nlast);
#endif
			}
			lo = slablo[blockIndex][iswap * SPACEDIM + idim];
			hi = slabhi[blockIndex][iswap * SPACEDIM + idim];

			if (periodic[idim] == 0) {
				sendNum[blockIndex][iswap] = 0;
				recvNum[blockIndex][iswap] = 0;
				reverseSendNum[blockIndex][iswap] = 0;
				reverseRecvNum[blockIndex][iswap] = 0;
				firstRecv[blockIndex][iswap] = -1;
				continue;
			}

			nsend = 0;
			nrecv = 0;

			if (procSend[blockIndex][iswap] > -1) {
				for (int iParticle = nfirst; iParticle < nlast; iParticle++) {
					if (idim == 0)
						xl = xp[blockIndex][iParticle];
					else if (idim == 1)
						xl = yp[blockIndex][iParticle];
#ifdef OPS_3D
					else if (idim == 2)
						xl = zp[blockIndex][iParticle];
#endif
					flag = decide(xl, lo, hi);

					if (flag) {
						if (nsend == maxSend[blockIndex][iswap]) {
							Reallocate(sendList[blockIndex][iswap], nsend+100, 1);
							maxSend[blockIndex][iswap] = nsend + 100;
						}
						sendList[blockIndex][iswap][nsend++] = iParticle;
					}
				}

				if (nsend * sizeForward > maxSending) {
					maxSending = nsend * sizeForward;
					Reallocate(buffSend, sizeForward * nsend, 1);
					Reallocate(buffRecv, sizeForward * nsend, 1);
				}

				if (nsend > 0) {
#ifdef OPS_2D
					packForward2D(blockIndex, nsend, iswap, sendList[blockIndex][iswap], buffSend, periodic, pbcFlag[blockIndex]);
#endif

#ifdef OPS_3D
					packForward3D(blockIndex, nsend, iswap, sendList[blockIndex][iswap], buffSend, periodic, pbcFlag[blockIndex]);
#endif
				}

				if (procSend[blockIndex][iswap] != meLocal)
					MPI_Send(&nsend, 1, MPI_INT, procSend[blockIndex][iswap],0, sb->comm);

			}

			if (procRecv[blockIndex][iswap] > -1) {
				if (meLocal == procRecv[blockIndex][iswap])
					nrecv = nsend;
				else
					MPI_Recv(&nrecv, 1, MPI_INT, procRecv[blockIndex][iswap], 0, sb->comm, &status);
			}

			if ((nrecv > 0) || (nsend > 0)) {
				if (procSend[blockIndex][iswap] == meLocal)
					buf = buffSend;
				else {
					if (nrecv * sizeForward > maxReceiving) {
						maxReceiving = nrecv * sizeForward;
						Reallocate(buffRecv, sizeForward * nrecv, 1);
						Reallocate(buffSend, sizeForward * nrecv, 1);
					}

					if (nrecv)
						MPI_Irecv(buffRecv, nrecv * sizeForward, MPI_DOUBLE, procRecv[blockIndex][iswap], 0,sb->comm, &request);

					if (nsend)
						MPI_Send(buffSend,nsend * sizeForward,MPI_DOUBLE,procSend[blockIndex][iswap],0,sb->comm);

					if (nrecv)
						MPI_Wait(&request, &status);

					buf = buffRecv;
				}

				//Unpack buffer
				if (nrecv) {
#ifdef OPS_2D
					unpackComm2D(blockIndex, nrecv, Nparticles[blockIndex] + Nperiodic[blockIndex], buf);
#endif

#ifdef OPS_3D
					unpackComm3D(blockIndex, nrecv, Nparticles[blockIndex] + Nperiodic[blockIndex], buf);
#endif
				}

			}

			sendNum[blockIndex][iswap] = nsend;
			recvNum[blockIndex][iswap] = nrecv;

			reverseSendNum[blockIndex][iswap] = nrecv;
			reverseRecvNum[blockIndex][iswap] = nsend;

			if (recvNum[blockIndex][iswap] > 0) {
				firstRecv[blockIndex][iswap] = Nparticles[blockIndex] + Nperiodic[blockIndex];
				Nperiodic[blockIndex] += nrecv;
			}
			else {
				firstRecv[blockIndex][iswap] = -1;
			}

		}
	}

#if DebugLevel >= 2
	printf("Rank %d: Actual Particles %d: Periodic particles %d\n", meLocal, Nparticles[0], Nperiodic[0]);
	for (int iswap = 0; iswap < 2 * SPACEDIM; iswap++)
		printf("Rank: %d, Swap: %d, sendNum: %d, recvNum: %d, firstReceive: %d\n",meLocal, iswap, sendNum[0][iswap],
				 recvNum[0][iswap], firstRecv[0][iswap]);

		ops_printf("Printing Particles\n");
		for (int iPar = 0; iPar < Nparticles[0] + Nperiodic[0]; iPar++) {
			if (iPar < Nparticles[0])
				printf("Rank %d, Particle %i (Actual) x = %f y=%f z=%f\n",meLocal, iPar, xp[0][iPar], yp[0][iPar],zp[0][iPar]);
			else
				printf("Rank %d, Particle %i (Periodic) x= %f y=%f z=%f\n", meLocal, iPar, xp[0][iPar], yp[0][iPar],zp[0][iPar]);
		}
#endif

}

void ForwardComm() {

	MPI_Request request;
	MPI_Status status;
	Real* buf;
	int idim, meLocal, sizeRcv, sizeSnd;
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		idim = -1;
		MPI_Comm_rank(sb->comm, &meLocal);
		for (int iswap = 0; iswap < 2 * SPACEDIM; iswap++) {

			if (iswap % 2 == 0)
				idim ++;
			if (periodic[idim] == 0) //Non-periodic in direction idim
				continue;

			if ((sendNum[blockIndex][iswap] > 0) || (recvNum[blockIndex][iswap] > 0)) {

				if (recvNum[blockIndex][iswap] > 0) { //Data will be received in this swap
					if (procRecv[blockIndex][iswap] != meLocal) {
						sizeRcv = recvNum[blockIndex][iswap] * sizeForward;
						MPI_Irecv(buffRecv, sizeRcv, MPI_DOUBLE, procRecv[blockIndex][iswap], 0, sb->comm, &request);
					}
				}
				//packing data to be send
				if (sendNum[blockIndex][iswap] > 0) {//Data will be send in this swap
#ifdef OPS_2D
					packForward2D(blockIndex, sendNum[blockIndex][iswap], iswap, sendList[blockIndex][iswap], buffSend, periodic, pbcFlag[blockIndex]);
#endif

#ifdef OPS_3D
					packForward3D(blockIndex, sendNum[blockIndex][iswap], iswap, sendList[blockIndex][iswap], buffSend, periodic, pbcFlag[blockIndex]);
#endif
					if (procSend[blockIndex][iswap] != meLocal) {
						sizeSnd = sendNum[blockIndex][iswap] * sizeForward;
						MPI_Send(buffSend, sizeSnd, MPI_DOUBLE, procSend[blockIndex][iswap], 0, sb->comm);
					}
				}

				if (recvNum[blockIndex][iswap] > 0) {
					if (procRecv[blockIndex][iswap] == meLocal)
						buf = buffSend;
					else {
						MPI_Wait(&request, &status);
						buf = buffRecv;
					}
#ifdef OPS_2D
					unpackComm2D(blockIndex, recvNum[blockIndex][iswap], firstRecv[blockIndex][iswap], buf);
#endif

#ifdef OPS_3D
					unpackComm3D(blockIndex, recvNum[blockIndex][iswap], firstRecv[blockIndex][iswap], buf);
#endif
				}


			}

		}
	}
}

void ReverseComm() {

	int idim, meLocal, size;
	Real * buff;
	MPI_Status status;
	MPI_Request request;

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif


		idim = SPACEDIM;
		MPI_Comm_rank(sb->comm, &meLocal);
		for (int iswap = 2 * SPACEDIM - 1; iswap >= 0; iswap--) {
			if (iswap %2 != 0)
				idim--;

			if (periodic[idim] == 0)
				continue;

			//printf("Rank %d Swap %d: Receiving number %d\n",ops_get_proc(), iswap, reverseRecvNum[blockIndex][iswap]);
			if (reverseRecvNum[blockIndex][iswap] > 0) {
				if (procSend[blockIndex][iswap] != meLocal) {
					size = reverseRecvNum[blockIndex][iswap] * sizeReverse;
					MPI_Irecv(buffRecv,size,MPI_DOUBLE, procSend[blockIndex][iswap],0,sb->comm,&request);

				}
			}

			if (reverseSendNum[blockIndex][iswap] > 0) {
#ifdef OPS_2D
				packReverse2D(blockIndex, reverseSendNum[blockIndex][iswap], firstRecv[blockIndex][iswap],buffSend);
#endif

#ifdef OPS_3D
				packReverse3D(blockIndex, reverseSendNum[blockIndex][iswap], firstRecv[blockIndex][iswap],buffSend);
#endif
				if (procRecv[blockIndex][iswap] != meLocal) {
					size = reverseSendNum[blockIndex][iswap] * sizeReverse;
					MPI_Send(buffSend, size ,MPI_DOUBLE, procRecv[blockIndex][iswap],0,sb->comm);
				}
			}

			if (reverseRecvNum[blockIndex][iswap] > 0) {
				if (procSend[blockIndex][iswap] != meLocal) {
					MPI_Wait(&request, &status);
					buff = buffRecv;
				}
				else
					buff = buffSend;
#ifdef OPS_2D
				unPackReverse2D(blockIndex, reverseRecvNum[blockIndex][iswap], sendList[blockIndex][iswap], buff);
#endif

#ifdef OPS_3D
				unPackReverse3D(blockIndex, reverseRecvNum[blockIndex][iswap], sendList[blockIndex][iswap], buff);
#endif
			}

		}
	}


}

#ifdef OPS_2D
void packForward2D(int blockIndex, int nsend, int iswap, int* sendList, Real* buff, int* periodic, int* pbc) {

	Real dx, dy;
	int iPar;
	dx = static_cast<Real>(pbc[iswap * SPACEDIM + 0]) * xper[0];
	dy = static_cast<Real>(pbc[iswap * SPACEDIM + 1]) * xper[1];

	int m = 0;
	for (int ibuff = 0; ibuff < nsend; ibuff++) {
		iPar = sendList[ibuff];
		buff[m++] = xp[blockIndex][iPar] + dx;
		buff[m++] = yp[blockIndex][iPar] + dy;
		buff[m++] = Radius[blockIndex][iPar];
		buff[m++] = up[blockIndex][iPar];
		buff[m++] = vp[blockIndex][iPar];
		buff[m++] = omegaZ[blockIndex][iPar];

	}

}


void unpackComm2D(int blockIndex, int nrecv, int first, Real *buf) {

	int m = 0;
	int last = first + nrecv;

	for (int iPar = first; iPar < last; iPar++) {
		if (Nmax[blockIndex] == iPar)
			reallocate2dMemory(100, blockIndex);
		xp[blockIndex][iPar] = buf[m++];
		yp[blockIndex][iPar] = buf[m++];
		Radius[blockIndex][iPar] = buf[m++];
		up[blockIndex][iPar] = buf[m++];
		vp[blockIndex][iPar] = buf[m++];
		omegaZ[blockIndex][iPar] = buf[m++];
	}

}

void packReverse2D(int blockIndex, int nParticles, int nfirst, Real* buff) {

	int m = 0;
	int nlast = nfirst + nParticles;

	for (int iPar = nfirst; iPar < nlast; iPar++) {
		buff[m++] = FDrag[blockIndex][3*iPar];
		buff[m++] = FDrag[blockIndex][3 * iPar + 1];
		buff[m++] = FDrag[blockIndex][3 * iPar + 2];

		FDrag[blockIndex][3 * iPar] = 0.0;
		FDrag[blockIndex][3 * iPar + 1] = 0.0;
		FDrag[blockIndex][3 * iPar + 2] = 0.0;
	}


}

#endif


#ifdef OPS_3D
void packForward3D(int blockIndex, int nsend, int iswap, int* sendList,Real* buff,int* periodic,int * pbc) {

	Real dx, dy, dz;
	int iPar;
	dx = static_cast<Real>(pbc[iswap * SPACEDIM + 0]) * xper[0];
	dy = static_cast<Real>(pbc[iswap * SPACEDIM + 1]) * xper[1];
	dz = static_cast<Real>(pbc[iswap * SPACEDIM + 2]) * xper[2];

	int m = 0;
	for (int ibuff = 0; ibuff < nsend; ibuff++) {
		iPar = sendList[ibuff];
		buff[m++] = xp[blockIndex][iPar] + dx;
		buff[m++] = yp[blockIndex][iPar] + dy;
		buff[m++] = zp[blockIndex][iPar] + dz;
		buff[m++] = Radius[blockIndex][iPar];

		buff[m++] = up[blockIndex][iPar];
		buff[m++] = vp[blockIndex][iPar];
		buff[m++] = wp[blockIndex][iPar];

		buff[m++] = omegaX[blockIndex][iPar];
		buff[m++] = omegaY[blockIndex][iPar];
		buff[m++] = omegaZ[blockIndex][iPar];
	}


}


void unpackComm3D(int blockIndex, int nrecv, int first, Real *buf) {

	int m = 0;
	int last = first + nrecv;

	for (int iPar = first; iPar < last; iPar++) {
		if (Nmax[blockIndex] == iPar)
			ReAllocateMemory(100, blockIndex);
		xp[blockIndex][iPar] = buf[m++];
		yp[blockIndex][iPar] = buf[m++];
		zp[blockIndex][iPar] = buf[m++];
		Radius[blockIndex][iPar] = buf[m++];
		up[blockIndex][iPar] = buf[m++];
		vp[blockIndex][iPar] = buf[m++];
		wp[blockIndex][iPar] = buf[m++];

		omegaX[blockIndex][iPar] = buf[m++];
		omegaY[blockIndex][iPar] = buf[m++];
		omegaZ[blockIndex][iPar] = buf[m++];
	}

}

void packReverse3D(int blockIndex, int nParticles, int nfirst, Real* buff) {

	int m = 0;
	int nlast = nfirst + nParticles;

	for (int iPar = nfirst; iPar < nlast; iPar++) {
		buff[m++] = FDrag[blockIndex][6 * iPar];
		buff[m++] = FDrag[blockIndex][6 * iPar + 1];
		buff[m++] = FDrag[blockIndex][6 * iPar + 2];
		buff[m++] = FDrag[blockIndex][6 * iPar + 3];
		buff[m++] = FDrag[blockIndex][6 * iPar + 4];
		buff[m++] = FDrag[blockIndex][6 * iPar + 5];

		FDrag[blockIndex][6 * iPar] = 0.0;
		FDrag[blockIndex][6 * iPar + 1] = 0.0;
		FDrag[blockIndex][6 * iPar + 2] = 0.0;
		FDrag[blockIndex][6 * iPar + 3] = 0.0;
		FDrag[blockIndex][6 * iPar + 4] = 0.0;
		FDrag[blockIndex][6 * iPar + 5] = 0.0;
	}


}

void unPackReverse3D(int blockIndex, int nTot, int* sendList, Real* buff) {
	int m = 0;
	int jPar;


	m = 0;
	for (int iPar = 0; iPar < nTot; iPar++) {
		jPar = sendList[iPar];

		FDrag[blockIndex][6 * jPar] += buff[m++];

		FDrag[blockIndex][6 * jPar + 1] += buff[m++];

		FDrag[blockIndex][6 * jPar + 2] += buff[m++];

		FDrag[blockIndex][6 * jPar + 3] += buff[m++];

		FDrag[blockIndex][6 * jPar + 4] += buff[m++];

		FDrag[blockIndex][6 * jPar + 5] += buff[m++];
	}
}


void SetupVertexType(VertexType* vertexType) {
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
}


#endif
