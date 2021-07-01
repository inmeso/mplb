/*
 * dem_particles.cpp
 *
 *  Created on: Sep 18, 2020
 *      Author: jpd38567
 */
#include "dem_particles.h"
#include <stdlib.h>
#include <stdio.h>

int* StenList; //Cells containing a particles
int* Nparticles; //Number of particles per block
int* Nmax; // Maximum particles that can be inserted in a block without memory re-allocation
int* Nperiodic; //Periodic particles list

Real** xp; //Particle location in x
Real** yp; //Particle location in y
Real** zp; //Particle location iz z
Real** up; //Particle velocity in u
Real** vp; //Particle velocity in v
Real** wp; //Particle velocity in w
Real** omegaX; //Rotational velocity in x
Real** omegaY; //Rotational velocity in y
Real** omegaZ; //Rotational velocity in w
Real** Radius; //Particle radius

Real** FDrag;  //Drag Force
int** ParticleList; //Particle stencil

Real** xp_old; //Old particle location in x
Real** yp_old; //Old particle location in y
Real** zp_old; //Old particle location in z


void particleAllocation(int blockNo) {

	Nparticles = (int *) malloc(blockNo * sizeof(int));
	Nmax = (int *) malloc(blockNo * sizeof(int));

	Nperiodic = (int *) malloc(blockNo * sizeof(int));

	xp = (Real **) malloc(blockNo * sizeof(Real *));
	yp = (Real **) malloc(blockNo * sizeof(Real *));

	up = (Real **) malloc(blockNo * sizeof(Real *));
	vp = (Real **) malloc(blockNo * sizeof(Real *));

	omegaZ = (Real **) malloc(blockNo * sizeof(Real *));
	Radius = (Real **) malloc(blockNo * sizeof(Real *));


	zp = (Real **) malloc(blockNo * sizeof(Real *));
	wp = (Real **) malloc(blockNo * sizeof(Real *));

	omegaX = (Real **) malloc(blockNo * sizeof(Real *));
	omegaY = (Real **) malloc(blockNo * sizeof(Real *));



	for (int blockIndex = 0; blockIndex < blockNo; blockIndex++) {
		Nmax[blockIndex] = 10;
		Nparticles[blockIndex] = 0;
		Nperiodic[blockIndex] = 0;
		allocateMemory(xp[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(yp[blockIndex], Nmax[blockIndex], 1);

		allocateMemory(up[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(vp[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(Radius[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(omegaZ[blockIndex], Nmax[blockIndex], 1);

		allocateMemory(zp[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(wp[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(omegaX[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(omegaY[blockIndex], Nmax[blockIndex], 1);

	}

}

void AllocateFluidParticleInteractionVeriables() {

	int blockNo = BlockNum();
	ParticleList = (int **) malloc(blockNo * sizeof(int *));
	FDrag = (Real **) malloc(blockNo * sizeof (Real *));
	xp_old = (Real **) malloc(blockNo * sizeof(Real *));
	yp_old = (Real **) malloc(blockNo * sizeof(Real *));


	zp_old = (Real **) malloc(blockNo * sizeof(Real *));


	for (int blockIndex = 0; blockIndex < blockNo; blockIndex++) {

		allocateMemory(xp_old[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(yp_old[blockIndex], Nmax[blockIndex], 1);
		allocateMemory(zp_old[blockIndex], Nmax[blockIndex], 1);

		allocateMemory(ParticleList[blockIndex], Nmax[blockIndex], 6);
		allocateMemory(FDrag[blockIndex], Nmax[blockIndex], 6);

	}

}


void ReallocateMemoryFSI(int blockIndex) {
	Reallocate(xp_old[blockIndex], Nmax[blockIndex], 1);
	Reallocate(yp_old[blockIndex], Nmax[blockIndex], 1);


	Reallocate(zp_old[blockIndex], Nmax[blockIndex], 1);
	Reallocate(ParticleList[blockIndex], Nmax[blockIndex], 6);
	Reallocate(FDrag[blockIndex], Nmax[blockIndex], 6);

}

void ReAllocateMemory(int Np, int blockIndex) {
	Nmax[blockIndex] += Np;

	Reallocate(xp[blockIndex], Nmax[blockIndex], 1);
	Reallocate(up[blockIndex], Nmax[blockIndex], 1);

	Reallocate(yp[blockIndex], Nmax[blockIndex], 1);
	Reallocate(vp[blockIndex], Nmax[blockIndex], 1);

	Reallocate(omegaZ[blockIndex], Nmax[blockIndex], 1);


	Reallocate(zp[blockIndex], Nmax[blockIndex], 1);
	Reallocate(wp[blockIndex], Nmax[blockIndex], 1);

	Reallocate(omegaX[blockIndex], Nmax[blockIndex], 1);
	Reallocate(omegaY[blockIndex], Nmax[blockIndex], 1);


	Reallocate(Radius[blockIndex], Nmax[blockIndex], 1);

	ReallocateMemoryFSI(blockIndex);

}




