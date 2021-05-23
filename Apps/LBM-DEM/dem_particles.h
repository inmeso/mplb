/*
 * DEM_particles.h
 *
 *  Created on: Sep 18, 2020
 *      Author: jpd38567
 */

#ifndef DEM_PARTICLES_H_
#define DEM_PARTICLES_H_

#include "box_handling.h"
#include "memory_handle.h"


extern int* StenList; //Cells containing a particles

extern int* Nparticles; //Number of particles per block
extern int* Nmax; // Maximum particles that can be inserted in a block without memory re-allocation
extern int* Nperiodic; //Periodic particles list

extern Real** xp; //Particle location in x
extern Real** yp; //Particle location in y
extern Real** zp; //Particle location iz z
extern Real** up; //Particle velocity in u
extern Real** vp; //Particle velocity in v
extern Real** wp; //Particle velocity in w
extern Real** omegaX; //Rotational velocity in x
extern Real** omegaY; //Rotational velocity in y
extern Real** omegaZ; //Rotational velocity in w
extern Real** Radius; //Particle radius

extern Real** FDrag;  //Drag Force
extern int** ParticleList; //Not remembering

extern Real** xp_old; //Old particle location in x
extern Real** yp_old; //Old particle location in y
extern Real** zp_old; //Old particle location in z

void particleAllocation(int blockNo);
void AllocateFluidParticleInteractionVeriables();

void ReallocateMemoryFSI(int blockIndex);
void ReAllocateMemory(int Np, int blockIndex);

#endif /* DEM_PARTICLES_H_ */
