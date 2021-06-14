/*
 * block_particles_wrapper.h
 *
 *  Created on: May 14, 2021
 *      Author: jpd38567
 */

#ifndef BLOCK_PARTICLES_HELPER_H_
#define BLOCK_PARTICLES_HELPER_H_

#include "block.h"
#include "scheme.h"
#include "block_particles.h"
#include "flowfield.h"
#include "type.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "fpi.h"
#include <string>
#include <vector>

void ParticleMapping(int flag);

int CheckParticleDistance();
void UpdateOldParticleLocation();
void ParticleEnvelopes();
void MappingParticlesToLBMGrid();
void UpdateParticlesToLBMGrid();


//Wrapper for FSI object functions
void PostVelocityFSIFunctions();
void PreCollisionFSIFunctions();
void PostStreamingFSIFunctions();
void CalculateParticleMomentum();

void FluidParticleCollisions();
void PreDefineCollision3DFSI(int* velID, int* loop, Real tauRef, CollisionType collisionType);

void CalculateDragForce(Real dt, SizeType currentStep);
void InitializeDragForce();
void InitializeFSILists();
void UpdateFPIVelocities3D();
void WriteFPIDataToHdf5(SizeType currentStep);

void DefineBlockParticles(int spacedim, Real cutoff, Real dx1,
		std::string particleType);

void DefineInteractionModel(std::vector<FSIType> FluidParticleInteractionType,
		std::vector<int> porosModel, std::vector<int> fsiCompoId, Real* forceUser,
		double gamma = 0, SizeType timeStep = 0);


//PArticle assigning wrappers
void AssignParticlesToBlocksSpheres(int Nparticles,Real* xTmp, Real* yTmp,Real* zTmp,
		Real* radTmp, Real*  uTmp,Real* vTmp, Real* wTmp, Real* oxTmp,
		Real* oyTmp, Real* ozTmp);

void DefineBlockOwnership();
void DefineLocalBoxBound();

#endif /* APPS_LBM_DEM_BLOCK_PARTICLES_WRAPPER_H_ */
