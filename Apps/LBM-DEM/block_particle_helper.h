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
 * @brief   Wrapper functions for fluid-particle interaction and particle data
 * @author  C. Tsigginos
 */


#ifndef BLOCK_PARTICLES_HELPER_H_
#define BLOCK_PARTICLES_HELPER_H_

#include "block.h"
#include "scheme.h"
#include "block_particles.h"
#include "flowfield.h"
#include "type.h"
#include "ops_lib_core.h"
#include "fpi_functions.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include <string>
#include <vector>

void ParticleMapping(int flag);

int CheckParticleDistance();
void UpdateOldParticleLocation();
void ParticleEnvelopes();



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

void UpdateParticleMappingDragForceRestart(SizeType currentStep);

void DefineBlockParticles(int spacedim, Real cutoff, Real dx1,
		std::string particleType);

void DefineMappingModel(std::vector<ParticleMappingModel> mappingModels,
						std::vector<int> compoId, std::vector<int> copyFrom,
						ParticleShapeDiscriptor particleShape, SizeType timeStep = 0);

void  DefineFsiModel(std::vector<FluidParticleModel> fsiModel, std::vector<int> compoId,
			std::vector<Real> variable, SizeType timeStep = 0);


//PArticle assigning wrappers
void AssignParticlesToBlocksSpheres(int Nparticles,Real* xTmp, Real* yTmp,Real* zTmp,
		Real* radTmp, Real*  uTmp,Real* vTmp, Real* wTmp, Real* oxTmp,
		Real* oyTmp, Real* ozTmp);

void DefineBlockOwnership();
void DefineLocalBoxBound();

#endif /* APPS_LBM_DEM_BLOCK_PARTICLES_WRAPPER_H_ */
