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



#include "block_particle_helper.h"
#include <map>
//#include "model_kernel.inc"
#include "dem_data.h"
#include "block_particles.h"
#include "fpi_functions.h"
#include "particle_mapping.h"
#include "mapping_particles.h"
#include <memory>

//#include "particle.h"
//#include <vector>

void DefineBlockParticles( int spacedim, Real cutoff, Real dx1, std::string particleType) {


	ParticleShapeDiscriptor particleTypeNum;

	//Compare strings to find particletype
	std::string s1{"spherical"};
	std::string s2{"quadratic"};
	std::string s3{"mesh"};
	if (particleType.compare(s1) == 0)
		particleTypeNum = spherical;
	else if (particleType.compare(s2) == 0)
		particleTypeNum = quadratic;
	else if (particleType.compare(s3) == 0)
		particleTypeNum = mesh;
	else {
		ops_printf("MPLB: This type of particle is not supported\n");
		ops_printf("MPLB: Reversing to spherical particles\n");
		particleTypeNum = spherical;
	}

	for (const auto& idBlock : g_Block()) {
		const Block& myBlock{idBlock.second};

		BlockParticles ParticleNew(spacedim, dx1, cutoff, myBlock, particleTypeNum,0 ,0);
		BlockParticleList.emplace(idBlock.first, ParticleNew);

	}


}

void DefineMappingModel(std::vector<ParticleMappingModel> mappingModels,
						std::vector<int> compoId, std::vector<int> copyFrom,
						ParticleShapeDiscriptor particleShape, SizeType timeStep) {


	DefineParticleMappingModels(mappingModels, compoId,  copyFrom,  particleShape,
								timeStep);
}

void  DefineFsiModel(std::vector<FluidParticleModel> fsiModel, std::vector<int> compoId,
			std::vector<Real> variable, SizeType timeStep) {

	DefineInteractionModels(fsiModel, compoId, variable, timeStep);
}
void UpdateParticleMappingDragForceRestart(SizeType currentStep) {

	Real dt = TimeStep();
	InitializeDragForce();

	ParticleEnvelopes();

	CalculateParticleMomentum();


	CalculateDragForce(dt, currentStep);

}


void ParticleMapping(int flag) {

	int checkDistance;

	checkDistance = CheckParticleDistance();


	if (flag == 1)
		checkDistance = 1;


	if (checkDistance == 1) {
		UpdateOldParticleLocation();

		ParticleEnvelopes();

		InitializeFSILists(); //FSI schemes

		MappingParticlesToLBMGrid(); //FSI Wrappers

	}
	else
		UpdateParticlesToLBMGrid();
}


int CheckParticleDistance() {
	int  distanceGlobal, distanceTmp;
	distanceGlobal = 0;
	for ( auto& blockParticle : BlockParticleList) {
		if (blockParticle.second.OwnedStatus()) {
			distanceTmp = blockParticle.second.CheckDistanceBlock();
			if (distanceTmp == 1) {
				distanceGlobal = 1;
				break;
			}
		}
	}

#ifdef OPS_MPI
	distanceTmp = distanceGlobal;
	if (ops_num_procs()>1)
		MPI_Allreduce(&distanceTmp, &distanceGlobal, 1, MPI_INT, MPI_MAX, OPS_MPI_GLOBAL);
#endif

#ifdef CPU
#if DebugLevel >= 2
	printf("Rank %d: distanceGlobal = %d\n",ops_get_proc(), distanceGlobal);
#endif
#endif

	return distanceGlobal;

}

void UpdateOldParticleLocation() {

	for (auto & blockParticle : BlockParticleList) {
		if (blockParticle.second.OwnedStatus())
			blockParticle.second.UpdateOldParticlePosition();
	}

}

void ParticleEnvelopes() {

	for (auto & blockParticle : BlockParticleList) {
		if (blockParticle.second.OwnedStatus())
			blockParticle.second.FindStencil();
	}

	/*	for (auto &blockParticle : BlockParticleList) {
			int iD = blockParticle.second.GetBlock().ID();
				if (blockParticle.second.OwnedStatus()) {
					int Np = blockParticle.second.NParticles;
					for (int iPar = 0; iPar < Np ; iPar++) {
						Particle& particle = blockParticle.second.particleList.at(iPar);
						printf("Rank %d at block %d: Stencil of Particle %d ([%f %f %f])[%d %d] x [%d %d] x [%d %d] cells\n",
								ops_get_proc(), iD, iPar, particle.xParticle[0], particle.xParticle[1], particle.xParticle[2],
								particle.stenList[0], particle.stenList[1], particle.stenList[2],
								particle.stenList[3], particle.stenList[4], particle.stenList[5]);
					}
				}
			}*/

}



//Modified
void PostVelocityFSIFunctions() {

	for (auto& fsi : fpiModels) {
		if (fsi.second->OwnedPostVelocity())
			FSIVelocityFunctions(fsi.second);
	}
}

//Modified
void PreCollisionFSIFunctions() {

	for (auto &fsi : fpiModels) {
		if (fsi.second->OwnedPreCollision())
			FsiPreCollisionFunction(fsi.second);
	}

}

//modified
void PostStreamingFSIFunctions() {

	for (auto &fsi : fpiModels) {
		if (fsi.second->OwnedPostStreaming())
			FsiPostCollisionFunction(fsi.second);
	}

}
//modified
void CalculateParticleMomentum() {

	for (auto &fsi : fpiModels) {
		if (fsi.second->OwnedDragForce())
			FsiCalculateDragForce(fsi.second);
	}
}

//modified
void InitializeFSILists() {

	for (auto &fsi : fpiModels) {
		if (fsi.second->OwnedInitialize());
		FsiInitializeFunction(fsi.second);
	}

	InitializeMappingLists();
}
//modified.
void FluidParticleCollisions() {

	int loop[2];

#if OPS_3D
	int velId[3];
#endif

#if OPS_2D
	int velId[2];
#endif

	CollisionType collisionModel;
	Real tauRef;
	int componentId, rhoId, Tid;
	for (auto &fsi : fpiModels) {

		if (fsi.second->OwnedCollisionModel()) {
			FsiCollisionFunction(fsi.second);
		}
		else {
			//TODO need to obtain the id of the
			ObtainData(fsi.second, velId,loop, tauRef, collisionModel, componentId, rhoId,
					Tid);

			PreDefinedCollision3D(velId, loop, tauRef, collisionModel, componentId, rhoId, Tid);
		}

	}
}

void CalculateDragForce(Real dt, SizeType currentStep) {

	for (auto& blockParticle : BlockParticleList) {
		blockParticle.second.EvaluateDragForce(dt);
	}
}

void InitializeDragForce() {

	for (auto& blockParticle : BlockParticleList) {
		if (blockParticle.second.OwnedStatus())
			blockParticle.second.InitializeDragForce();

	}


}

void UpdateFPIVelocities3D() {

	UpdateMacroVars3D();
	PostVelocityFSIFunctions();


}


void WriteFPIDataToHdf5(SizeType currentStep) {

	for (auto &fsi : fpiModels) {
		fsi.second->WriteToHdf5(CaseName(), currentStep);
	}

	WriteParticleMappingToHDF5(currentStep);

}

void AssignParticlesToBlocksSpheres(int Nparticles,Real* xTmp, Real* yTmp,Real* zTmp,
		Real* radTmp, Real*  uTmp,Real* vTmp, Real* wTmp, Real* oxTmp,
		Real* oyTmp, Real* ozTmp) {

	int spaceDim;
#ifdef OPS_3D
	 spaceDim = 3;
#endif

#ifdef OPS_2D
	 spaceDim = 2;
#endif
	 int idx;
	 std::vector<Real> shape(1.0);
	 std::vector<Real> extra;
	 shape.reserve(1);
	 int ipx = 0;
	 Real xParticle[spaceDim], uParticle[spaceDim], omParticle[spaceDim];
	 for (int iPar = 0; iPar < Nparticles; iPar++) {
			 ipx += 1;
		 	 xParticle[0] = xTmp[iPar];
			 xParticle[1] = yTmp[iPar];
			 xParticle[2] = zTmp[iPar];
			 uParticle[0] = uTmp[iPar];
			 uParticle[1] = vTmp[iPar];
			 uParticle[2] = wTmp[iPar];
			 omParticle[0] = oxTmp[iPar];
			 omParticle[1] = oyTmp[iPar];
			 omParticle[2] = ozTmp[iPar];
			 shape.at(0) = radTmp[iPar];
			 int ix = 0;
			 for (auto &blockParticle : BlockParticleList) {
				 ix = ix + 1;
				 idx = blockParticle.second.InsertParticle(xParticle, radTmp[iPar], shape,
							uParticle, omParticle, extra);
				 printf("Rank %d: Idx = %d for particle %d\n", ops_get_proc(), idx, iPar );
			 }
	 }


}

void DefineBlockOwnership() {

	for (auto& blockParticle : BlockParticleList)
		blockParticle.second.GetOwnership();

}

/*void DefineLocalBoxBound() {

	for (auto& blockParticle : BlockParticleList) {
		if (blockParticle.second.owned)
			blockParticle.second.FindBoxLocalBound();
	}
}
*/
