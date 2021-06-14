/*
 * Block_particles_wrapper.cpp
 *
 *  Created on: May 14, 2021
 *      Author: jpd38567
 */


#include "block_particle_helper.h"
#include <map>
//#include "model_kernel.inc"
#include "dem_data.h"
#include "block_particles.h"
std::map<int, FsiBase*> fluidPartInteractionModels;

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

void DefineInteractionModel(std::vector<FSIType> FluidParticleInteractionType,
		 std::vector<int> fsiCompoId, Real* forceUser,
		double gamma, SizeType timeStep) {

	int noComp = FluidParticleInteractionType.size();
	int noFSICompo = fsiCompoId.size();
	FsiBase* model;
	FSIType fsiModel;
	int idCompo;
	std::vector<SolFracType> porosModel;
	int nelem;

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

	ParticleShapeDiscriptor particleType = (ParticleShapeDiscriptor) BlockParticleList
			.at(0).GetParticleShape();
	//find particle mapping schemes
	int iter  = 0;
	for (int iComp = 0; iComp < noComp; iComp++) {
		fsiModel =  FluidParticleInteractionType.at(iComp);
		if (fsiModel == Model_Prati || fsiModel == Model_PSM) {
			if (iter == 0) {
				if (particleType == spherical)
					porosModel.push_back(Mode_Spherical);
				else
					porosModel.push_back(Mode_Grid);
				iter = 1;
			}
			else
				porosModel.push_back(Mode_Copy);
		}
		else if (fsiModel == Model_None)
			porosModel.push_back(Mode_None);
		else {
			ops_printf("This type of model is not supported\n");
			exit(EXIT_FAILURE);
		}


	}

	if (particleType == spherical)
		nelem = 2;
	else
		nelem = 8;

	for (int iComp = 0; iComp < noComp; iComp++) {
		fsiModel = (FSIType) FluidParticleInteractionType.at(iComp);



		switch (fsiModel) {
			case Model_None:
				model = new FsiBase( components.at(iComp), SpaceDim(), forceUser, false,  porosModel.at(iComp), gamma);
				break;
			case Model_PSM:
				model = new Psm(components.at(iComp), SpaceDim(),forceUser, true , porosModel.at(iComp), gamma, nelem, particleType);
				break;
			case Model_Prati:
				model = new Prati(components.at(iComp), SpaceDim(), forceUser, true, porosModel.at(iComp), gamma, nelem, particleType);
				break;
			default:
				ops_printf("The chosen model is not supported\n");
				exit(EXIT_FAILURE);
		}
		idCompo = components.at(iComp).id;
		fluidPartInteractionModels.emplace(idCompo, model); //TODO DEfine it

	}

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
		distanceTmp = blockParticle.second.CheckDistanceBlock();
		if (distanceTmp == 1) {
			distanceGlobal = 1;
			break;
		}

	}

#ifdef OPS_MPI
	distanceTmp = distanceGlobal;
	MPI_Allreduce(&distanceTmp, &distanceGlobal, 1, MPI_INT, MPI_MAX, OPS_MPI_GLOBAL);
#endif
	return distanceGlobal;
}

void UpdateOldParticleLocation() {

	for (auto & blockParticle : BlockParticleList) {
		blockParticle.second.UpdateOldParticlePosition();
	}
}

void ParticleEnvelopes() {

	for (auto & blockParticle : BlockParticleList) {
		blockParticle.second.FindStencil();
	}
}

void MappingParticlesToLBMGrid() {

	for (auto& fsi : fluidPartInteractionModels) {
		fsi.second->MappingFunction(true);
	}



}

void UpdateParticlesToLBMGrid() {
	for (auto& fsi : fluidPartInteractionModels) {
		fsi.second->MappingFunction(false);
	}
}

void PostVelocityFSIFunctions() {

	for (auto &fsi : fluidPartInteractionModels) {
		fsi.second->PostVelocityCalculation();
	}
}

void PreCollisionFSIFunctions() {

	for (auto &fsi : fluidPartInteractionModels) {
		fsi.second->PreCollision();
	}

}

void PostStreamingFSIFunctions() {

	for (auto &fsi : fluidPartInteractionModels) {
		fsi.second->PostStreaming();
	}

}

void CalculateParticleMomentum() {

	for (auto &fsi :  fluidPartInteractionModels) {
		fsi.second->CalculateDragForce();
	}
}

void InitializeFSILists() {

	for (auto &fsi : fluidPartInteractionModels) {
		fsi.second->InitializeVariables();
	}
}

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
	for (auto &fsi : fluidPartInteractionModels) {
		if (fsi.second->collisionOwned) {
			fsi.second->ModelCollision();
		}
		else {
			//TODO need to obtain the id of the
			fsi.second->ObtainID(velId, loop, tauRef, collisionModel,componentId, rhoId, Tid);
			PreDefineCollision3DComponent(velId, loop, tauRef, collisionModel, componentId, rhoId, Tid);
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
		blockParticle.second.InitializeDragForce();

	}
}

void UpdateFPIVelocities3D() {

	UpdateMacroVars3D();
	PostVelocityFSIFunctions();


}

void WriteFPIDataToHdf5(SizeType currentStep) {

	for (auto &fsi : fluidPartInteractionModels) {
		fsi.second->WriteToHdf5(CaseName(), currentStep);
	}

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
	 std::vector<Real> shape;
	 shape.reserve(1);
	 Real xParticle[spaceDim], uParticle[spaceDim], omParticle[spaceDim];
	 for (int iPar = 0; iPar < Nparticles; iPar++) {
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
			 for (auto &blockParticle : BlockParticleList) {
				 idx = blockParticle.second.InsertParticle(xParticle, radTmp[iPar], shape,
							uParticle, omParticle);
			 }


	 }
}

void DefineBlockOwnership() {

	for (auto& blockParticle : BlockParticleList)
		blockParticle.second.GetOwnership();

}

void DefineLocalBoxBound() {

	for (auto& blockParticle : BlockParticleList) {
		if (blockParticle.second.owned)
			blockParticle.second.FindBoxLocalBound();
	}
}
