/*
 * poros_grid.cpp
 *
 *  Created on: May 29, 2021
 *      Author: jpd38567
 */

#include "poros_grid.h"
#include "flowfield.h"
#include "ops_seq_v2.h"
#include "block_particles.h"
#include <map>
#include <cassert>

#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif

#include "poros_grid.inc"
PorosGrid::PorosGrid(int particleshape, int spacedim) :
	ParticleToGridBase( particleshape,  spacedim) {

	//Assign variables
	RealField poros{"SolFracGrid"};
	mappingRealVariableList.push_back(poros);
	RealField vParts{"velPartGrid"};
	mappingRealVariableList.push_back(vParts);
	IntField idParts{"idParticlesGrid"};
	mappingIntVariableList.push_back(idParts);
	RealField xAverage{"xAvgGrid"};
	mappingRealVariableList.push_back(xAverage);

}

void PorosGrid::DefineVariables(int noelem, SizeType timeStep) {


	noElem = noelem;
	int size = noElem * spaceDim;
	mappingRealVariableList.at(0).SetDataDim(noElem);
	mappingRealVariableList.at(1).SetDataDim(size);
	mappingRealVariableList.at(2).SetDataDim(size);

	for (auto& mappingVariable : mappingRealVariableList) {
		if (timeStep == 0) {
			mappingVariable.CreateFieldFromScratch(g_Block());
		}
		else {
			mappingVariable.CreateFieldFromFile(CaseName(), g_Block(), timeStep);

		}

	}

	mappingIntVariableList.at(0).SetDataDim(noElem);
	if (timeStep == 0)
		mappingIntVariableList.at(0).CreateFieldFromScratch(g_Block());
	else
		mappingIntVariableList.at(0).CreateFieldFromFile(CaseName(), g_Block(), timeStep);

}

void PorosGrid::MappingFunction(bool flag) {

	if (flag)
		ParticleProjection();
	else
		UpdateProjection();

}

void PorosGrid::ParticleProjection() {

	int size = noElem * spaceDim;
	int idParticle;
	int blockIndex;
	const Real Dx {GetDx()};
	Real xPos[spaceDim], Radius, uParticle[spaceDim], omParticle[spaceDim];
	int stenList[2 * spaceDim];
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		if ( !ParticleCurrentBlock.owned ) continue;
		blockIndex = ParticleCurrentBlock.GetBlock().ID();
		int nlocal =  ParticleCurrentBlock.NParticles +
				ParticleCurrentBlock.NPeriodic;
		for (int iPart = 0; iPart < nlocal; iPart++) {
			xPos[0] = ParticleCurrentBlock.particleList.at(iPart).xParticle[0];
			xPos[1] = ParticleCurrentBlock.particleList.at(iPart).xParticle[1];
			xPos[2] = ParticleCurrentBlock.particleList.at(iPart).xParticle[2];
			Radius = ParticleCurrentBlock.particleList.at(iPart).particleShape->Rparticle;

			uParticle[0] = ParticleCurrentBlock.particleList.at(iPart).uParticle[0];
			uParticle[1] = ParticleCurrentBlock.particleList.at(iPart).uParticle[1];
			uParticle[2] = ParticleCurrentBlock.particleList.at(iPart).uParticle[2];

			omParticle[0] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[0];
			omParticle[1] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[1];
			omParticle[2] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[2];
			idParticle = iPart;
			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
				stenList[iDir] = ParticleCurrentBlock.particleList.at(iPart).stenList[iDir];

			switch (particleDiscriptor) {
				case (ParticleSpherical):
					 ops_par_loop(KerSolidFracGridSphere, "KerSolidFracGridSphere",
							 ParticleCurrentBlock.GetBlock().Get(), spaceDim,  stenList,
							 ops_arg_dat(mappingIntVariableList.at(0).at(blockIndex), noElem,
									 	 LOCALSTENCIL,"int", OPS_RW),
							 ops_arg_dat(mappingRealVariableList.at(0).at(blockIndex), noElem,
								 LOCALSTENCIL,"double", OPS_RW),
							 ops_arg_dat(mappingRealVariableList.at(1).at(blockIndex), size,
								 LOCALSTENCIL,"double", OPS_RW),
							 ops_arg_dat(mappingRealVariableList.at(2).at(blockIndex), size,
								 LOCALSTENCIL,"double", OPS_RW),
							 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
				    		 	 LOCALSTENCIL,"double", OPS_READ),
						     ops_arg_gbl(xPos, spaceDim, "double", OPS_READ),
						     ops_arg_gbl(&Radius, 1, "double", OPS_READ),
						     ops_arg_gbl(&idParticle, 1, "double", OPS_READ),
						     ops_arg_gbl(uParticle, spaceDim, "double", OPS_READ),
						     ops_arg_gbl(omParticle, spaceDim, "double", OPS_READ),
						     ops_arg_gbl(&Dx, 1, "double", OPS_READ),
						     ops_arg_gbl(&Npoints, 1, "double", OPS_READ),
							 ops_arg_gbl(&noElem, 1, "double", OPS_READ),
						     ops_arg_gbl(&spaceDim, 1, "double", OPS_READ));
					break;
				default:
					break;
			}

		}
	}



}

void PorosGrid::UpdateProjection() {

	int size = noElem * spaceDim;
	Real velP[spaceDim], omP[spaceDim];
	Real xPar[spaceDim], radius;
	int stenList[2 * spaceDim];
	int nlocal;
	int idParticle;
	Real dx = GetDx();
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		if ( !ParticleCurrentBlock.owned ) continue;
		int blockIndex = ParticleCurrentBlock.GetBlock().ID();
		int nlocal =  ParticleCurrentBlock.NParticles +
				ParticleCurrentBlock.NPeriodic;
		for (int iPart = 0; iPart < nlocal; iPart++) {
			xPar[0] = ParticleCurrentBlock.particleList.at(iPart).xParticle[0];
			xPar[1] = ParticleCurrentBlock.particleList.at(iPart).xParticle[1];
			xPar[2] = ParticleCurrentBlock.particleList.at(iPart).xParticle[2];
			radius = ParticleCurrentBlock.particleList.at(iPart).particleShape->Rparticle;

			velP[0] = ParticleCurrentBlock.particleList.at(iPart).uParticle[0];
			velP[1] = ParticleCurrentBlock.particleList.at(iPart).uParticle[1];
			velP[2] = ParticleCurrentBlock.particleList.at(iPart).uParticle[2];

			omP[0] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[0];
			omP[1] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[1];
			omP[2] = ParticleCurrentBlock.particleList.at(iPart).omegaParticle[2];
			idParticle = iPart;
			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
				stenList[iDir] = ParticleCurrentBlock.particleList.at(iPart).stenList[iDir];

			switch (particleDiscriptor) {
					case (ParticleSpherical):
					ops_par_loop(KerSolidVelocityUpdate,"KerSolidVelocityUpdate",
						 ParticleCurrentBlock.GetBlock().Get(), spaceDim, stenList,
						 ops_arg_dat(mappingRealVariableList.at(1).at(blockIndex), size,
								     LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(mappingIntVariableList.at(0).at(blockIndex), noElem,
							 	 	 LOCALSTENCIL,"int", OPS_READ),
						 ops_arg_dat(mappingRealVariableList.at(2).at(blockIndex), size,
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(&idParticle, 1, "int", OPS_READ),
						 ops_arg_gbl(xPar, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&radius, 1, "double", OPS_READ),
						 ops_arg_gbl(velP, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(omP, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&dx, 1, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
						 ops_arg_gbl(&noElem, 1, "int", OPS_READ));

					break;
					default:
						ops_printf("Only spherical particles are implemented\n");
			}
		}
	}
}

void PorosGrid::InitializeVariables() {
	int size = noElem * spaceDim;
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		std::vector<int> iterRng;
		iterRng.assign(ParticleCurrentBlock.GetBlock().WholeRange().begin(),
					   ParticleCurrentBlock.GetBlock().WholeRange().end());

		const int blockIndex{ParticleCurrentBlock.GetBlock().ID()};
		ops_par_loop(KerInitializePorousGrid,"KerInitializePorousGrid",
				ParticleCurrentBlock.GetBlock().Get(),spaceDim, iterRng.data(),
				ops_arg_dat(mappingRealVariableList.at(0).at(blockIndex), size,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingRealVariableList.at(1).at(blockIndex), noElem,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingRealVariableList.at(2).at(blockIndex), size,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingIntVariableList.at(0).at(blockIndex), size,
							LOCALSTENCIL, "int", OPS_WRITE),
				ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
				ops_arg_gbl(&noElem, 1, "int", OPS_READ));




	}

}



