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

/*! @brief  Definition of mapping spherical particles with solid fraction
 * @author C. Tsigginos
 * @details Functions for evaluating solid fraction for spherical particles
 *
 */


#include "poros_spherical.h"
#include "flowfield.h"
#include "ops_seq_v2.h"
#include "block_particles.h"
#include <map>
#include <cassert>
#include "type.h"

#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif


#include "poros_spherical.inc"
PorosSpherical::PorosSpherical(int particleshape, int spacedim) :
	ParticleToGridBase( particleshape,  spacedim) {

	//Verifications for particle shape
	if (particleDiscriptor != ParticleSpherical) {
		ops_printf("ERROR: Class PorosSpherical supports only spherical"
				   "particles\n");
		assert(particleDiscriptor == ParticleSpherical);

	}

	//Assign variables
	RealField poros{"SolFracSph"};
	mappingRealVariableList.push_back(poros);
	RealField vParts{"velPartSph"};
	mappingRealVariableList.push_back(vParts);
	IntField idParts{"idParticlesSph"};
	mappingIntVariableList.push_back(idParts);
	RealField xAverage{"xAvgSph"};
	mappingRealVariableList.push_back(xAverage);

}

void PorosSpherical::DefineVariables(int noelem, SizeType timeStep) {


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
	for (auto& mappingVariable : mappingIntVariableList) {
		if (timeStep == 0) {
			mappingVariable.CreateFieldFromScratch(g_Block());
		}
		else {
			mappingVariable.CreateFieldFromFile(CaseName(), g_Block(), timeStep);

		}
	}

}
void PorosSpherical::MappingFunction(bool flag) {

	if (flag)
		ParticleProjection();
	else
		UpdateProjection();


}

void PorosSpherical::ParticleProjection() {

	ops_printf("I entered on particle projection\n");

	int size = noElem * spaceDim;
	int idParticle;
	int blockIndex;
	const Real Dx {GetDx()};
	Real xPos[spaceDim], Radius, uParticle[spaceDim], omParticle[spaceDim];
	int iterate[2 * spaceDim];
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		if ( !ParticleCurrentBlock.OwnedStatus() ) continue;
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
				iterate[iDir] = ParticleCurrentBlock.particleList.at(iPart).stenList[iDir];

			ops_par_loop(KerSolidFracSphere, "KerSolidFracSphere", ParticleCurrentBlock.GetBlock().Get(),
						 spaceDim,  iterate,
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
						 ops_arg_gbl(&noElem, 1, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "double", OPS_READ));

		}
	}



}

void PorosSpherical::UpdateProjection() {

	int size = noElem * spaceDim;
	Real velP[spaceDim], omP[spaceDim];
	Real xPar[spaceDim], radius;
	int stenList[2 * spaceDim];
	int nlocal;
	int idParticle;
	Real dx = GetDx();
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		if ( !ParticleCurrentBlock.OwnedStatus() ) continue;
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

			ops_par_loop(KerSolidVelocitySphere,"KerSolidVelocitySphere",
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
		}
	}

}

void PorosSpherical::InitializeVariables() {
	int size = noElem * spaceDim;
	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		std::vector<int> iterRng;
		iterRng.assign(ParticleCurrentBlock.GetBlock().WholeRange().begin(),
					   ParticleCurrentBlock.GetBlock().WholeRange().end());

		const int blockIndex{ParticleCurrentBlock.GetBlock().ID()};
		ops_par_loop(KerInitializePorousSpherical,"KerInitializePorousSpherical",
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

void PorosSpherical::PrintMappingVariables() {

	int size = noElem * spaceDim;

	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		std::vector<int> iterRng;
		iterRng.assign(ParticleCurrentBlock.GetBlock().WholeRange().begin(),
					   ParticleCurrentBlock.GetBlock().WholeRange().end());

		const int blockIndex{ParticleCurrentBlock.GetBlock().ID()};
		ops_par_loop(KerPrintPorousData,"KerPrintPorousData",
				ParticleCurrentBlock.GetBlock().Get(),spaceDim, iterRng.data(),
				ops_arg_dat(mappingRealVariableList.at(0).at(blockIndex), size,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingRealVariableList.at(1).at(blockIndex), noElem,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingRealVariableList.at(2).at(blockIndex), size,
							LOCALSTENCIL, "double", OPS_WRITE),
				ops_arg_dat(mappingIntVariableList.at(0).at(blockIndex), size,
							LOCALSTENCIL, "int", OPS_WRITE),
				ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
				    		 	 LOCALSTENCIL,"double", OPS_READ),
				ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
				ops_arg_gbl(&noElem, 1, "int", OPS_READ));
	}
}

