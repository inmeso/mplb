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

/*! @brief   Implement member functions for PRATI scheme
 * @author C. Tsigginos
 * @details Implementing functions for Partially saturate method.
*/

#include "psm.h"
#include "flowfield.h"
#include "flowfield_host_device.h"
#include "scheme.h"
#include <vector>
#include <map>
#include "ops_seq_v2.h"
#include "psm.inc"

Psm::Psm(Component componentUser, int spacedim, Real* forceUser, bool owned,
		SolFracType porosModel, Real gammaUser, int nelem, int ParticleType) :
		FsiBase(componentUser, spacedim, forceUser, owned, porosModel),
		Fd{"FdPSM"} {

	noElem = nelem;

	if ((porosModel!= Mode_Spherical) && (porosModel != Mode_Grid)) {
		ops_printf("ERROR: Implemented porosity model  %d not consistent with PSM scheme\n", porosModel);
		exit(EXIT_FAILURE);
	}

	if (porosModel == Mode_Spherical) {
		ops_printf("Porosity model for spherical particles for component %d\n",
					componentUser.id);
		poros = new PorosSpherical(ParticleType, spaceDim);;
	}
	else if (porosModel == Mode_Grid) {
		ops_printf("Grid model for porosity calculation for component %d\n",
						componentUser.id);
		poros = new PorosGrid(ParticleType, spaceDim);

	}
}

void Psm::DefineVariables(SizeType timestep) {

	int size = noElem * spaceDim;

	Fd.SetDataDim(size);
	if (timestep == 0) {
		Fd.CreateFieldFromScratch(g_Block());
		poros->DefineVariables(noElem);
	}
	else {
		Fd.CreateFieldFromFile(CaseName(), g_Block(), timestep);
		poros->DefineVariables(noElem, timestep);
	}

}

void Psm::ModelCollision() {

	int size = noElem * spaceDim;



	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());

		Real tauRef = compo.tauRef;
		const Real* pdt {pTimeStep()};
		const int blockIndex{block.ID()};
		ops_par_loop(KerCollisionPSM3D,"KerCollisionPSM3D", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_dat(Fd[blockIndex], size, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
							 	 "double", OPS_READ),
					 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
                                 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id)
								 .at(blockIndex), 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex),
								 1, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(poros->GetIntFieldVariable(0).at(blockIndex), noElem,
							     LOCALSTENCIL, "int", OPS_READ),
					 ops_arg_dat(poros->GetRealFieldVariable(0).at(blockIndex), noElem,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(poros->GetRealFieldVariable(1).at(blockIndex), size,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_gbl(pdt, 1, "double", OPS_READ),
					 ops_arg_gbl(&tauRef, 1, "int", OPS_READ),
					 ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
					 ops_arg_gbl(force, spaceDim, "double", OPS_READ),
					 ops_arg_gbl(&noElem, 1, "int", OPS_READ),
					 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
					 ops_arg_gbl(compo.index, 2, "int", OPS_READ));
	}



}

void Psm::InitializeVariables() {

	int size = noElem * spaceDim;

	poros->InitializeVariables();

	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());

		const int blockIndex{block.ID()};
		ops_par_loop(KerInitialize,"KerInitialize", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(Fd[blockIndex], noElem, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_gbl(&size, 1, "int", OPS_READ));

	}
}


void Psm::CalculateDragForce() {


	Real tauRef = compo.tauRef;
	const Real* pdt {pTimeStep()};
	const Real dx {GetDx()}; //TODO ADD to flowfield
	Real xPos[spaceDim], FdLocal[spaceDim], TdLocal[spaceDim];
	int size = noElem * spaceDim;
	int StenList[2 * spaceDim];
	int blockIndex;
	int idParticle;
	for (auto& idBlock : BlockParticleList) {
		BlockParticles& ParticlesCurrentBlock = idBlock.second;
		if ( !ParticlesCurrentBlock.owned ) continue;
		blockIndex = ParticlesCurrentBlock.GetBlock().ID();
		int nLocal = ParticlesCurrentBlock.NParticles
				+ ParticlesCurrentBlock.NPeriodic;
		for (int iPart = 0; iPart < nLocal; iPart++) {

			xPos[0] = ParticlesCurrentBlock.particleList[iPart].xParticle[0];
			xPos[1] = ParticlesCurrentBlock.particleList[iPart].xParticle[1];
			xPos[2] = ParticlesCurrentBlock.particleList[iPart].xParticle[2];

			idParticle = iPart;

			for (int iDir = 0; iDir < spaceDim; iDir++) {
				TdLocal[iDir] = 0.0;
				FdLocal[iDir] = 0.0;
			}

			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
				StenList[iDir] = ParticlesCurrentBlock.particleList[iPart].stenList[iDir];


			ops_par_loop(KerDragPSM,"KerDragPSM", ParticlesCurrentBlock.GetBlock().Get(),
						 spaceDim, StenList,
						 ops_arg_dat(poros->GetIntFieldVariable(0).at(blockIndex),
								 	 noElem, LOCALSTENCIL, "int", OPS_READ),
						 ops_arg_dat(poros->GetRealFieldVariable(0).at(blockIndex),
								 	 noElem, LOCALSTENCIL,"double", OPS_READ),
					     ops_arg_dat(poros->GetRealFieldVariable(2).at(blockIndex),
					    		 	 size, LOCALSTENCIL,"double", OPS_READ),
						 ops_arg_dat(Fd[blockIndex], size, LOCALSTENCIL,
								 	 "double", OPS_READ),
						 ops_arg_gbl(FdLocal, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(TdLocal, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&idParticle, 1, "int", OPS_READ),
						 ops_arg_gbl(xPos, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(pdt, 1, "double", OPS_READ),
						 ops_arg_gbl(&tauRef, 1, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
						 ops_arg_gbl(&noElem, 1, "int", OPS_READ));

			for (int iDir = 0; iDir < spaceDim; iDir++) {
				ParticlesCurrentBlock.particleList[iPart].FDrag[iDir] += FdLocal[iDir] * dx * dx * dx;
				ParticlesCurrentBlock.particleList[iPart].TDrag[iDir] += TdLocal[iDir] * dx * dx * dx;
			}

			printf("Rank %d particle [%f %f %f] Fd = [%e %e %e] Td=[%e %e %e]\n",
					ops_get_proc(), xPos[0], xPos[1], xPos[2],
					FdLocal[0] * dx * dx * dx, FdLocal[1] * dx * dx *dx,
					FdLocal[2] * dx * dx * dx, TdLocal[0] * dx * dx * dx,
					TdLocal[1] * dx * dx * dx, TdLocal[2] * dx * dx * dx);
		}
	}


}

void Psm::MappingFunction(bool flag) {

	if (flag) {
		poros->ParticleProjection();
	}
	else
		poros->UpdateProjection();

}

void Psm::WriteToHdf5(const std::string& caseName, const SizeType timeStep) {

	Fd.WriteToHDF5(caseName, timeStep);

	poros->WriteToHdf5(caseName, timeStep);

}


