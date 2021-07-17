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
 * @brief   Functions for the implementation of the partially saturated method
 * @author  C. Tsigginos
 * @details Functions for the implementation of the partially saturated method
 */

#include "psm.h"

#include "flowfield_host_device.h"
#include "ops_seq_v2.h"
#include "force_fsi.h"

#include "psm_kernel.inc"


void FsiInitializePSM(std::shared_ptr<FpiData>& fpiPtr) {

	int spaceDim = fpiPtr->GetSpaceDim();
	int size = fpiPtr->SizeAtRealType(0);

	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
		const int blockIndex{block.ID()};
		ops_par_loop(KerInitializePSM,"KerInitializePSM", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(fpiPtr->GetRealField(0, blockIndex), size,
							 	 LOCALSTENCIL,"double", OPS_WRITE),
					 ops_arg_gbl(&size, 1, "int", OPS_READ));

	}

}

void FsiForcePSM(std::shared_ptr<FpiData>& fpiPtr) {};

void FsiCollisionsPSM(std::shared_ptr<FpiData>& fpiPtr) {

	int spaceDim, noElem, componentID,size;
	int sizeR1, sizeR2, sizeI1;
	spaceDim = fpiPtr->GetSpaceDim();

	Real tauRef;
	Component compo{fpiPtr->FluidComponent()}; //TO IF ICANT I need to

	Forces forceParams{ForceModel()};
	int sizeForce = forceParams.forceDefinition.size();
	Real force[sizeForce];

	for (int iDir = 0; iDir < sizeForce; iDir++)
		force[iDir] = forceParams.forceDefinition[iDir];

	int forceFlag = (int) forceParams.model;


	tauRef = compo.tauRef;
	componentID = compo.id;
	const Real* pdt{pTimeStep()};
	noElem = fpiPtr->GetNElem();
	size = fpiPtr->SizeAtRealType(0);
	std::shared_ptr<MappingParticles> mappingModel = GetMappingModels().at(componentID);
	sizeR1 = mappingModel->SizeAtRealType(0);
	sizeR2 = mappingModel->SizeAtRealType(1);
	sizeI1 = mappingModel->SizeAtIntType(0);

	for (const auto& idBlock : g_Block()) {
		const Block& block{idBlock.second};
		std::vector<int> iterRng;
		iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
		const int blockIndex{block.ID()};
		ops_par_loop(KerCollisionPSM, "KerCollisionPSM", block.Get(),
					 spaceDim, iterRng.data(),
					 ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
							 	 "double", OPS_WRITE),
					 ops_arg_dat(fpiPtr->GetRealField(0, blockIndex), size, LOCALSTENCIL,
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
					 ops_arg_dat(mappingModel->GetIntField(0, blockIndex), sizeI1,
							 	 LOCALSTENCIL, "int", OPS_READ),
					 ops_arg_dat(mappingModel->GetRealField(0, blockIndex), sizeR1,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_dat(mappingModel->GetRealField(1, blockIndex), sizeR2,
							 	 LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_gbl(pdt, 1, "double", OPS_READ),
					 ops_arg_gbl(&tauRef, 1, "int", OPS_READ),
					 ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
					 ops_arg_gbl(force, sizeForce, "double", OPS_READ),
					 ops_arg_gbl(&sizeForce, 1, "int", OPS_READ),
					 ops_arg_gbl(&noElem, 1, "int", OPS_READ),
					 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
					 ops_arg_gbl(compo.index, 2, "int", OPS_READ));
	}


}

void CalculateDragPSM(std::shared_ptr<FpiData>& fpiPtr) {

	int spaceDim, size, noElem, componentID;
	int blockIndex, idParticle;
	int sizeI1, sizeR1, sizeR3;
	spaceDim = fpiPtr->GetSpaceDim();
	int stencil[2*spaceDim];
	Real tauRef, dx, xPos[spaceDim], FdLocal[spaceDim], TdLocal[spaceDim];
	const Real* pdt{pTimeStep()};
	noElem = fpiPtr->GetNElem();
	size = fpiPtr->SizeAtRealType(0);
	Component compo{fpiPtr->FluidComponent()};
	tauRef = compo.tauRef;
	componentID = compo.id;
	dx = GetDx();

	//GetPorosModel
	std::shared_ptr<MappingParticles> mappingModel = GetMappingModels().at(componentID);
	sizeR1 = mappingModel->SizeAtRealType(0);
	sizeR3 = mappingModel->SizeAtRealType(2);
	sizeI1 = mappingModel->SizeAtIntType(0);

	for (auto& idBlock : BlockParticleList) {
		BlockParticles& particlesCurrentBlock = idBlock.second;
		if (!particlesCurrentBlock.owned) continue;
		blockIndex = particlesCurrentBlock.GetBlock().ID();
		int nLocal = particlesCurrentBlock.NParticles
				+ particlesCurrentBlock.NPeriodic;
		for (int iPart = 0; iPart < nLocal; iPart++) {
			for (int iDir = 0; iDir < spaceDim; iDir++)  {
				xPos[iDir] = particlesCurrentBlock.particleList[iPart].xParticle[iDir];
				TdLocal[iDir] = 0.0;
				FdLocal[iDir] = 0.0;
			}

			idParticle = iPart;
			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
					stencil[iDir] = particlesCurrentBlock.particleList[iPart].stenList[iDir];

			ops_par_loop(KerDragPSM, "KerDragPSM", particlesCurrentBlock.GetBlock().Get(),
					 	 spaceDim, stencil,
						 ops_arg_dat(mappingModel->GetIntField(0, blockIndex), sizeI1,
								 	 LOCALSTENCIL, "int", OPS_READ),
						 ops_arg_dat(mappingModel->GetRealField(0, blockIndex), sizeR1,
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(mappingModel->GetRealField(2, blockIndex), sizeR3,
						 	 	 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(fpiPtr->GetRealField(0, blockIndex), size,
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(FdLocal, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(TdLocal, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&idParticle, 1, "int", OPS_READ),
						 ops_arg_gbl(xPos, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(pdt, 1, "double", OPS_READ),
						 ops_arg_gbl(&tauRef, 1, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
						 ops_arg_gbl(&noElem, 1, "int", OPS_READ));


			for (int iDir = 0; iDir < spaceDim; iDir++) {
				particlesCurrentBlock.particleList[iPart].FDrag[iDir] += FdLocal[iDir] * dx * dx * dx;
				particlesCurrentBlock.particleList[iPart].TDrag[iDir] += TdLocal[iDir] * dx * dx * dx;
			}
		}
	}

}


