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
 * @brief   Wrapper functions for the initialization of MappingParticles ops variables
 * @author  C. Tsigginos
 * @details wrapper for the initialization of ops dat lists
 */


#include "particle_mapping.h"
#include "grid_mapping.h"
#include "flowfield.h"
#include "scheme.h"
#include "particle_mapping_kernel.inc"
#include <memory.h>
#include "block_particles.h"
#include "grid_mapping.h"
void InitializePorousLists(std::shared_ptr<MappingParticles>& mappingPtr, int component) {

	int sizeReal[3], sizeInt[1];
	int spaceDim = mappingPtr->GetSpaceDim();
	for (int iType = 0; iType < 3; iType++)
		sizeReal[iType] = mappingPtr->SizeAtRealType(iType);

	sizeInt[0] = mappingPtr->SizeAtIntType(0);

	for (const auto& idBlock : BlockParticleList) {
		BlockParticles ParticleCurrentBlock = idBlock.second;
		std::vector<int> iterRng;
		iterRng.assign(ParticleCurrentBlock.GetBlock().WholeRange().begin(),
					   ParticleCurrentBlock.GetBlock().WholeRange().end());

		const int blockIndex{ParticleCurrentBlock.GetBlock().ID()};
		ops_par_loop(KerInitializePorousGrid, "KerInitializePorousGrid",
					 ParticleCurrentBlock.GetBlock().Get(), spaceDim, iterRng.data(),
					 ops_arg_dat(mappingPtr->GetRealField(0, blockIndex), sizeReal[0],
							 	 LOCALSTENCIL, "double", OPS_WRITE),
					 ops_arg_dat(mappingPtr->GetRealField(1, blockIndex), sizeReal[1],
							 	 LOCALSTENCIL, "double", OPS_WRITE),
				     ops_arg_dat(mappingPtr->GetRealField(2, blockIndex), sizeReal[2],
				    		 	 LOCALSTENCIL, "double", OPS_WRITE),
					 ops_arg_dat(mappingPtr->GetIntField(0, blockIndex), sizeInt[0],
							 	 LOCALSTENCIL, "int", OPS_WRITE),
					 ops_arg_gbl(sizeReal, 3, "int", OPS_READ),
					 ops_arg_gbl(sizeInt, 1, "int", OPS_READ));
	}



}

/*
 * 	for (const auto& idBlock : BlockParticleList) {
		BlockParticles particleCurrentBlock = idBlock.second;
		if (!particleCurrentBlock.owned) continue;
		int blockIndex = particleCurrentBlock.GetBlock().ID();
		nlocal = particleCurrentBlock.NParticles +
				 particleCurrentBlock.NPeriodic;


	}
 */
