
/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]s
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

/*! @brief   Define functions and variables for handling particles
 * @author  C. Tsigginos
 *
 * @details: Functions for handling particle data, function for iterating functions
 * 			 the correct DEM-LBM scheme.
 */

#include "dem_handle.h"
#include "flowfield_host_device.h"

std::map<int, BlockParticles> BlockParticleList;

void DefineBlockParticles(const BlockGroup& blocks, int spacedim, Real cutoff, Real dx1) {

	for (const auto& idBlock : blocks) {
		const Block& myBlock{idBlock.second};
		BlockParticles ParticleNew(spacedim, dx1, cutoff,myBlock);
		BlockParticleList.emplace(idBlock.first, ParticleNew);
	}


}

void findBlockOwnership() {

	 for (auto& idBlock : BlockParticleList) {

		 const int blockIdx{idBlock.first};
		 bool flag{true};
#ifdef OPS_MPI
		 sub_block_list sb = OPS_sub_block_list[blockIdx];
		 	if (sb->owned)
		 		flag = true;
		 	else
		 		flag = false;

#endif
		 	idBlock.second.getOwnership(flag);
	 }


}

void DefineGlobalBlockBox() {
	int spaceDim1 = SpaceDim();
	Real xMinMaxTmp[spaceDim1];
	Real xb[2 * spaceDim1], xMin[spaceDim1], xMax[spaceDim1], xTemp[2 * spaceDim1];
	for (auto& idBlock : BlockParticleList) {
		BlockParticles& Particles = idBlock.second;
		const int blockIdx{idBlock.first};
		Particles. getLocalBound(xb);

		for (int iDim = 0; iDim < spaceDim1; iDim++) {
			xMin[iDim] = xb[2 * iDim];
			xMax[iDim] = xb[2 * iDim + 1];
		}

#ifdef OPS_MPI


		MPI_Allreduce(xMin, xMinMaxTmp, spaceDim1, MPI_DOUBLE, MPI_MIN, OPS_MPI_GLOBAL);

		for (int iDim = 0; iDim < spaceDim1; iDim++) {
			xTemp[2 * iDim] = xMinMaxTmp[iDim];
		}

		MPI_Allreduce(xMin, xMinMaxTmp, spaceDim1, MPI_DOUBLE, MPI_MAX, OPS_MPI_GLOBAL);
		for (int iDim = 0; iDim < spaceDim1; iDim++) {
			xTemp[2 * iDim + 1] = xMinMaxTmp[iDim + 1];
		}


#else
		for (int iDim = 0; iDim < 2 * spaceDim1; iDim++)
			xTemp[iDim] = xb[iDim];
#endif
		Particles.setGlobalBound(xTemp);


	}
}
