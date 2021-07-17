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
 * @brief   Wrapper function for evaluating local box bound
 * @author  C. Tsigginos
 * @details Wrapper for the kernel that evaluates the local box bound
 * 			outside the BlockParticle Class
 */
#include "block_particle_helper.h"
#include <stdlib.h>
#include <limits>

#include "block_particle_helper_kernel.inc"
#include "ops_seq_v2.h"

void DefineLocalBoxBound() {

	for (auto& blockParticle : BlockParticleList) {
		if (blockParticle.second.owned) {

			int spaceDim = SpaceDim();
			int start[spaceDim], end[spaceDim], range[2 * spaceDim], disp[spaceDim];
			int Nf[2 * spaceDim]; //TODO ADD IT TO BlockPArticles
			Real xBoundLocal[2 * spaceDim];
			int size = 2 * spaceDim;
			int iterRng[size];
			Real xb[spaceDim];
			Real dx = GetDx();
			for (int nDim = 0; nDim < spaceDim; nDim++) {
				xBoundLocal[2 * nDim] = std::numeric_limits<Real>::max();
				xBoundLocal[2 * nDim + 1] = -1.0 * std::numeric_limits<Real>::max();
			}

			int blockIndex =blockParticle.second.GetBlock().ID();
			std::vector<int> iterRng1;
			iterRng1.assign(blockParticle.second.GetBlock().WholeRange().begin(),
							blockParticle.second.GetBlock().WholeRange().end());

			for (int iDim = 0; iDim < 2 * spaceDim; iDim++)
				range[iDim] = iterRng1.at(iDim);

			ops_get_abs_owned_range(blockParticle.second.GetBlock().Get(), range, start,
					 	 	 	 	 end, disp);

			for (int iDir = 0; iDir < spaceDim; iDir++) {
				Nf[2 * iDir] = start[iDir];
				Nf[2 * iDir + 1] = end[iDir];
			}



			for (int iDim = 0; iDim < spaceDim; iDim++) {
				iterRng[2 * iDim] = Nf[2 *iDim];
				iterRng[2 * iDim + 1] = Nf[2 * iDim] + 1;
			}



			ops_par_loop(KerCarBound,"KerCarBound", blockParticle.second.GetBlock().Get(),
						 spaceDim, iterRng,
						 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
								 	 LOCALSTENCIL, "double", OPS_READ),
						  ops_arg_gbl(xb, spaceDim, "double", OPS_READ),
						  ops_arg_gbl(&spaceDim, 1, "int", OPS_READ));

			for (int iDim = 0; iDim < spaceDim; iDim++) {
				xBoundLocal[2 * iDim] = xb[iDim] - 0.5 * dx;
			}


			//2nd point
			for (int iDir = 0; iDir < spaceDim; iDir++) {
				iterRng[2 * iDir] = Nf[2 * iDir + 1] - 1;
				iterRng[2 * iDir + 1] = Nf[2 * iDir + 1];
			}



			ops_par_loop(KerCarBound,"KerCarBound", blockParticle.second.GetBlock().Get(),
					 	 spaceDim, iterRng,
						 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(xb, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ));

			for (int iDim = 0; iDim < spaceDim; iDim++)
				xBoundLocal[2 * iDim + 1] = xb[iDim] + 0.5 * dx;

			blockParticle.second.SetLocalBound(xBoundLocal);
			blockParticle.second.SetNfLocal(Nf);
		}



	}
#ifdef CPU
#if DebugLevel >= 2
	Real xBound[6];
	for (auto& particleBlock : BlockParticleList) {
		particleBlock.second.GetLocalBound(xBound);
		printf("Rank %d, Block %d xLocal =[%f %f] x[ %f %f] x [%f %f]\n",
					ops_get_proc(), particleBlock.second.GetBlock().ID(),
					xBound[0], xBound[1], xBound[2], xBound[3],
					xBound[4], xBound[5]);
	}
#endif
#endif

}


