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
 * @brief  Wrapper functions for mapping particles into LBM grid via the subcell method
 * @author  C. Tsigginos
 * @details The scheme is used in conjuction with the PSM or Prati schemes.
*/
#include "grid_mapping.h"
#include "block_particles.h"

#include "flowfield.h"
#include "flowfield_host_device.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "grid_mapping_kernel.inc"

void ParticleProjectionSphereGrid(std::shared_ptr<MappingParticles>& mappingPtr,
		int component) {

	int sizeReal[3], sizeInt, idParticle, blockIndex, nPoints, nelem;
	int spaceDim = mappingPtr->GetSpaceDim();
	int stenList[2 * spaceDim];
	Real Dx, xPos[spaceDim], Radius, uPart[spaceDim], omPart[spaceDim];
	for (int iDim = 0; iDim < 3; iDim++)
		sizeReal[iDim] = mappingPtr->SizeAtRealType(iDim);
	Dx = GetDx();
	nPoints = 9;
	sizeInt = mappingPtr->SizeAtIntType(0);
	nelem = mappingPtr->NumberOfElements();

	for (const auto& idBlock : BlockParticleList) {
		BlockParticles particleCurrentBlock = idBlock.second;
		if (!particleCurrentBlock.owned) continue;
		blockIndex = particleCurrentBlock.GetBlock().ID();
		int nlocal =  particleCurrentBlock.NParticles +
					  particleCurrentBlock.NPeriodic;

		for (int iPart = 0; iPart < nlocal; iPart++) {

			for (int iDim = 0; iDim < spaceDim; iDim++) {
				xPos[iDim] = particleCurrentBlock.particleList.at(iPart).xParticle[iDim];
				uPart[iDim] = particleCurrentBlock.particleList.at(iPart).uParticle[iDim];
				omPart[iDim] = particleCurrentBlock.particleList.at(iPart).omegaParticle[iDim];
			}
			idParticle = iPart;
			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
				stenList[iDir] = particleCurrentBlock.particleList.at(iPart).stenList[iDir];

			ops_par_loop(KerMappingSphericalGrid,"KerMappingSphericalGrid",
					     particleCurrentBlock.GetBlock().Get(), spaceDim,  stenList,
						 ops_arg_dat(mappingPtr->GetIntField(0, blockIndex), sizeInt,
								 	 LOCALSTENCIL,"int", OPS_RW),
						 ops_arg_dat(mappingPtr->GetRealField(0, blockIndex), sizeReal[0],
									 LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(mappingPtr->GetRealField(1, blockIndex), sizeReal[1],
									 LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(mappingPtr->GetRealField(2, blockIndex), sizeReal[2],
								     LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(xPos, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&Radius, 1, "double", OPS_READ),
						 ops_arg_gbl(&idParticle, 1, "int", OPS_READ),
						 ops_arg_gbl(uPart, spaceDim , "double", OPS_READ),
						 ops_arg_gbl(omPart, spaceDim, "double", OPS_READ),
					     ops_arg_gbl(&Dx, 1, "double", OPS_READ),
					     ops_arg_gbl(&nPoints, 1, "double", OPS_READ),
						 ops_arg_gbl(&nelem, 1, "double", OPS_READ),
					     ops_arg_gbl(&spaceDim, 1, "double", OPS_READ));

		}
	}

}

void UpdateProjectectionSphereGrid(std::shared_ptr<MappingParticles>& mappingPtr,
		   int component) {

	int sizeReal[3], sizeInt, idParticle, blockIndex, nPoints, nelem;
	int spaceDim = mappingPtr->GetSpaceDim();
	int stenList[2 * spaceDim];
	Real dx, xPar[spaceDim], radius, velP[spaceDim], omP[spaceDim];
	dx = GetDx();

	for (int iDim = 0; iDim < 3; iDim++)
		sizeReal[iDim] = mappingPtr->SizeAtRealType(iDim);
	sizeInt = mappingPtr->SizeAtIntType(0);
	nelem = mappingPtr->NumberOfElements();

	for (auto& idBlock : BlockParticleList) {
		BlockParticles particleCurrentBlock = idBlock.second;
		if ( !particleCurrentBlock.owned ) continue;
		int blockIndex = particleCurrentBlock.GetBlock().ID();
		int nlocal =  particleCurrentBlock.NParticles +
				particleCurrentBlock.NPeriodic;

		for (int iPart = 0; iPart < nlocal; iPart++) {
			for (int iDim = 0; iDim < spaceDim; iDim++) {
				xPar[iDim] = particleCurrentBlock.particleList.at(iPart).xParticle[iDim];
				velP[iDim] = particleCurrentBlock.particleList.at(iPart).uParticle[iDim];
				omP[iDim] = particleCurrentBlock.particleList.at(iPart).omegaParticle[iDim];
			}
			radius = particleCurrentBlock.particleList.at(iPart).particleShape->Rparticle;
			idParticle = iPart;
			for (int iDir = 0; iDir < 2 * spaceDim; iDir++)
				stenList[iDir] = particleCurrentBlock.particleList.at(iPart).stenList[iDir];

			ops_par_loop(KerSolidVelocityUpdateGridSphere, "KerSolidVelocityUpdateGridSphere",
						 particleCurrentBlock.GetBlock().Get(), spaceDim, stenList,
						 ops_arg_dat(mappingPtr->GetRealField(1, blockIndex), sizeReal[1],
									 LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(mappingPtr->GetIntField(0, blockIndex), sizeInt,
								 	 LOCALSTENCIL, "int", OPS_READ),
						 ops_arg_dat(mappingPtr->GetRealField(2, blockIndex), sizeReal[2],
								 	 LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(&idParticle, 1, "int", OPS_READ),
						 ops_arg_gbl(xPar, spaceDim, "double", OPS_READ),
					     ops_arg_gbl(&radius, 1, "double", OPS_READ),
						 ops_arg_gbl(velP, spaceDim, "double", OPS_READ),
					     ops_arg_gbl(omP, spaceDim, "double", OPS_READ),
						 ops_arg_gbl(&dx, 1, "double", OPS_READ),
						 ops_arg_gbl(&spaceDim, 1, "int", OPS_READ),
						 ops_arg_gbl(&nelem, 1, "int", OPS_READ));

		}
	}

}
