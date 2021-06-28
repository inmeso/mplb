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
 * @brief   Class for handling particles at various blocks
 * @author  C. Tsigginos
 * @details: Handling of particles that owned by a given block
 */

#include "block_particles.h"

#include <stdlib.h>
#include <limits>
#include "ops_seq_v2.h"
#include "block_particles.inc"

ListBlockParticles BlockParticleList;

BlockParticles::BlockParticles(int spacedim, Real dx1, Real cutoff,const Block& myBlock,
		ParticleShapeDiscriptor particleType, int nInputExtraVariables,
		int nOutputExtraVariables) : currentBlock(myBlock){

	spaceDim = spacedim;
	int size = 2 * spaceDim;
//	xBoundGlobal = new Real [size];
//	xBoundLocal = new Real [size];

	owned = true;
	dx = dx1;

	if ((spaceDim < 1) && (spaceDim > 3)) {

		ops_printf("ERROR: Spacedim is not consistent");
		exit(EXIT_FAILURE);
	}


	cutOff = cutoff;
	if (cutOff <= 0.0) {
		ops_printf("Error Particle object: Cut-off %12.9e is non-positive number\n", cutOff);
		exit(EXIT_FAILURE);
	}


	if (dx <= 0.0) {
		ops_printf("Error BlockParticles: Grid size is not positive\n");
		exit(EXIT_FAILURE);
	}





	NParticles  = 0;
	NPeriodic = 0;
	particleShape = particleType;
	//sanity checks
	switch (particleShape) {
		case spherical:
			ops_printf("Spherical particles will be generated for Block %d\n", myBlock.ID());
			break;
		case quadratic:
			ops_printf("Quadratic particles will be generated for Block %d\n", myBlock.ID());
			break;
		case mesh:
			ops_printf("Polygonal particles will be generated for Block %d\n", myBlock.ID());
			break;
		default:
			ops_printf("This types of particles is not supported.\n Exiting\n");
			exit(EXIT_FAILURE);
	}

	Nmax = 0;
	NReserve = 100;
	particleList.reserve(NReserve);

	if (nInputExtraVariables > 0) {
		 hasExtraInputVariables = true;
		 extraInputSize = nInputExtraVariables;
	}
	else {
		 hasExtraInputVariables = false;
		 extraInputSize = 0;
	}

	if (nOutputExtraVariables > 0) {
		 hasExtraOutputVariables = true;
		 extraOutputSize = nOutputExtraVariables;
	}
	else {
		hasExtraOutputVariables = false;
		extraOutputSize = 0;
	}

}

BlockParticles::~BlockParticles() {

//	delete[] xBoundGlobal;
//	delete[] xBoundLocal;
//	delete[] Nf;

}

int BlockParticles::CheckDistanceBlock() {

	int nTotal = NParticles + NPeriodic;
	Real distance;
	int cut = 0;
	for (int iPart = 0; iPart < nTotal; iPart++) {
		distance  = 0.0;
		for (int iDim = 0; iDim < spaceDim; iDim++)
			distance  += (particleList[iPart].xParticle[iDim] - particleList[iPart].xOld[iDim]) *
			(particleList[iPart].xParticle[iDim] - particleList[iPart].xOld[iDim]);

		if (distance > cutOff * cutOff) {
			cut = 1;
			break;
		}
	}
	return cut;
}

void BlockParticles::UpdateOldParticlePosition() {

	int nTotal = NParticles + NPeriodic;
	for (int iPart = 0; iPart < nTotal; iPart++) {
		for (int iDim = 0; iDim < spaceDim; iDim++)
			particleList.at(iPart).UpdateOldParticlePositions();
	}
}

void BlockParticles::InitializeDragForce() {

	for (int iPart = 0; iPart < NParticles; iPart++)
		particleList[iPart].InitializeDrag();
}

void BlockParticles::EvaluateDragForce(Real dt)    {

	for (int iPart = 0; iPart < NParticles; iPart++)
		particleList[iPart].EvaluateDrag(dt);
}

int BlockParticles::UpdateParticle(Real* xpos, Real radius, std::vector<Real> shape) {


	NParticles += 1;
	if (NParticles > Nmax) {
		Nmax += 1;
		particleList.push_back(Particle(spaceDim, radius, xpos, shape, particleShape,
				extraInputSize, extraOutputSize));
		if (Nmax > NReserve) {
			NReserve = Nmax + 100;

#ifdef CPU
#if DebugLevel > 0
			printf("Rank %d at Block %d reserves totally %d particles\n", ops_get_proc(), currentBlock.ID(), Nmax+100);
#endif
#endif
			particleList.reserve(NReserve);
		}
	}
	else {
		particleList[NParticles-1].UpdateParticleLocation(xpos);
		particleList[NParticles-1].UpdateParticleShape(radius, shape);
	}

	return NParticles-1;

}
void BlockParticles::UpdateParticleVelocities(Real* uPart, Real* omegaT) {

	particleList.at(NParticles - 1).UpdateParticleVelocities(uPart, omegaT);
}


int BlockParticles::InsertParticle(Real* xpos, Real radius, std::vector<Real> shape,
		Real* uPart, Real* omegaT, std::vector<Real>& inputData) {

	int idParticle;
	bool flag = SphereParallepipedIntersection(xpos, radius);

	if (flag) {
#ifdef CPU
#if DebugLevel >=2
		printf("Rank %d: Particle [%f %f %f] will be inserted in the particleList\n",
				ops_get_proc(), xpos[0], xpos[1], xpos[2]);
#endif
#endif
		idParticle = UpdateParticle(xpos, radius, shape);
		UpdateParticleVelocities(uPart, omegaT);
		if (hasExtraInputVariables)
			GetAdditionalInputVariables(inputData, NParticles-1);
		return idParticle;
	}


	return -1;
}


void BlockParticles::FindStencil() {

	int nTotal = NPeriodic + NParticles;

	for (int iPart = 0; iPart < nTotal; iPart++)
		particleList[iPart].UpdateStencil(xBoundGlobal, Nf,  dx);

}

bool BlockParticles:: SphereParallepipedIntersection(Real* xpos,Real radius) {

	Real testX[spaceDim];
	Real distance;




	for (int iDim = 0; iDim < spaceDim; iDim++) {

		testX[iDim] = xpos[iDim];
		if (xpos[iDim] < xBoundLocal[2 * iDim])
			testX[iDim] = xBoundLocal[2 * iDim];
		else if (xpos[iDim] > xBoundLocal[2 * iDim + 1])
			testX[iDim] = xBoundLocal[2 * iDim + 1];
	}

	distance = 0.0;
	for (int iDim = 0; iDim < spaceDim; iDim++)
		distance += (xpos[iDim] - testX[iDim]) * (xpos[iDim] - testX[iDim]);

	if (distance <= radius * radius)
		return true;

	return false;

}

void BlockParticles::FindBoxLocalBound() {

	int start[spaceDim], end[spaceDim], range[2 * spaceDim], disp[spaceDim];
	int size = 2 * spaceDim;
	int iterRng[size];
	Real xb[spaceDim];
	for (int nDim = 0; nDim < spaceDim; nDim++) {
		xBoundLocal[2 * nDim] = std::numeric_limits<Real>::max();
		xBoundLocal[2 * nDim + 1] = -1.0 * std::numeric_limits<Real>::max();
	}
	int blockIndex = currentBlock.ID();

	if (owned) {
		//Find block ranges
		std::vector<int> iterRng1;
		iterRng1.assign(currentBlock.WholeRange().begin(), currentBlock.WholeRange().end());

		for (int iDim = 0; iDim < 2 * spaceDim; iDim++)
			range[iDim] = iterRng1.at(iDim);

		 ops_get_abs_owned_range(currentBlock.Get(), range, start, end, disp);



		for (int iDir = 0; iDir < spaceDim; iDir++) {
			Nf[2 * iDir] = start[iDir];
			Nf[2 * iDir + 1] = end[iDir];
		}



		for (int iDim = 0; iDim < spaceDim; iDim++) {
			iterRng[2 * iDim] = Nf[2 *iDim];
			iterRng[2 * iDim + 1] = Nf[2 * iDim] + 1;
		}


		int spacedim = spaceDim;
		ops_par_loop(KerCarBound,"KerCarBound", currentBlock.Get(), spaceDim, iterRng,
				ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
							LOCALSTENCIL, "double", OPS_READ),
				ops_arg_gbl(xb, spaceDim, "double", OPS_READ),
				ops_arg_gbl(&spacedim, 1, "int", OPS_READ));

		for (int iDim = 0; iDim < spaceDim; iDim++) {
			xBoundLocal[2 * iDim] = xb[iDim] - 0.5 * dx;
		}


		//2nd point
		for (int iDir = 0; iDir < spaceDim; iDir++) {
			iterRng[2 * iDir] = Nf[2 * iDir + 1] - 1;
			iterRng[2 * iDir + 1] = Nf[2 * iDir + 1];
		}



		ops_par_loop(KerCarBound, "KerCarBound", currentBlock.Get(), spaceDim, iterRng,
				ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
							LOCALSTENCIL, "double", OPS_READ),
				ops_arg_gbl(xb, spaceDim, "double", OPS_READ),
				ops_arg_gbl(&spacedim, 1, "int", OPS_READ));

		for (int iDim = 0; iDim < spaceDim; iDim++)
			xBoundLocal[2 * iDim + 1] = xb[iDim] + 0.5 * dx;
#ifdef CPU
#if DebugLevel >= 2
		for (int iDim = 0; iDim < spaceDim; iDim++)
			printf("Block %d: Current local Bound [%f %f]\n",currentBlock.ID(), xBoundLocal[2*iDim], xBoundLocal[2 * iDim+1]);
#endif
#endif
	}

}


void BlockParticles::ExtractBound(Real* xMin, Real* xMax) {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		xMin[iDim] = xBoundLocal[2 * iDim];
		xMax[iDim] = xBoundLocal[2 * iDim + 1];
	}


}

void BlockParticles::GetLocalBound(Real* xBound) {

	for (int iDim = 0; iDim <2 * spaceDim; iDim++)
		xBound[iDim] = xBoundLocal[iDim];


}

void BlockParticles::SetGlobalBound(Real* xBound) {

	for (int iDim = 0; iDim < 2* spaceDim; iDim++)
		xBoundGlobal[iDim] =  xBound[iDim];
#ifdef CPU
#if DebugLevel >= 2
	printf("Block %d at rank %d Global box ",currentBlock.ID(), ops_get_proc());
	for (int iDim = 0; iDim < spaceDim; iDim++)
		printf("[%f %f] ", xBoundGlobal[2 * iDim], xBoundGlobal[2* iDim + 1]);
#endif
#endif
}

void BlockParticles::GetOwnership() {

	owned = true;
#ifdef OPS_MPI
	int blockIndex = currentBlock.ID();

	 sub_block_list sb = OPS_sub_block_list[blockIndex];
	 	if (sb->owned)
	 		owned = true;
	 	else
	 		owned = false;
#endif
}

void BlockParticles::ExtractDragForce(Real* Fd, Real* Td,int iParticle) {

	particleList.at(iParticle).PushDrag(Fd, Td);

}

void BlockParticles::ExtractPositions(Real* xPos,int iParticle) {

	particleList.at(iParticle).GetParticlePositions(xPos);
}

void BlockParticles::GetAdditionalOutputVariables(std::vector<Real>& output,int iParticle) {

	particleList.at(iParticle).GetOutputVariables(output);
}

void BlockParticles::GetAdditionalInputVariables(std::vector<Real>& input,int idParticle) {

	particleList.at(idParticle).SetInputVariables(input);
}
