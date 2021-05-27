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

BlockParticles::BlockParticles(int spacedim, Real dx1, Real cutoff, Block myBlock) {

	spaceDim = spacedim;

	if ((spaceDim < 1) && (spaceDim > 3)) {

		ops_printf("ERROR: Spacedim is not consistent");
		exit(EXIT_FAILURE);
	}


	cutOff = cutoff;
	if (cutOff <= 0.0) {
		ops_printf("Error: Cut-off is non-positive number\n");
		exit(EXIT_FAILURE);
	}

	dx = dx1;

	xBoundGlobal = new Real[2 * spaceDim];
	xBoundLocal = new Real[2 * spaceDim];
	Nf = new int[2 * spaceDim];
	owned = false;

	currentBlock = &myBlock;

	NParticles  = 0;
	NPeriodic = 0;

	Nmax = 100;
	particleList.reserve(Nmax);

}

BlockParticles::~BlockParticles() {

	delete[] xBoundGlobal;
	delete[] xBoundLocal;
	delete[] Nf;

}

int BlockParticles::checkDistanceBlock() {

	int nTotal = NParticles + NPeriodic;
	Real distance;
	int cut = 0;
	for (int iPart = 0; iPart < nTotal; iPart) {
		distance  = 0.0;
		for (int iDim = 0; iDim < spaceDim; iDim++)
			distance  += (particleList[iPart].xParticle[iDim] - particleList[iPart].xOld[iPart]) *
			(particleList[iPart].xParticle[iDim] - particleList[iPart].xOld[iPart]);
		if (distance < cutOff * cutOff) {
			cut = 1;
			break;
		}
	}

	return cut;
}

void BlockParticles::initializeDragForce() {

	for (int iPart = 0; iPart < NParticles; iPart++)
		particleList[iPart].initializeDrag();
}

void BlockParticles::evaluateDragForce(Real dt)    {

	for (int iPart = 0; iPart < NParticles; iPart++)
		particleList[iPart].evaluateDrag(dt);
}

void BlockParticles::updateParticle(Real* xpos, Real radius, std::vector<Real> shape) {


	NParticles += 1;
	if (NParticles >= Nmax) {
		Nmax += 100;
		particleList.reserve(Nmax);
	}

	particleList[NParticles-1].updateParticleLocation(xpos);
	particleList[NParticles-1].updateParticleShape(radius, shape);

}
void BlockParticles::updateParticleVelocities(Real* uPart, Real* omegaT) {

	particleList[NParticles - 1].updateParticleVelocities(uPart, omegaT);
}


void BlockParticles::insertParticle(Real* xpos, Real radius, std::vector<Real> shape,
		Real* uPart, Real* omegaT) {

	bool flag = SphereParallepipedIntersection(xpos, radius);

	if (flag == true) {
		updateParticle(xpos, radius, shape);
		updateParticleVelocities(uPart, omegaT);
	}

}


void BlockParticles::findStencil() {

	int nTotal = NPeriodic + NParticles;

	for (int iPart = 0; iPart < nTotal; iPart++)
		particleList[iPart].updateStencil(xBoundGlobal, Nf,  dx);

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

	int start[spaceDim], end[spaceDim], range[spaceDim], disp[spaceDim];
	int iterRng[2 * spaceDim];
	Real xb[spaceDim];
	for (int nDim = 0; nDim < spaceDim; nDim++) {
		xBoundLocal[2 * nDim] = std::numeric_limits<Real>::max();
		xBoundLocal[2 * nDim + 1] = -1.0 * std::numeric_limits<Real>::max();
	}
	int blockIndex = currentBlock->ID();

	if (owned) {
		//Find block ranges

		ops_get_abs_owned_range(currentBlock->Get(), range, start, end, disp);
		for (int iDir = 0; iDir < spaceDim; iDir++) {
			Nf[2 * iDir] = start[iDir];
			Nf[2 * iDir + 1] = end[iDir];
		}

		for (int iDim = 0; iDim < spaceDim; iDim++) {
			iterRng[2 * iDim] = Nf[2 *iDim];
			iterRng[2 * iDim + 1] = Nf[2 * iDim + 1] + 1;
		}
		int spacedim = spaceDim;
		ops_par_loop(KerCarBound,"KerCarBound", currentBlock->Get(), spaceDim, iterRng,
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



		ops_par_loop(KerCarBound, "KerCarBound", currentBlock->Get(), spaceDim, iterRng,
				ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
							LOCALSTENCIL, "double", OPS_READ),
				ops_arg_gbl(xb, spaceDim, "double", OPS_READ),
				ops_arg_gbl(&spacedim, 1, "int", OPS_READ));

		for (int iDim = 0; iDim < spaceDim; iDim++)
			xBoundLocal[2 * iDim + 1] = xb[iDim] + 0.5 * dx;




	}

}

void BlockParticles::getLocalBound(Real* xBound) {

	for (int iDim = 0; iDim <2 * spaceDim; iDim++)
		xBound[iDim] = xBoundLocal[iDim];


}

void BlockParticles::setGlobalBound(Real* xBound) {

	for (int iDim = 0; iDim < 2* spaceDim; iDim++)
		xBoundGlobal[iDim] = 2 * xBound[iDim];
}

void BlockParticles::getOwnership(bool flag) {

	owned = flag;
}
