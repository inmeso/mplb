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

/*! @brief  Head files for handling a single particle
 * @author C. Tsigginos
 *  @details Define the Particle Class which handles the particle information
 */


#include "particle.h"
#include <stdlib.h>
#include <math.h>
#include "ops_seq_v2.h"

Particle::Particle(int Dimension, Real radius, Real* xp, std::vector<Real> Shape, std::string particleImport) {

	spaceDim = Dimension;
	if (spaceDim < 2 && spaceDim > 3) {
		ops_printf("Error: Space dimensions not consistent\n");
		exit(EXIT_FAILURE);
	}

	xParticle = new Real[spaceDim];
	FDrag = new Real[spaceDim];
	TDrag = new Real[spaceDim];
	stenList = new int[2 * spaceDim];
	xOld = new Real[spaceDim];
	uParticle = new Real[spaceDim];
	omegaParticle = new Real[spaceDim];

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		xParticle[iDim] = xp[iDim];
		xOld[iDim] = 0.0;
	}

	//TODO Assign the Particle Shape
	particleShapeType = particleImport;
	std::string Spherical("spherical");
	std::string Quadratic("quadratic");
	std::string Mesh("mesh");


	if (particleShapeType.compare(Spherical)==0) {
		particleShape = new ParticleShape(radius, Spherical);
	}
	else if (particleShapeType.compare(Quadratic)==0) {
		particleShape = new ParticleShapeQuadratic(radius, Quadratic, Shape);
	}
	else if (particleShapeType.compare(Mesh) == 0) {
		particleShape = new ParticleShapeMesh(radius, Mesh, Shape);
	}
	else {
		ops_printf("Error: This options is not supported\n");
		exit(EXIT_FAILURE);

	}


}

Particle::~Particle() {

	delete[] xParticle;
	delete particleShape;
	delete[] FDrag;
	delete[] TDrag;
	delete[] stenList;
	delete[] omegaParticle;
	delete[] uParticle;
}

void Particle::initializeDrag() {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		FDrag[iDim] = 0.0;
		TDrag[iDim] = 0.0;
	}

}

void Particle::updateParticleLocation(Real* xp) {

	for (int iDim = 0; iDim < spaceDim; iDim++)
		xParticle[iDim] = xp[iDim];
}

void Particle::updateParticleShape(Real radius, std::vector<Real> shape) {

	particleShape->updateShape(radius, shape);
	particleShape->rotate();
}

void Particle::addDrag(Real* Fp, Real* Tp) {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		FDrag[iDim] += Fp[iDim];
		TDrag[iDim] += Tp[iDim];
	}

}

void Particle::evaluateDrag(Real dt) {

	for (int iDim = 0; iDim < spaceDim;iDim++) {
		FDrag[iDim] /= dt;
		TDrag[iDim] /= dt;
	}

}

void Particle::pushDrag(Real* Fd, Real* Td) {

	for (int iDim = 0; iDim < spaceDim ; iDim++) {
		Fd[iDim] = FDrag[iDim];
		Td[iDim] = TDrag[iDim];
	}
}

void Particle::updateXOld(Real *x) {

	for (int iDim =0; iDim < spaceDim; iDim++)
		xOld[iDim] = x[iDim];
}

void Particle::updateParticleVelocities(Real* uP, Real* omP) {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		uParticle[iDim] = uP[iDim];
		omegaParticle[iDim] = omP[iDim];
	}

}

void Particle::updateStencil(Real* xBounds, int *Nf, Real dx) {

	Real xMin[spaceDim];
	int iLocal[spaceDim];
	int iRadius;
	for (int iDim = 0; iDim < spaceDim; iDim++)
		xMin[iDim] = xBounds[2 * spaceDim + 2 * iDim];

	iRadius = (int) ceil(particleShape->Rparticle / dx);
	for (int iDim = 0; iDim < spaceDim; iDim++)
		iLocal[iDim] = (int) floor((xParticle[iDim] - xMin[iDim]) / dx);

	//Identify the particle stencil with no corrections
	for (int iDim = 0; iDim < spaceDim; iDim++) {
		stenList[2 * iDim] = iLocal[iDim] - iRadius;
		stenList[2 * iDim + 1] = iLocal[iDim] + iRadius + 1;
	}

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		if (stenList[2 * iDim] < Nf[2 * iDim])
			stenList[2 * iDim] = Nf[2 * iDim];

		if (stenList[2 * iDim + 1] > Nf[2 * iDim + 1])
			stenList[2 * iDim + 1] = Nf[2 * iDim + 1];
	}



}
