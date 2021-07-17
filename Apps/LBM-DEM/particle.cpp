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
 *  @details Class for handling the data stored per given particle
 */


#include "particle.h"
#include <stdlib.h>
#include <math.h>
#include "ops_seq_v2.h"
#include <new>

Particle::Particle(int Dimension, Real radius, Real* xp, std::vector<Real> Shape,
		ParticleShapeDiscriptor particleImport, int inputVariablesSize,
		int outputVariableSize) {

	spaceDim = Dimension;
	if (spaceDim < 2 && spaceDim > 3) {
		ops_printf("Error: Space dimensions not consistent\n");
		exit(EXIT_FAILURE);
	}

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		xParticle[iDim] = xp[iDim];
		xOld[iDim] = 0.0;
	}

	//TODO Assign the Particle Shape
	particleShapeType = particleImport;

	switch(particleShapeType) {
		case spherical:
			particleShape.reset(new ParticleShape(radius, spherical));
			break;
		case quadratic:
			particleShape.reset(new ParticleShapeQuadratic(radius, quadratic, Shape));
			break;
		case mesh:
			particleShape.reset(new ParticleShape(radius, mesh)); //TODO ADD ACTUAL MESH PARTICLE
			break;
		default:
			ops_printf("ERROR: This type of particle type is not supported\n");
			exit(EXIT_FAILURE);
	}

	//Assing extra variables
	if (inputVariablesSize > 0) {
		nInputExtraVariables = inputVariablesSize;
		inputVariables.reserve(nInputExtraVariables);
	}
	else
		nInputExtraVariables = 0;

	if (outputVariableSize > 0) {
		nOutputExtraVariables = outputVariableSize;
		outputVariables.reserve(nOutputExtraVariables);
	}
	else
		nOutputExtraVariables = 0;


}

Particle::~Particle() {


}

void Particle::InitializeDrag() {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		FDrag[iDim] = 0.0;
		TDrag[iDim] = 0.0;
	}

}

std::string Particle::GetParticleShape() {

	switch (particleShapeType) {
		case spherical:
			return "sphere";
		case quadratic:
			return "quadratic";
		case mesh:
			return "mesh";

	}

	return "none";
}

void Particle::UpdateParticleLocation(Real* xp) {

	for (int iDim = 0; iDim < spaceDim; iDim++)
		xParticle[iDim] = xp[iDim];
}

void Particle::UpdateParticleShape(Real radius, std::vector<Real> shape) {

	particleShape->UpdateShape(radius, shape);
	particleShape->Rotate();
}

void Particle::AddDrag(Real* Fp, Real* Tp) {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		FDrag[iDim] += Fp[iDim];
		TDrag[iDim] += Tp[iDim];
	}

}

void Particle::EvaluateDrag(Real dt) {

	for (int iDim = 0; iDim < spaceDim;iDim++) {
		FDrag[iDim] /= dt;
		TDrag[iDim] /= dt;
	}

	printf("Rank: %d [%f %f %f] Fd=[%12.9e %12.9e %12.9e] Td=[%12.9e %12.9e %12.9e]\n",
			ops_get_proc(), xParticle[0], xParticle[1], xParticle[2],
			FDrag[0], FDrag[1], FDrag[2], TDrag[0], TDrag[1], TDrag[2]);
}

void Particle::PushDrag(Real* Fd, Real* Td) {

	for (int iDim = 0; iDim < spaceDim ; iDim++) {
		Fd[iDim] = FDrag[iDim];
		Td[iDim] = TDrag[iDim];
	}
}

void Particle::UpdateOldParticlePositions() {

	for (int iDim =0; iDim < spaceDim; iDim++)
		xOld[iDim] = xParticle[iDim];
}

void Particle::UpdateParticleVelocities(Real* uP, Real* omP) {

	for (int iDim = 0; iDim < spaceDim; iDim++) {
		uParticle[iDim] = uP[iDim];
		omegaParticle[iDim] = omP[iDim];
	}

}

void Particle::UpdateStencil(Real* xBounds, int *Nf, Real dx) {

	Real xMin[spaceDim];
	int iLocal[spaceDim];
	int iRadius;
	for (int iDim = 0; iDim < spaceDim; iDim++)
		xMin[iDim] = xBounds[2 * iDim];


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


void Particle::GetParticlePositions(Real* xPos) {

	for (int iDim = 0; iDim < spaceDim; iDim++)
		xPos[iDim] = xParticle[iDim];
}

void Particle::SetInputVariables(std::vector<Real>& inputData) {

	if (inputData.size() > nInputExtraVariables) {
		nInputExtraVariables = inputData.size();
		inputVariables.reserve(nInputExtraVariables);
	}

	for (int iDir = 0; iDir < nInputExtraVariables; iDir++)
		inputVariables.at(iDir) = inputData.at(iDir);
}

void Particle::GetInputVariables(std::vector<Real>& inputData) {

	for (int iDir = 0; iDir < nInputExtraVariables; iDir++)
		inputData.at(iDir) = inputVariables.at(iDir);
}


void Particle::GetOutputVariables(std::vector<Real>& outputData) {

	for (int iDir = 0; iDir < nOutputExtraVariables; iDir++) {
		outputData.at(iDir) = outputVariables.at(iDir);
	}

}

void Particle::SetOutputVariables(std::vector<Real>& outputData) {

	if (outputData.size() > nOutputExtraVariables) {
		nOutputExtraVariables = outputData.size();
		outputVariables.reserve(nOutputExtraVariables);
	}

	for (int iDir = 0; iDir < nOutputExtraVariables; iDir++)
		outputVariables.at(iDir) = outputData.at(iDir);

}
