/*
 * mui_wrapper.cpp
 *
 *  Created on: Jun 7, 2021
 *      Author: jpd38567
 */



#include "mui_wrapper.h"

MuiInterface* interface;

void DefineMuiInterface(ListBlockParticles* blockParticles,
		std::vector<std::string>  input,  std::vector<std::string> particleShape,
		std::vector<std::string> output, Real skin) {

	interface = new MuiInterface(blockParticles, input, particleShape, output, skin);

}

void SetMuiDomains(SizeType maxStep) {

	interface->SetDomains(maxStep);

}

void ExtractSimulationData(SizeType timeStep, SizeType& start, SizeType& maxIter, Real& alpha, int* flags,
		ParticleShapeDiscriptor& particleShape) {

	interface->ExtractData(timeStep, start, maxIter, alpha, flags, particleShape);

}

void ExtractParticleData(SizeType timeStep) {

	interface->ExtractParticles(timeStep);
}

void SendParticleData(SizeType timeStep) {

	interface->SendParticles(timeStep);
}
