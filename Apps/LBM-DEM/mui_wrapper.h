/*
 * mui_wrapper.h
 *
 *  Created on: Jun 7, 2021
 *      Author: jpd38567
 */

#ifndef MUI_WRAPPER_H_
#define MUI_WRAPPER_H_

#include "mui_interface.h"
#include <vector>
#include <string>
#include "block_particles.h"
#include "type.h"

void DefineMuiInterface(ListBlockParticles* blockParticles,
		std::vector<std::string>  input, std::vector<std::string> particleShape,
		std::vector<std::string> output, Real skin = 1.0);

void SetMuiDomains(SizeType maxStep);
void ExtractSimulationData(SizeType timeStep,SizeType& firstStep, SizeType& maxIter,
		Real& alpha, int* flags, ParticleShapeDiscriptor& particleShape);
void ExtractParticleData(SizeType timeStep);
void SendParticleData(SizeType timeStep);
void UpdateRegions(SizeType startStep, SizeType lastStep);
#endif /* APPS_LBM_DEM_MUI_WRAPPER_H_ */
