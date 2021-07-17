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
 * @brief   Functions for running fluid-particle simulations
 * @author  C. Tsigginos
*/

#ifndef FPI_FUNCTIONS_H_
#define FPI_FUNCTIONS_H_

#include "fpi_data.h"
#include <memory>
#include <vector>
#include <map>
#include <string>
#include "type.h"
#include "block_particles.h"
#include "dem_data.h"
#include "particle_mapping.h"
#include "mapping_particles.h"
#include "fsi_models.h"


void DefineInteractionModels(std::vector<FluidParticleModel> particleModels,
						     std::vector<int> compoId, std::vector<Real> variables,
							 SizeType timeStep=0);

void DefineNone(int iComponent, int spaceDim, FluidParticleModel model);
void DefinePsm(int iComponent, int spaceDim, FluidParticleModel model,
				std::vector<Real> variables, SizeType timeStep = 0);

void FSIVelocityFunctions(std::shared_ptr<FpiData>& fpiPtr);
void FsiForceFunction(std::shared_ptr<FpiData>& fpiPtr);
void FsiPreCollisionFunction(std::shared_ptr<FpiData>& fpiPtr);
void FsiCollisionFunction(std::shared_ptr<FpiData>& fpiPtr);
void FsiPostCollisionFunction(std::shared_ptr<FpiData>& fpiPtr);
void FsiCalculateDragForce(std::shared_ptr<FpiData>& fpiPtr);
void ObtainData(std::shared_ptr<FpiData>& fpiPtr, int* idVel, int* loop, Real& tauCompo,
		CollisionType& collisModel, int& idComponent,int& rhoId,int &thId);
void FsiInitializeFunction(std::shared_ptr<FpiData>& fpiPtr);

#endif /* APPS_LBM_DEM_NOP_FPI_FUNCTIONS_H_ */
