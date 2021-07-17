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
 * @brief   Functions for defining and handling particle mapping models
 * @author  C. Tsigginos
 */


#ifndef PARTICLE_MAPPING_H_
#define PARTICLE_MAPPING_H_

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "dem_data.h"
#include "mapping_particles.h"
#include "model.h"
#include "flowfield.h"
#include "mapping_models.h"

void DefineParticleMappingModels(std::vector<ParticleMappingModel> particleMappingModels,
		 std::vector<int> compoId, std::vector<int> copyFrom, ParticleShapeDiscriptor particleShape,
		 SizeType timeStep);

void InitializeMappingLists();
void MappingParticlesToLBMGrid();
void UpdateParticlesToLBMGrid();
void WriteParticleMappingToHDF5(SizeType timestep);
void InitializePorousLists(std::shared_ptr<MappingParticles>& mappingPtr,
						   int component);
std::map<int, std::shared_ptr<MappingParticles>>& GetMappingModels();
#endif /* APPS_LBM_DEM_NOP_PARTICLE_MAPPING_H_ */
