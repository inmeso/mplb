
/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]s
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

/*! @brief   Define functions and variables for DEM-LBM simulations
 * @author  C. Tsigginos
 */

#ifndef DEM_HANDLE_H_
#define DEM_HANDLE_H_

#include "block_particles.h"
#include "flowfield.h"
#include <map>
#include "block.h"
#include "fpi.h"
#include <vector>
#include "model.h"
#include "block_particle_helper.h"
#include "dem_data.h"
#include <string>




void SetupBlockParticles(InteractionData data, SizeType maxStep);
void DefineGlobalBlockBox();
void CreateMuiInterface(InteractionData data, std::vector<FSIType> models, Real skin);
void StreamCollisionFSI3D(int flag);
void IterateFSI(InteractionData data, int savingFlag);
void IterateFSI(Real convergenceRate,const SizeType checkPointPeriod,const SizeType maxIters);

void SetupDEMLBM(InteractionData& data);
void ReadParticleDataSpherical();

void SetDemLbMParams(InteractionData* data, bool flag, bool muiOn, Real convergeRate,
		SizeType checkperiod, SizeType timeStep, SizeType checkPeriodSteady,
		SizeType maximumIterations, std::string particleType);

void SetupParticleBoxes(InteractionData data);
void SetupRestartSimulation(SizeType timestep, SizeType endStep);
#endif /* APPS_LBM_DEM_DEM_HANDLE_H_ */
