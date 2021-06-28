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

/*! @brief   Base calls for mapping particles in the LBM grid
 * @author C. Tsigginos
 * @details Base consturctor and output files
 */

#include "particle_to_grid_base.h"

ParticleToGridBase::ParticleToGridBase(int particleshape, int spacedim) {

	spaceDim = spacedim;

	particleDiscriptor = (ParticleType) particleshape;
	ops_printf("Particle type passed is %d\n", particleshape);
	ops_printf("Particle type assigned is %d\n", (ParticleType) particleDiscriptor);
	noElem = 1;

	requiresCopy = false;

}

void ParticleToGridBase::WriteToHdf5(const std::string& caseName, const SizeType timeStep) {

	for (auto& element : mappingRealVariableList )
		element.WriteToHDF5(caseName, timeStep);

	for (auto& element : mappingIntVariableList)
		element.WriteToHDF5(caseName, timeStep);

}

RealField& ParticleToGridBase::GetRealFieldVariable(int index) {

	return mappingRealVariableList.at(index);
}

IntField& ParticleToGridBase::GetIntFieldVariable(int index) {

	return mappingIntVariableList.at(index);
}


