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

/*! @brief Base class for mapping particles into the LBM grid,
 *  @author C. Tsigginos
 **/



#ifndef PARTICLE_TO_GRID_BASE_H_
#define PARTICLE_TO_GRID_BASE_H_

#include "field.h"
#include "type.h"
#include <vector>
#include <string.h>

class ParticleToGridBase {

	protected:
		enum ParticleType { ParticleNone = 0, ParticleSpherical = 1,
			ParticleSuperQuadratic= 2, ParticleMesh = 3};

		ParticleType particleDiscriptor;
		int spaceDim;
		int requiresCopy;
		int noElem;
		std::vector<RealField> mappingRealVariableList;
		std::vector<IntField> mappingIntVariableList;

	public:

		ParticleToGridBase(int particleshape, int spacedim);
		virtual ~ParticleToGridBase() { };
		virtual int  particleShape() {return 0;}
		virtual void DefineVariables(int noElem, SizeType timestep = 0) { };
		virtual void ParticleProjection() { }
		virtual void UpdateProjection() { }
		virtual void WriteToHdf5(const std::string& caseName, const SizeType timeStep);
		virtual void InitializeVariables() { }
		virtual RealField& GetRealFieldVariable(int index);
		virtual IntField& GetIntFieldVariable(int index);
};



#endif
