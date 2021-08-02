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
 */
#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "type.h"
#include "particle_shape.h"
#include "dem_data.h"
#include <vector>
#include <string>
#include <memory>
class Particle {

	public:
		Real xParticle[3]; //Particle Center of Mass
		std::shared_ptr<ParticleShape> particleShape;
		Real FDrag[3];
		Real TDrag[3];
		Real xOld[3];
		Real uParticle[3];
		Real omegaParticle[3];
		int stenList[6];
		ParticleShapeDiscriptor particleShapeType;
		int spaceDim;


		Particle(int Dimension, Real radius, Real* xp, std::vector<Real> Shape,
				ParticleShapeDiscriptor particleImport, int inputVariablesSize,
				int outputVariableSize);
		~Particle();
		void InitializeDrag();
		void AddDrag(Real* Fp, Real* Tp);
		void EvaluateDrag(Real dt);
		void GetParticlePositions(Real* xPos);
		std::string GetParticleShape();
		void UpdateParticleLocation(Real* xp);
		void UpdateParticleShape(Real Rmax, std::vector<Real> Shape);
		void UpdateParticleVelocities(Real* uPart, Real* omegaPart);
		void UpdateExtraVariable(std::vector<Real> inputdata) { };
		void UpdateOutputVariable(std::vector<Real>& outputData) { };
		void PushDrag(Real* Fd, Real* Td);
		void UpdateStencil(Real* xBounds, int *Nf, Real dx);
		void UpdateOldParticlePositions();
		void SetInputVariables(std::vector<Real>& inputData);
		void GetInputVariables(std::vector<Real>& inputData);
		void GetOutputVariables(std::vector<Real>& outputData);
		void SetOutputVariables(std::vector<Real>& outputData);

	private:
		int nInputExtraVariables;
		int nOutputExtraVariables;
		std::vector<Real> inputVariables;  //User defined input variables
		std::vector<Real> outputVariables; //User defined output variables
};






#endif /* APPS_LBM_DEM_PARTICLE_H_ */
