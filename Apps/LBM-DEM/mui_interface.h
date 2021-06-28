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

/*! @brief  Head files for the base Mui class
 * @author C. Tsigginos
 */

#ifndef MUI_INTERFACE_BASE_H_
#define MUI_INTERFACE_BASE_H_

#include <limits>
#include <vector>
#include <map>
#include <string>
#include "block_particles.h"
#include "mui.h"
#include "type.h"
#include "dem_data.h"
class MuiInterface {

	public:
		MuiInterface(ListBlockParticles* blockParticles,
				std::vector<std::string> input, std::vector<std::string> particleShape,
				std::vector<std::string> output, Real skin = 1.0, SizeType timeStep = 0);
	    ~MuiInterface();
		void SetDomains(SizeType maxIteration);
		void UpdateDomains(SizeType steps);
		void ExtractData(SizeType currentStep, SizeType &firstStep, SizeType& maxStep,
				Real& alpha, int* flags, ParticleShapeDiscriptor& particleShape);
		void ExtractParticles(SizeType timeStep);
		void SendParticles(SizeType timesStep);
		void ForgetData(SizeType timeStep);
		void UpdateRegionFirstLast(SizeType firstStep, SizeType endStep);


	private:
		mui::uniface3d* interface;			//mui object
		SizeType maxStep;					//maximum step that the assinged region is ok
		Real Rmax;							//maximum particle size for defining region
		void DefineProcBox(Real* xmin, Real* xmax);
		ListBlockParticles*  blockPointer;  //point to list block particles
		int spaceDim;						//size of spatial space
		SizeType startStep;					//Start point of simulation
		std::vector<std::string> inputData;
		std::vector<std::string> outputData;
		std::vector<std::string> particleShapeData;
		int inputDataSize;
		int outputDataSize;
		int inputParticleSize;
		void SendExtraParameters();
		void ReceiveExtraParameters();
		void ReceiveParticleDetails();
};



#endif /* APPS_LBM_DEM_MUI_INTERFACE_BASE_H_ */
