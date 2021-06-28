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
 * @brief   Wrapper functions for the data exchange operations
 * @author  C. Tsigginos
 * @details: Wrapper functions for intersolver communication
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

void UpdateRegions(SizeType startStep, SizeType lastStep) {

	interface->UpdateRegionFirstLast(startStep, lastStep);
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
