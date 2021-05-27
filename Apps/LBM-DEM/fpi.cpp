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

/*! @brief Define main functions
 * @author C. Tsigginos
 * @details Declaring the base definitions for the fluid-particle interaction model
 */

#include "fpi.h"

void DefineInteractionModel(std::vector<FSIType> FluidParticleType,
		std::vector<SolidFracType>  SolFracType, vector<int> fsiCompoId,
		double gamma, SizeType timeStep) {

	int noComp = FluidParticleType.size();
	int noSol = SolFracType.size();
	int noFSICompo = fsiCompoId.size();


	if (noComp != NUMCOMPONENTS) {
		ops_printf("Error: Number of FSI components inconsistent with actual components\n");
		assert(noComp != NUMCOMPONENTS);
	}

	if (noComp != NUMCOMPONENTS) {
		opr_printf("Error: Inconsistent number of solid fraction components with actual components\n");
		assert(noComp != NUMCOMPONENTS);
	}

	if (noFSICompo != NUMCOMPONENTS) {
		ops_printf("Error: Number of fsiCompoId differs from actual number of components\n");
		assert(noFSICompo != NUMCOMPONENTS)
	}

	for (int iComp = 0; iComp < noComp; iComp++) {
		if (FluidParticleType[iComp] == Model_None) {
			/* TODO Add general case */


		}
		else if (FluidParticleType[iComp] == Model_PSM) {
			/* TODO add PSM case */
			/* TODO Assign solid fraction if needed
		}
		else if (FluidParticleType[iComp] == Model_Prati) {
			/* TODO add Model Prati*/
			/* TODO assign solid fraction */
		}
		else {
			ops_printf("Entered here because new model is not defined properly\n");
			ops_printf("The empty model will be added automatically\n");
		}

	}






}
