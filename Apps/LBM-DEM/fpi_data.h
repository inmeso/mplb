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
 * @brief   Class for storing fluid-particle interaction ops and user defined variables
 * @author  C. Tsigginos
 */

#ifndef NOP_FPI_DATA_H_
#define NOP_FPI_DATA_H_

#include "flowfield.h"
#include "type.h"
#include "dem_data.h"
#include <memory>
#include <vector>
#include "block.h"
#include "model.h"

class FpiData {

	private:
		Component compo;
		int spaceDim;
		int noElem;
		std::vector<Real> userVariables;
		bool ownedPostVelocity;
		bool ownedForceModel;
		bool ownedCollisionModel;
		bool ownedPreCollision;
		bool ownedPostStreaming;
		bool ownedDragForce;
		bool ownedInitialize;
		bool thermalFlag;
		FluidParticleModel model;

		std::vector<RealField> mappingRealVariableList;
		std::vector<IntField> mappingIntVariableList;

		std::vector<int> mappingIntSize;
		std::vector<int> mappingRealSize;

		std::vector<std::string> nameRealVariables;
		std::vector<std::string> nameIntVariables;


	public:
		FpiData(Component& compoUser, int spacdim, FluidParticleModel model, int noelem,
				std::vector<std::string> nameOfRealVars, std::vector<int> sizeRealVars,
				std::vector<std::string> nameOfIntVars,  std::vector<int> sizeIntVars,
				std::vector<bool> flagParameters, std::vector<Real> variables,
				SizeType timeStep = 0);

		Component& FluidComponent() {return compo;};
		bool OwnedPostVelocity() {return ownedPostVelocity;};
		bool OwnedForceModel() {return ownedForceModel;};
		bool OwnedPreCollision() {return ownedPreCollision;};
		bool OwnedCollisionModel() {return ownedCollisionModel;};
		bool OwnedPostStreaming() {return ownedPostStreaming;};
		bool OwnedDragForce() {return ownedDragForce;};
		bool OwnedInitialize() {return ownedInitialize;};
		bool ThermalFlag() {return thermalFlag;};
		int GetNElem() {return noElem;};
		int GetSpaceDim() {return spaceDim;};
		FluidParticleModel GetModel() {return model;};
		ops_dat&  GetRealField(int variable, int blockIndex);
		ops_dat&  GetIntField(int variable, int blockIndex);
		std::string ReturnRealVariablesName(int variable) {return nameRealVariables.at(variable);};
		std::string ReturnIntVariablesName(int variable) {return nameIntVariables.at(variable); };
		int  SizeAtIntType(int variable) { return mappingIntSize.at(variable);};
		int  SizeAtRealType(int variable) { return mappingRealSize.at(variable);};
		Real GetRealVariableAt(int elem) { return userVariables.at(elem);};
		void WriteToHdf5(const std::string& caseName, const SizeType timeStep);
};

extern std::map<int, std::shared_ptr<FpiData>> fpiModels;

#endif /* APPS_LBM_DEM_NOP_FPI_DATA_H_ */
