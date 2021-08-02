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
 * @details: Members functions of FPiData class. store data for FSI models
 */

#include "fpi_data.h"

std::map<int, std::shared_ptr<FpiData>> fpiModels;

FpiData::FpiData(Component& compoUser, int spacdim, FluidParticleModel modelU, int noelem,
		std::vector<std::string> nameOfRealVars, std::vector<int> sizeRealVars,
		std::vector<std::string> nameOfIntVars,  std::vector<int> sizeIntVars,
		std::vector<bool> flagParameters, std::vector<Real> variables,
		SizeType timeStep) : compo{compoUser} {

	model = modelU;
	noElem = noelem;
	spaceDim = spacdim;
	model = modelU;
	//Define bools
	if (flagParameters.size() != 8) {
		ops_printf("ERROR: Number of user flags is larger than the defined number of flags\n");
		exit(EXIT_FAILURE);
	}

	ownedPostVelocity = flagParameters.at(0);
	ownedForceModel = flagParameters.at(1);
	ownedPreCollision = flagParameters.at(2);
	ownedCollisionModel = flagParameters.at(3);
	ownedPostStreaming = flagParameters.at(4);
	ownedDragForce = flagParameters.at(5);
	ownedInitialize = flagParameters.at(6);
	thermalFlag = flagParameters.at(7);
	//Define user required variables
	for (const auto& variable : variables)
		userVariables.push_back(variable);

	//Define ops_model parameters
	if (sizeRealVars.size() != nameOfRealVars.size()) {
		ops_printf("ERROR: Inconsistent sizes for Real ops_dat variables for"
				   " generating OPS_dat variables for fluid-particle interaction.\n");
		exit(EXIT_FAILURE);
	}

	for (int iVar = 0; iVar < nameOfRealVars.size(); iVar++) {
		RealField variable{nameOfRealVars.at(iVar), sizeRealVars.at(iVar)};

		mappingRealVariableList.push_back(variable);
		nameRealVariables.push_back(nameOfRealVars.at(iVar));
		mappingRealSize.push_back(sizeRealVars.at(iVar));
	}

	for (auto& mappingVariable : mappingRealVariableList) {
		if (timeStep == 0)
			mappingVariable.CreateFieldFromScratch(g_Block());
		else
			mappingVariable.CreateFieldFromFile(CaseName(), g_Block(), timeStep);
	}

	if (sizeIntVars.size() != nameOfIntVars.size()) {
			ops_printf("ERROR: Inconsistent sizes for int ops_dat variables"
					" for fluid-particle interaction models .\n");
			exit(EXIT_FAILURE);
	}

	for (int iVar = 0; iVar < nameOfIntVars.size(); iVar++) {
		IntField variable{nameOfIntVars.at(iVar),  sizeIntVars.at(iVar)};
		mappingIntVariableList.push_back(variable);
		nameIntVariables.push_back(nameOfIntVars.at(iVar));
		mappingIntSize.push_back(sizeIntVars.at(iVar));

	}

	for (auto& mappingVariable : mappingIntVariableList) {
		if (timeStep == 0)
			mappingVariable.CreateFieldFromScratch(g_Block());
		else
			mappingVariable.CreateFieldFromFile(CaseName(), g_Block(), timeStep);
	}

}

ops_dat& FpiData::GetRealField(int variable, int blockIndex) {

	return mappingRealVariableList.at(variable).at(blockIndex);
}

ops_dat& FpiData::GetIntField(int variable, int blockIndex) {

	return mappingIntVariableList.at(variable).at(blockIndex);
}

void FpiData::WriteToHdf5(const std::string& caseName, SizeType timeStep) {

	for (const auto& variableReal : mappingRealVariableList)
		variableReal.WriteToHDF5(caseName, timeStep);

	for (const auto& variableInt : mappingIntVariableList)
		variableInt.WriteToHDF5(caseName, timeStep);


}


