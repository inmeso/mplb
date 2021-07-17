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

/*! @brief   Abstract storage of particle mapping data
 * @author  C. Tsigginos
 * @details Implementing functions related to create particle mapping
 * 			storage function.
 */

#ifndef MAPPING_PARTICLES_H_
#define MAPPING_PARTICLES_H_

#include "field.h"
#include "type.h"
#include "flowfield.h"
#include "dem_data.h"
#include <vector>
#include <string>
#include <map>
#include <memory>



class MappingParticles {

private:


	ParticleShapeDiscriptor particleDiscriptor;
	ParticleMappingModel mappingModel;
	int spaceDim;
	int requiresCopy;
	int noElem;
	std::vector<RealField> mappingRealVariableList;
	std::vector<IntField> mappingIntVariableList;
	std::vector<int> mappingIntSize;
	std::vector<int> mappingRealSize;
	std::vector<std::string> nameRealVariables;
	std::vector<std::string> nameIntVariables;

public:
	MappingParticles(int spacedim, int noElements, ParticleShapeDiscriptor particleShape,
					ParticleMappingModel mappingMod,
					std::vector<std::string> nameOfRealVars, std::vector<int> sizeRealVars,
					std::vector<std::string> nameOfIntVars,  std::vector<int> sizeIntVars,
					SizeType timestep = 0);

	~MappingParticles() { };
	ParticleShapeDiscriptor  ParticleShape() {return particleDiscriptor;};
	ParticleMappingModel     MappingModel() {return mappingModel;};
	int 	  GetSpaceDim() {return spaceDim; };
	ops_dat&  GetRealField(int variable, int blockIndex);
	std::string ReturnRealVariablesName(int variable) {return nameRealVariables.at(variable);};
	std::string ReturnIntVariablesName(int variable) {return nameIntVariables.at(variable); };
	ops_dat&  GetIntField(int variable, int blockIndex);
	int  SizeAtIntType(int variable) { return mappingIntSize.at(variable);};
	int  SizeAtRealType(int variable) { return mappingRealSize.at(variable);};
	void WriteToHdf5(const std::string& caseName, const SizeType timeStep);
	int NumberOfElements() {return noElem;};
};

using mappingModelList = std::map<int, std::shared_ptr<MappingParticles>>;
extern mappingModelList particleMappingModels;



#endif /* APPS_LBM_DEM_NOP_MAPPING_PARTICLES_H_ */
