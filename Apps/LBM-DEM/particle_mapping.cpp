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
 * @brief   Functions for defining and handling particle mapping models
 * @author  C. Tsigginos
 * @details Function define particle mapping per component, mapping particles to
 * 			LBM grid, update grid lists
 */

#include "particle_mapping.h"
#include "flowfield.h"
#include <memory>
#include "block_particles.h"
#include "grid_mapping.h"


void DefineParticleMappingModels(std::vector<ParticleMappingModel> mappingModels,
								 std::vector<int> compoId, std::vector<int> copyFrom, ParticleShapeDiscriptor particleShape,
								 SizeType timeStep) {

	int noComp = mappingModels.size();
	int noIDComp = compoId.size();
	int spaceDim = SpaceDim();

	std::map<int, Component> components = g_Components();
	int numComponents = ComponentNum();

	if (noIDComp != noComp) {
		ops_printf("ERROR: The size of input data is inconsistent\n"
				   "noID = %d noComp = %d", noComp, noIDComp);
		exit(EXIT_FAILURE);
	}
	//Checking if copyModel is given after the first
	for (int elem = 0; elem < noComp; elem++) {
		if (mappingModels.at(elem) == copyData) {
			if (elem < copyFrom.at(elem)) {
				int ax = copyFrom.at(elem);
				int ay = -1;
				for (int copyElem = 0; copyElem < noComp; copyElem++) {
					if (compoId.at(copyElem) == ax) {
						ay = copyElem;
						break;
						if (mappingModels.at(copyElem) == copyData) {
							ops_printf("ERROR: mapping model cannot copy from another copy model\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				if (elem < ay) {
					mappingModels.at(elem) = mappingModels.at(ay);
					mappingModels.at(ay) = copyData;
					copyFrom.at(elem) = -1;
					copyFrom.at(ay) = compoId.at(elem);

					int compoTmp = compoId.at(elem);
					compoId.at(elem) = compoId.at(ay);
					compoId.at(ay) = compoTmp;
				}

			}

		}
	}

	//Creating maps
	for (int elem = 0; elem < noComp; elem++) {
		int iComponent = compoId.at(elem);
		switch (mappingModels[elem]) {
			case (noMapping):  {
				int noElem = 2;
				std::vector<std::string> nameOfRealVars{};
				std::vector<std::string> nameOfIntVars{};
				std::vector<int> sizeReals{};
				std::vector<int> sizeInt{};

				std::shared_ptr<MappingParticles>  model(new MappingParticles(spaceDim, noElem, particleShape,
								mappingModels.at(elem),  nameOfRealVars, sizeReals,
								nameOfIntVars,  sizeInt));
				particleMappingModels.insert(std::make_pair(iComponent, model));

			} break;
			case (sphericalMapping): {
				int noElem = 2;
				std::string real1{"sfP"};
				std::string real2{"vP"};
				std::string real3{"xAvg"};
				std::string intS{"idP"};

				real1 += std::to_string(iComponent);
				real2 += std::to_string(iComponent);
				real3 += std::to_string(iComponent);
				intS += std::to_string(iComponent);
				std::vector<std::string> nameOfRealVars{real1, real2, real3};
				std::vector<std::string> nameOfIntVars{intS};
				std::vector<int> sizeReals{noElem, noElem * spaceDim, noElem * spaceDim};
				std::vector<int> sizeInt{noElem};


				std::shared_ptr<MappingParticles> model(new MappingParticles(spaceDim, noElem, particleShape,
									mappingModels.at(elem),  nameOfRealVars, sizeReals,
										nameOfIntVars,  sizeInt));


				particleMappingModels.insert(std::make_pair(iComponent, model));
			} break;
			case (gridSpherical) : {
				int noElem = 2;
				std::string real1{"sfP"};
				std::string real2{"vP"};
				std::string real3{"xAvg"};
				std::string intS{"idP"};

				real1 += std::to_string(iComponent);
				real2 += std::to_string(iComponent);
				real3 += std::to_string(iComponent);
				intS += std::to_string(iComponent);
				std::vector<std::string> nameOfRealVars{real1, real2, real3};
				std::vector<std::string> nameOfIntVars{intS};
				std::vector<int> sizeReals{noElem, noElem * spaceDim, noElem * spaceDim};
				std::vector<int> sizeInt{noElem};


				std::shared_ptr<MappingParticles> model(new MappingParticles(spaceDim, noElem, particleShape,
									mappingModels.at(elem),  nameOfRealVars, sizeReals,
											nameOfIntVars,  sizeInt));


				particleMappingModels.insert(std::make_pair(iComponent, model));


			}  break;
			case gridQuadratic : {
				int noElem = 8;
				std::string real1{"sfP"};
				std::string real2{"vP"};
				std::string real3{"xAvg"};
				std::string intS{"idP"};

				real1 += std::to_string(iComponent);
				real2 += std::to_string(iComponent);
				real3 += std::to_string(iComponent);
				intS += std::to_string(iComponent);
				std::vector<std::string> nameOfRealVars{real1, real2, real3};
				std::vector<std::string> nameOfIntVars{intS};
				std::vector<int> sizeReals{noElem, noElem * spaceDim, noElem * spaceDim};
				std::vector<int> sizeInt{noElem};


				std::shared_ptr<MappingParticles> model(new MappingParticles(spaceDim, noElem, particleShape,
														mappingModels.at(elem),  nameOfRealVars, sizeReals,
														nameOfIntVars,  sizeInt));


				particleMappingModels.insert(std::make_pair(iComponent, model));



			} break;
			case gridMesh : {
				int noElem = 8;
				std::string real1{"sfP"};
				std::string real2{"vP"};
				std::string real3{"xAvg"};
				std::string intS{"idP"};

				real1 += std::to_string(iComponent);
				real2 += std::to_string(iComponent);
				real3 += std::to_string(iComponent);
				intS += std::to_string(iComponent);
				std::vector<std::string> nameOfRealVars{real1, real2, real3};
				std::vector<std::string> nameOfIntVars{intS};
				std::vector<int> sizeReals{noElem, noElem * spaceDim, noElem * spaceDim};
				std::vector<int> sizeInt{noElem};



				std::shared_ptr<MappingParticles> model(new MappingParticles(spaceDim, noElem, particleShape,
														mappingModels.at(elem),  nameOfRealVars, sizeReals,
														nameOfIntVars,  sizeInt));


				particleMappingModels.insert(std::make_pair(iComponent, model));
			} break;
			case copyData: {
				particleMappingModels.insert(std::make_pair(iComponent, particleMappingModels.at(copyFrom.at(elem))));
			} break;
			default: {
				int noElem = 2;
					std::vector<std::string> nameOfRealVars{"sfP", "vP", "xAvg"};
					std::vector<std::string> nameOfIntVars{"idP"};
					std::vector<int> sizeReals{noElem * spaceDim, noElem * spaceDim, noElem * spaceDim};
					std::vector<int> sizeInt{noElem * spaceDim};
					ParticleMappingModel  tmpModel{sphericalMapping};

					std::shared_ptr<MappingParticles> model(new MappingParticles(spaceDim, noElem, particleShape,
									tmpModel,  nameOfRealVars, sizeReals,
												nameOfIntVars,  sizeInt));


					particleMappingModels.insert(std::make_pair(iComponent, model));
			}
		}


	}

}

std::map<int, std::shared_ptr<MappingParticles>>& GetMappingModels() {

	return particleMappingModels;
}

void InitializeMappingLists() {

	int component;
	for (auto& mappingModel : particleMappingModels) {
		component = mappingModel.first;
		switch(mappingModel.second->MappingModel()) {
			case noMapping:
				break;
			case copyData:
				break;
			default: { //The default is the porous media type model
				InitializePorousLists(mappingModel.second, component);

			}
		}

	}
}

void MappingParticlesToLBMGrid() {

	int component;

	for (auto& mappingModel : particleMappingModels) {
		component = mappingModel.first;
		switch(mappingModel.second->MappingModel()) {
			case noMapping:
				break;
			case copyData:
				break;
			case gridQuadratic:
				ops_printf("ERROR: This model is not imported yet\n");
				break;
			case gridMesh:
				ops_printf("ERROR: This model is not imported yet\n");
				break;
			case gridSpherical:
				ParticleProjectionSphereGrid(mappingModel.second, component);
				break;
			case sphericalMapping:
				ParticleProjectionSphere(mappingModel.second, component);
				break;
			default:
				ParticleProjectionSphere(mappingModel.second, component);

		}
	}

}

void WriteParticleMappingToHDF5(SizeType currentStep) {

	for (auto& mappingParticle : GetMappingModels())
		mappingParticle.second->WriteToHdf5(CaseName(), currentStep);
}

void UpdateParticlesToLBMGrid() {

	int component;

	for (auto &mappingModel : particleMappingModels) {
		component = mappingModel.first;
		switch (mappingModel.second->MappingModel()) {
			case noMapping:
				break;
			case copyData:
				break;
			case gridQuadratic:
				ops_printf("ERROR: This model is not imported yet\n");
				break;
			case gridMesh:
				ops_printf("ERROR: This model is not imported yet\n");
				break;
			case gridSpherical:
				UpdateProjectectionSphereGrid(mappingModel.second, component);
				break;
			case sphericalMapping:
				UpdateProjectionSphere(mappingModel.second, component);
				break;
			default:
				UpdateProjectionSphere(mappingModel.second, component);
		}
	}

}
