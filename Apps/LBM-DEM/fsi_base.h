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

/*! @brief Base class for fluid-particle interaction modeling
 *  @author C. Tsigginos
 **/

#ifndef FSI_BASE_H_
#define FSI_BASE_H_

#include "flowfield.h"
#include <string>
#include "block.h"
#include "type.h"
#include "model.h"


enum FSIType { Model_Prati = 2, Model_PSM = 1, Model_None = 0};
enum SolFracType {Mode_None = 0, Mode_Spherical = 1, Mode_Grid = 2, Mode_Copy = 3};


class FsiBase {

	protected:
		Component compo;  			//Component associated with this FSI model
		Real gamma; 				//User-defined parameter
		SolFracType porosModel;		//Model for mapping particles to the grid
		int spaceDim;
		Real* force;				//Base force
		int forceFlag;				//Flag for force update
		int isThermalModel = 0; 	//Thermal flag
	public:
		bool collisionOwned;  //(true) model invoked own collision

		FsiBase(Component& compoUser, int spacedim, Real* forceUser, bool owned = false,
				SolFracType porosModel = Mode_None, Real gammaUser = 0.0);
		virtual ~FsiBase();
		virtual void ModelCollision() {} //inputs to be determined
		virtual void PostVelocityCalculation() { } //inputs to be determined
		virtual void PreCollision() { } //inputs to be de determined
		virtual void PostStreaming() {} //inputs to be determined
		virtual void MappingFunction(bool flag) {}
		virtual void InitializeVariables() {}
		virtual void CalculateDragForce() { } // inputs to be determined
		virtual void SetupSimulation() { } //inputs to be determined
		virtual void RestartSimulation() { } // inputs to be determined
		virtual void WriteToHdf5(const std::string& caseName, const SizeType timeStep) {}
		virtual void DefineConstants() {} //To be added
		virtual void DefineVariables(SizeType timestep = 0) {} //Base functions for defining model required variables
		virtual void ObtainID(int* velId, int* loop, Real& tauCompo,
				CollisionType& collisModel, int& idComponent, int& rhoId,int &thId);
		virtual int GetThermalFlag() {return isThermalModel;}


};



#endif /* APPS_LBM_DEM_FLPART_INTERACTION_BASE_H_ */
