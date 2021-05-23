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
 * @brief   Class for handling data exchange between LBM-DEM code
 * @author  Chrysovalantis Tsigginos
 * @details Default class for handling data exchange between MP-LBM and DEM codes. The current version
 * handles only momentum exchange data
 */

#ifndef MUI_INTERFACE_H_
#define MUI_INTERFACE_H_


#include "type.h"
#include "flowfield.h"
#include "box_handling.h"
#include "dem_particles.h"
#include "periodic.h"
#include "mui.h"
#include <vector>

using namespace std;

class muiInterface {
	public:
		muiInterface(Real skin = 1.0); //constructor of 3D class
		~muiInterface();
		void setDomains(int maxIteration); //transferring LBM domains
		void updateDomains(int steps);  //Update domain on the fly
		void extractData(long int timestep,long int& maxStep, Real& alpha, int* Flags);
		void extractDataPeriodic(Real* xper, Real* xcutOff);
		void updateParticles(int timeStep, vector<Real> &xTemp, vector<Real> &yTemp, vector<Real> &zTemp, vector<Real> &radTmp,
						   vector<Real> &uTemp, vector<Real> &vTemp, vector<Real> &wTemp,
						   vector<Real> &omXTemp, vector<Real> &omYTemp, vector<Real> &omZTemp);

		void sendParticles(int timestep);
		void forgetData(int timestep);

	private:
		void DefineProcBox(Real* xMin, Real* xmax);
		mui::uniface3d* interface;
		int maxStep;
		Real Rmax;

};




#endif /* APPS_LBM_DEM_CLASS_INTERFACE_MUI_INTERFACE_H_ */
