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


#include "mui_interface.h"
#include <limits>
#include <vector>
using namespace std;

muiInterface::muiInterface(Real skin) {

	interface = new mui::uniface3d("mpi://LBM/ifs");
	maxStep = 0;
	Rmax = skin;
}

muiInterface::~muiInterface() {

	delete interface;
}

void muiInterface::setDomains(int maxIteration) {

	maxStep = 2 * maxIteration;

	//Find minimum-maximum proc size
	Real xmin[SPACEDIM], xmax[SPACEDIM];

	DefineProcBox(xmin, xmax);


	mui::geometry::box3d send_recv_region({xmin[0], xmin[1], xmin[2]}, {xmax[0], xmax[1], xmax[2]});


	interface->announce_send_span(0, maxStep, send_recv_region);
	interface->announce_recv_span(0, maxStep, send_recv_region);

	printf("Rank %d announces region [%f %f] x[%f %f] x[%f %f]\n",ops_get_proc(), xmin[0], xmax[0], xmin[1], xmax[1], xmin[2], xmax[2]);
}

void muiInterface::updateDomains(int step) {

	if (step == maxStep) {
		maxStep = 2 * step;

		Real xmin[SPACEDIM], xmax[SPACEDIM];


		DefineProcBox(xmin, xmax);


		mui::geometry::box3d send_recv_region({xmin[0], xmin[1], xmin[2]}, {xmax[0], xmax[1], xmax[2]});


		interface->announce_send_span(step, maxStep, send_recv_region);
		interface->announce_recv_span(step, maxStep, send_recv_region);

		ops_printf("OPS-LBM: Update region");
		printf("Rank %d announces region [%f %f] x[%f %f] x[%f %f]\n",ops_get_proc(), xmin[0], xmax[0], xmin[1], xmax[1], xmin[2], xmax[2]);

	}

}

void muiInterface::DefineProcBox(Real* xMin, Real* xMax) {

	xMin[0] = std::numeric_limits<Real>::max();
	xMin[1] = xMin[0];
	xMin[2] = xMin[0];

	xMax[0] = -1.0 * std::numeric_limits<Real>::max();
	xMax[1] = xMax[0];
	xMax[2] = xMax[0];

	for (int blockIndex = 0 ; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		if  (xBoundLocal[6 * blockIndex] < xMin[0])
			xMin[0] = xBoundLocal[6 * blockIndex];

		if (xBoundLocal[6 * blockIndex + 2] < xMin[1])
			xMin[1] = xBoundLocal[6 * blockIndex + 2];

		if (xBoundLocal[6 * blockIndex + 4] < xMin[2])
			xMin[2] = xBoundLocal[6 * blockIndex + 4];

		if (xBoundLocal[6 * blockIndex + 1] > xMax[0])
			xMax[0] = xBoundLocal[6 * blockIndex + 1];

		if (xBoundLocal[6 * blockIndex + 3] > xMax[1])
			xMax[1] = xBoundLocal[6 * blockIndex + 3];

		if (xBoundLocal[6 * blockIndex + 5] > xMax[2])
			xMax[2] = xBoundLocal[6 * blockIndex + 5];

	}

	for (int iDim = 0; iDim < SPACEDIM; iDim++) {
		xMin[iDim] -= Rmax;
		xMax[iDim] += Rmax;
	}

}

void muiInterface::extractData(long int timestep,long int& maxStep, Real& alpha, int* Flags) {

	std::vector<mui::point3d> posMUI;
	mui::chrono_sampler_exact3d time_sampler;
	posMUI = interface->fetch_points<double, mui::chrono_sampler_exact3d>("radius", timestep, time_sampler);

	auto maxIter1 = interface->fetch<long int>("Nsteps");
	auto alpha1 = interface->fetch<Real>("alpha1");

	maxStep = maxIter1;
	alpha = alpha1;

	auto perFlagX = interface->fetch<int>("xPerFlag");
	auto perFlagY = interface->fetch<int>("yPerFlag");
	auto perFlagZ = interface->fetch<int>("zPerFlag");

	Flags[0] = perFlagX;
	Flags[1] = perFlagY;
	Flags[2] = perFlagZ;

}

void muiInterface::extractDataPeriodic(Real* xper, Real* xcutOff) {



	auto xper1 = interface->fetch<Real>("xper");
	auto yper1 = interface->fetch<Real>("yper");
	auto zper1 = interface->fetch<Real>("zper");


	auto xCutOff = interface->fetch<Real>("xCutOff");
	auto yCutOff = interface->fetch<Real>("yCutOff");
	auto zCutOff = interface->fetch<Real>("zCutoff");

	xcutOff[0] = xCutOff;
	xcutOff[1] = yCutOff;
	xcutOff[2] = zCutOff;

	xper[0] = xper1;
	xper[1] = yper1;
	xper[2] = zper1;


}

void muiInterface::updateParticles(int timeStep, vector<Real> &xTemp, vector<Real> &yTemp, vector<Real> &zTemp, vector<Real> &radTmp,
						   vector<Real> &uTemp, vector<Real> &vTemp, vector<Real> &wTemp,
						   vector<Real> &omXTemp, vector<Real> &omYTemp, vector<Real> &omZTemp) {

	std::vector<mui::point3d> posMUI;
	mui::chrono_sampler_exact3d time_sampler;
	mui::sampler_exact3d<Real> s1;

	posMUI = interface->fetch_points<Real, mui::chrono_sampler_exact3d>("radius", timeStep, time_sampler);

	radTmp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("radius", timeStep, time_sampler);

	uTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("u", timeStep, time_sampler);

	vTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("v", timeStep, time_sampler);

	wTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("w", timeStep, time_sampler);

	omXTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("ox", timeStep, time_sampler);

	omYTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("oy", timeStep, time_sampler);

	omZTemp = interface->fetch_values<double, mui::chrono_sampler_exact3d>("oz", timeStep, time_sampler);


	//Copy particles to correct location
	for (auto iVec = 0; iVec < posMUI.size(); ++iVec) {
		auto pointTmp = posMUI[iVec];//check  access


		Real xtemp = pointTmp[0];
		Real ytemp = pointTmp[1];
		Real ztemp = pointTmp[2];

		xTemp.push_back(xtemp);
		yTemp.push_back(ytemp);
		zTemp.push_back(ztemp);
    }

	interface->forget(timeStep);
}

void muiInterface::forgetData(int timestep) {

	interface->forget(timestep);
}

void muiInterface::sendParticles(int timestep) {

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
	#ifdef OPS_MPI
			sub_block_list sb = OPS_sub_block_list[blockIndex];
			if (!sb->owned) continue;
	#endif
			for (int iPar = 0; iPar < Nparticles[blockIndex]; iPar++) {
				interface->push("Fdx",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar]);
				interface->push("Fdy",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar + 1]);
				interface->push("Fdz",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar + 2]);
				interface->push("Mdx",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar + 3]);
				interface->push("Mdy",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar + 4]);
				interface->push("Mdz",{xp[blockIndex][iPar],yp[blockIndex][iPar],zp[blockIndex][iPar]},FDrag[blockIndex][6 * iPar + 5]);

			//printf("Iter %d, Rank %d: Fd [%12.9e %12.9e %12.9e]\n",t, ops_get_proc(), FDrag[blockIndex][6 * iPar],FDrag[blockIndex][6 * iPar +1], FDrag[blockIndex][6 * iPar +2] );
			}
	}

	interface->commit(timestep);
}

