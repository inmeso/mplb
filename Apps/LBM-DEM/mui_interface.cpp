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
 * @brief   Class for data exchange between MPLB and DEM solver
 * @author  C. Tsigginos
 * @details: Functions for intersolver communications via MUI
 */


#include "mui_interface.h"
#include <map>
#include <vector>

MuiInterface::MuiInterface( ListBlockParticles* blockParticles,
		std::vector<std::string> input, std::vector<std::string> particleShape,
		std::vector<std::string> output, Real skin, SizeType timeStep)
{

	interface = new mui::uniface3d("mpi://LBM/ifs");
	maxStep = 0;
	Rmax = skin;
	startStep = timeStep;
	blockPointer = blockParticles;
	spaceDim = 2;
#ifdef OPS_3D
	spaceDim = 3;
#endif

	inputData = input;
	outputData = output;
	particleShapeData = particleShape;

	inputDataSize = inputData.size();
	outputDataSize = outputData.size();
	inputParticleSize = particleShapeData.size();

}

MuiInterface::~MuiInterface() {

	delete interface;
}

void MuiInterface::SetDomains(SizeType maxIteration) {

	maxStep = 2 * maxIteration + startStep;
	Real xmin[spaceDim], xmax[spaceDim];
	long int firstStep = static_cast<long int>(startStep);
	long int lastStep = static_cast<long int>(maxStep);

	DefineProcBox(xmin, xmax);

	mui::geometry::box3d send_recv_region({xmin[0], xmin[1], xmin[2]},
			{xmax[0], xmax[1], xmax[2]});

	interface->announce_send_span(firstStep, lastStep, send_recv_region);
	interface->announce_recv_span(firstStep, lastStep, send_recv_region);

}

void MuiInterface::DefineProcBox(Real* xMin, Real* xMax) {

	Real xMinTmp[spaceDim], xMaxTmp[spaceDim];

	xMin[0] = std::numeric_limits<Real>::max();
	xMin[1] = xMin[0];
	xMin[2] = xMin[0];

	xMax[0] = -1.0 * xMin[0];
	xMax[1] = xMax[0];
	xMax[2] = xMax[0];


	for (auto particleBlock  = blockPointer->begin(); particleBlock != blockPointer->end();
			++particleBlock) {
		particleBlock->second.ExtractBound(xMinTmp, xMaxTmp);

		for (int iDir = 0; iDir < spaceDim; iDir++) {
			if (xMinTmp[iDir] < xMin[iDir])
				xMin[iDir] = xMinTmp[iDir];
			if (xMaxTmp[iDir] > xMax[iDir])
				xMax[iDir] = xMaxTmp[iDir];
		}
	}

#ifdef CPU
#if DebugLevel >= 2
	printf("Send MUI region of rank %d: [(%e %e %e), (%e %e %e])]\n", ops_get_proc(),
			xMin[0], xMin[1], xMin[2], xMax[0], xMax[1], xMax[2]);
#endif
#endif
}

void MuiInterface::ExtractData(SizeType currentStep, SizeType &firstStep, SizeType& maxIter,
		Real& alpha, int* flags, ParticleShapeDiscriptor& particleShape) {

	std::vector<mui::point3d> posMUI;
	mui::chrono_sampler_exact3d time_sampler;
	posMUI = interface->fetch_points<double, mui::chrono_sampler_exact3d>("radius", currentStep, time_sampler);

	auto maxIter1 = interface->fetch<long int>("Nsteps");
	auto alpha1 = interface->fetch<Real>("alpha1");
	auto particleCase = interface->fetch<int>("particleShape");
	auto firstStep1 = interface->fetch<long int>("firstStep");
	firstStep = static_cast<SizeType>(firstStep1);
	maxIter = static_cast<SizeType>(maxIter1) - firstStep;
	alpha = alpha1;
	particleShape = (ParticleShapeDiscriptor) particleCase;
}

void MuiInterface::ExtractParticles(SizeType timestep) {


	Real uTmp[spaceDim], oTmp[spaceDim];
	Real xTemp[spaceDim], radTmp;
	Real partTemp;
	std::vector<Real> particleShape;
	particleShape.reserve(inputParticleSize);
	std::vector<Real> inputVariables;
	inputVariables.reserve(inputDataSize);
	long int timeframe = static_cast<long int>(timestep);
	std::vector<mui::point3d> posMUI;
	int iDir;

	mui::chrono_sampler_exact3d timeSampler;
	mui::sampler_exact3d<Real> s1;

	posMUI = interface->fetch_points<Real, mui::chrono_sampler_exact3d>("radius",
			timeframe, timeSampler);

	for (auto particleBlock  = blockPointer->begin(); particleBlock != blockPointer->end(); ++particleBlock)
		particleBlock->second.ClearParticles();

	//Copy particles to correct location
	for (auto iVec = 0; iVec < posMUI.size(); ++iVec) {
		auto pointTmp = posMUI[iVec];//check  access

		xTemp[0] = pointTmp[0];
		xTemp[1] = pointTmp[1];
		xTemp[2] = pointTmp[2];

		//fetch other parameters
		radTmp = interface->fetch("radius", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
		uTmp[0] = interface->fetch("u", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
		uTmp[1] = interface->fetch("v", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
		uTmp[2] = interface->fetch("w", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);

		oTmp[0] = interface->fetch("ox", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
		oTmp[1] = interface->fetch("oy", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
		oTmp[2] = interface->fetch("oz", {xTemp[0], xTemp[1], xTemp[2]}, timeframe,
				s1, timeSampler);
#ifdef CPU
#if DebugLevel >= 2
		printf("Rank %d: [%f %f %f] R= %f u =[%f %f %f] omega=[%f %f %f]\n",
				ops_get_proc(), xTemp[0], xTemp[1], xTemp[2], radTmp, uTmp[0],
				uTmp[1], uTmp[2], oTmp[0], oTmp[1], oTmp[2]);
#endif
#endif
		if (inputParticleSize > 0) {
			ops_printf("Additional particle data are passed to MPLB via MUI\n");

			iDir = 0;
			for (auto& input : particleShapeData) {
				iDir += 1;
				auto dataTmp = interface->fetch(input, {xTemp[0], xTemp[1], xTemp[2]}, timeframe, s1, timeSampler);
				particleShape.at(iDir) = dataTmp;
			}
		}

		if (inputDataSize > 0) {
			ops_printf("Additional input are passed to MPLB via MUI\n");
			iDir = 0;
			for (auto& input: inputData) {
				iDir += 1;
				auto dataTmp = interface->fetch(input, {xTemp[0], xTemp[1], xTemp[2]}, timeframe, s1, timeSampler);
				inputVariables.at(iDir) = dataTmp;
			}
		}


		for (auto particleBlock  = blockPointer->begin(); particleBlock != blockPointer->end(); ++particleBlock) {

			int idParticle = particleBlock->second.InsertParticle(xTemp, radTmp, particleShape,
					uTmp, oTmp, inputVariables);
			//TODO ADD the additional parameters.
			if (inputParticleSize > 0 && idParticle > -1)
				particleBlock->second.GetAdditionalInputVariables(inputVariables, idParticle);
		}


	}

	interface->forget(timeframe);

}

void MuiInterface::SendParticles(SizeType timeStep) {

	if (timeStep > maxStep)
		SetDomains(10 * maxStep );
	long int timeframe = static_cast<long int>(timeStep);
	Real Fd[spaceDim], Td[spaceDim], xPos[spaceDim];
	std::vector<Real> output;
	output.reserve(outputDataSize);
	int Nparticles;
	for (auto particleBlock =  blockPointer->begin(); particleBlock != blockPointer->end(); ++particleBlock) {
		Nparticles = particleBlock->second.NParticles;
		for (int iParticle = 0; iParticle < Nparticles; iParticle++) {
			particleBlock->second.ExtractDragForce(Fd, Td, iParticle);
			particleBlock->second.ExtractPositions(xPos, iParticle);

			printf("Rank %d at timestep %d: Sends for particle [%f %f %f], Fd=[%12.9e %12.9e %12.9e] Td=[%12.9e %12.9e %12.9e\n",
					ops_get_proc(), timeStep, xPos[0], xPos[1], xPos[2], Fd[0], Fd[1], Fd[2], Td[0], Td[1], Td[2]);
			interface->push("Fdx", {xPos[0], xPos[1], xPos[2]}, Fd[0]);
			interface->push("Fdy", {xPos[0], xPos[1], xPos[2]}, Fd[1]);
			interface->push("Fdz", {xPos[0], xPos[1], xPos[2]}, Fd[2]);
			interface->push("Mdx", {xPos[0], xPos[1], xPos[2]}, Td[0]);
			interface->push("Mdy", {xPos[0], xPos[1], xPos[2]}, Td[1]);
			interface->push("Mdz", {xPos[0], xPos[1], xPos[2]}, Td[2]);

			if (outputDataSize > 0) {
				particleBlock->second.GetAdditionalOutputVariables(output, iParticle);
				for (int iDir = 0; iDir < outputDataSize; iDir++)
					interface->push(outputData[iDir], {xPos[0], xPos[1], xPos[2]}, output[iDir]);
			}
		}

	}

	interface->commit(timeframe);
}

void MuiInterface::UpdateRegionFirstLast(SizeType firstStep, SizeType endStep) {

	int a1  = 0;
	int a2 = 0;
	if (firstStep > startStep) {
		a1 = 1;
		startStep = firstStep;
	}
	if (endStep > maxStep) {
		a2 = 1;
		maxStep = endStep;

	}

	if  ((a1 == 1) || (a2 == 1)) {
		SizeType iters= maxStep - firstStep;
		SetDomains(iters);
	}

}
