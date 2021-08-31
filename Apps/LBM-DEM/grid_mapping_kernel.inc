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
 * @brief   Kernels for mapping particles via the subcell method
 * @author  C. Tsigginos
 * @details The kernels are used in conjuction with PSM or Prati schemes
 */

#ifndef GRID_MAPPING_KERNEL_INC_
#define GRID_MAPPING_KERNEL_INC_

#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "model.h"
#include "flowfield_host_device.h"
#include "model_host_device.h"
#include "grid_mapping.h"
#include "grid_mapping_host_device.h"

#ifdef OPS_3D
void KerMappingSphericalGrid3D(ACC<int>& idp, ACC<Real>& sfp, ACC<Real>& uPar,
							 ACC<Real>& xAvg,const ACC<Real>& xF, const Real* xPart,
							 const Real* Radius, const int* idParticle,
							 const Real* vParticle, const Real* omParticle,
							 const Real* dx, const int* Ngrid, const int* noelem,
							 const int* spacedim) {


	const int SPACEDIM = (*spacedim);
	const int NOELEM = (*noelem);
	const int noGrid = (*Ngrid);

	Real xfl[SPACEDIM];
	Real xe, ye, ze;
	Real dd;
	Real dx1, dy1, dz1;
	Real xavg[SPACEDIM];
	int ip, jp, kp;
	int status, nEdges, noIndex;
	int itot = 0;
	int sx[26] = { 1,  1,  1,  1, -1, -1, -1, -1,  1, -1,  0,  0,  0,  0,  1, -1,
				1, -1,  0,  0,  0,  0,  1,  1, -1, -1};
	int sy[26] = { 1,  1, -1, -1,  1,  1, -1, -1,  0,  0,  1, -1,
				0,  0,  1,  1,-1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
	int sz[26] = { 1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,
				0,  0, 0,  0,  1,  1, -1, -1,  1, -1,  1, -1};

	nEdges = 26;
	dx1 = (*dx);
	int no = 0;

	xfl[0] = xF(0, 0, 0, 0);
	xfl[1] = xF(1, 0, 0, 0);
	xfl[2] = xF(2, 0, 0, 0);

	no  = 0;
	for (int iPar = 0; iPar < NOELEM; iPar++)
		no += (idp(iPar, 0, 0, 0)>-1);


	for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
		xe = xfl[0] + 0.5 * (dx1) * static_cast<Real>(sx[edgeIndex]);
		ye = xfl[1] + 0.5 * (dx1) * static_cast<Real>(sy[edgeIndex]);
		ze = xfl[2] + 0.5 * (dx1) * static_cast<Real>(sz[edgeIndex]);

		dd = (xe - xPart[0]) * (xe - xPart[0]) + (ye - xPart[1]) * (ye - xPart[1]) + (ze - xPart[2]) * (ze - xPart[2]);

		if (dd <= (*Radius) * (*Radius))
			itot += 1;
	}


	if (itot == nEdges) {
		no += 1;
		noIndex = no - 1;
		idp(noIndex, 0, 0, 0) = (*idParticle);
		sfp(noIndex, 0, 0, 0) = 1.0;
		dx1 = xfl[0] - xPart[0];
		dy1 = xfl[1] - xPart[1];
		dz1 = xfl[2] - xPart[2];

		uPar(SPACEDIM * noIndex, 0, 0, 0) = vParticle[0] + omParticle[1] * dz1
				- omParticle[2] * dy1;
		uPar(SPACEDIM * noIndex +1, 0, 0, 0) = vParticle[1] + omParticle[2] * dx1
				- omParticle[0] * dz1;
		uPar(SPACEDIM * noIndex +2, 0, 0, 0) = vParticle[2] + omParticle[0] * dy1
				- omParticle[1] * dx1;

		xAvg(SPACEDIM * noIndex, 0, 0, 0) = xfl[0];
		xAvg(SPACEDIM * noIndex + 1, 0, 0, 0) = xfl[1];
		xAvg(SPACEDIM * noIndex + 2, 0, 0, 0) = xfl[2];

	}
	else if (itot > 0) {
		no += 1;
		noIndex = no - 1;
		idp(noIndex, 0, 0, 0) = (*idParticle);
		sfp(noIndex, 0, 0, 0) = ComputeCircularSolidFractionGrid(xfl, xPart, (*Radius),
				xavg, (*dx), noGrid);

		xAvg(SPACEDIM * noIndex, 0, 0, 0) = xavg[0];
		xAvg(SPACEDIM * noIndex + 1, 0, 0, 0) = xavg[1];
		xAvg(SPACEDIM * noIndex + 2, 0, 0, 0) = xavg[2];

		dx1 = xavg[0] - xPart[0];
		dy1 = xavg[1] - xPart[1];
		dz1 = xavg[2] - xPart[2];

		uPar(SPACEDIM * noIndex, 0, 0, 0) = vParticle[0] + omParticle[1] * dz1
				- omParticle[2] * dy1;
		uPar(SPACEDIM * noIndex +1, 0, 0, 0) = vParticle[1] + omParticle[2] * dx1
				- omParticle[0] * dz1;
		uPar(SPACEDIM * noIndex +2, 0, 0, 0) = vParticle[2] + omParticle[0] * dy1
				- omParticle[1] * dx1;


	}

}

void KerSolidVelocityUpdateGridSphere3D(ACC<Real>& vP, const ACC<int>& id,
		const ACC<Real>& xAvg, const int* idParticle,
		const Real* xPos, const Real* radPart, const Real* velP, const Real* omP,
		const Real* dx, const Real* spacedim, const Real* noelem) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelem);
	Real sfPart;
	Real xf[SPACEDIM], vEst[SPACEDIM];
	Real dx1, dy1, dz1;
	for (int iPart = 0; iPart < NELEM; iPart++) {
		if (id(iPart, 0, 0, 0) == (*idParticle)) {
			xf[0] = xAvg(iPart * SPACEDIM, 0, 0, 0);
			xf[1] = xAvg(iPart * SPACEDIM + 1, 0, 0, 0);
			xf[2] = xAvg(iPart * SPACEDIM + 2, 0, 0, 0);



			dx1 = xf[0] - xPos[0];
			dy1 = xf[1] - xPos[1];
			dz1 = xf[2] - xPos[2];

			vP(SPACEDIM * iPart, 0, 0, 0) = velP[0] + omP[1] * dz1 - omP[2] * dy1;
			vP(SPACEDIM * iPart + 1, 0, 0, 0) =  velP[1] + omP[2] * dx1 - omP[0] * dz1;
			vP(SPACEDIM * iPart + 2, 0, 0, 0) = velP[2] + omP[0] * dy1 - omP[1] * dx1;


			break;
		}
	}
}
#endif

#ifdef OPS_2D
void KerMappingSphericalGrid(ACC<int>& idp, ACC<Real>& sfp, ACC<Real>& uPar,
							 ACC<Real>& xAvg,const ACC<Real>& xF, const Real* xPart,
							 const Real* Radius, const int* idParticle,
							 const Real* vParticle, const Real* omParticle,
							 const Real* dx, const int* Ngrid, const int* noelem,
							 const int* spacedim) {


	const int SPACEDIM = (*spacedim);
	const int NOELEM = (*noelem);
	const int noGrid = (*Ngrid);

	Real xfl[SPACEDIM];
	Real xe, ye;
	Real dd;
	Real dx1, dy1, dz1;
	Real xavg[SPACEDIM];
	int ip, jp;
	int status, nEdges, noIndex;
	int itot = 0;
	int sx[8] = { -1,  1, 1, -1,  0, 1, 0, -1};
	int sy[8] = { -1, -1, 1,  1, -1, 0, 1,  0};


	nEdges = 8;
	dx1 = (*dx);
	int no = 0;

	xfl[0] = xF(0, 0, 0);
	xfl[1] = xF(1, 0, 0);
	xfl[2] = xF(2, 0, 0);

	no  = 0;
	for (int iPar = 0; iPar < NOELEM; iPar++)
		no += (idp(iPar, 0, 0)>-1);


	for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
		xe = xfl[0] + 0.5 * (dx1) * static_cast<Real>(sx[edgeIndex]);
		ye = xfl[1] + 0.5 * (dx1) * static_cast<Real>(sy[edgeIndex]);

		dd = (xe - xPart[0]) * (xe - xPart[0]) + (ye - xPart[1]) * (ye - xPart[1]);

		if (dd <= (*Radius) * (*Radius))
			itot += 1;
	}


	if (itot == nEdges) {
		no += 1;
		noIndex = no - 1;
		idp(noIndex, 0, 0) = (*idParticle);
		sfp(noIndex, 0, 0) = 1.0;
		dx1 = xfl[0] - xPart[0];
		dy1 = xfl[1] - xPart[1];


		uPar(SPACEDIM * noIndex, 0, 0) = vParticle[0]  - omParticle[2] * dy1;
		uPar(SPACEDIM * noIndex +1, 0, 0) = vParticle[1] + omParticle[2] * dx1;


		xAvg(SPACEDIM * noIndex, 0, 0) = xfl[0];
		xAvg(SPACEDIM * noIndex + 1, 0, 0) = xfl[1];
		xAvg(SPACEDIM * noIndex + 2, 0, 0) = xfl[2];

	}
	else if (itot > 0) {
		no += 1;
		noIndex = no - 1;
		idp(noIndex, 0, 0) = (*idParticle);
		sfp(noIndex, 0, 0) = ComputeCircularSolidFractionGrid2D(xfl, xPart, (*Radius),
				xavg, (*dx), noGrid);

		xAvg(SPACEDIM * noIndex, 0, 0) = xavg[0];
		xAvg(SPACEDIM * noIndex + 1, 0, 0) = xavg[1];
		xAvg(SPACEDIM * noIndex + 2, 0, 0) = xavg[2];

		dx1 = xavg[0] - xPart[0];
		dy1 = xavg[1] - xPart[1];


		uPar(SPACEDIM * noIndex, 0, 0) = vParticle[0] - omParticle[2] * dy1;
		uPar(SPACEDIM * noIndex +1, 0, 0) = vParticle[1] + omParticle[2] * dx1;


	}

}

void KerSolidVelocityUpdateGridSphere(ACC<Real>& vP, const ACC<int>& id,
		const ACC<Real>& xAvg, const int* idParticle,
		const Real* xPos, const Real* radPart, const Real* velP, const Real* omP,
		const Real* dx, const Real* spacedim, const Real* noelem) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelem);
	Real sfPart;
	Real xf[SPACEDIM], vEst[SPACEDIM];
	Real dx1, dy1;
	for (int iPart = 0; iPart < NELEM; iPart++) {
		if (id(iPart, 0, 0) == (*idParticle)) {
			xf[0] = xAvg(iPart * SPACEDIM, 0, 0);
			xf[1] = xAvg(iPart * SPACEDIM + 1, 0, 0);




			dx1 = xf[0] - xPos[0];
			dy1 = xf[1] - xPos[1];

			vP(SPACEDIM * iPart, 0, 0) = velP[0] - (*omP) * dy1;
			vP(SPACEDIM * iPart + 1, 0, 0) =  velP[1] + (*omP) * dx1;



			break;
		}
	}
}
#endif
#endif