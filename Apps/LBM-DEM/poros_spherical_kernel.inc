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
 * @brief   Kernels for poros_spherical model
 * @author  C. Tsigginos
 * @details Kernels for creating and updating particle
 * 			mapping model based on a spherical control volume
 */


#ifndef POROS_SPHERICAL_KERNEL_INC_
#define POROS_SPHERICAL_KERNEL_INC_

#include "type.h"
#include "model.h"
#include "flowfield_host_device.h"
#include "model_host_device.h"
#include "poros_spherical.h"
#include "poros_spherical_host_device.h"

void KerSolidFracSphere(ACC<int>& id, ACC<Real>& sfp, ACC<Real>& vp,
		ACC<Real>& xAvg, const ACC<Real>& xf, const Real* xPos, const Real* Radius,
		const int* idP, const Real* vPart, const Real* omPart, const Real* dx,
		const int* nelem, const int* spacedim) {

	const int SPACEDIM = (*spacedim);
	const int NELEM = (*nelem);
	int no = 0;
	int noIndex;
	Real vf[SPACEDIM], xfl[SPACEDIM], xavg[SPACEDIM];
	Real sfEst;
	Real dx1, dy1, dz1;
	xfl[0] = xf(0, 0, 0, 0);
	xfl[1] = xf(1, 0, 0, 0);
	xfl[2] = xf(2, 0, 0, 0);

	Real Ravg = (*dx) * sqrt(1/PI);

	sfEst = CalculateSolidFractionSpheres(xfl, Ravg, xPos, (*Radius), xavg); //TODO

	//find noTot;
	for (int iPar = 0; iPar < NELEM; iPar++)
		no += (id(iPar, 0, 0, 0)>-1);


	if (sfEst > 0.0) {
		no+=1;
		noIndex = no - 1;
		id(noIndex, 0, 0, 0) = (*idP);
		sfp(noIndex, 0, 0, 0) = sfEst;
		xAvg(noIndex * SPACEDIM, 0, 0, 0) = xavg[0];
		xAvg(noIndex * SPACEDIM + 1, 0, 0, 0) = xavg[1];
		xAvg(noIndex * SPACEDIM + 2, 0, 0, 0 ) = xavg[2];

		dx1 = xavg[0] - xPos[0];
		dy1 = xavg[1] - xPos[1];
		dz1 = xavg[2] - xPos[2];
		vp(SPACEDIM * noIndex, 0, 0, 0) = vPart[0] + omPart[1] * dz1 - omPart[2] * dy1;
		vp(SPACEDIM * noIndex + 1, 0, 0, 0) =  vPart[1] + omPart[2] * dx1 - omPart[0] * dz1;
		vp(SPACEDIM * noIndex + 2, 0, 0, 0) = vPart[2] + omPart[0] * dy1 - omPart[1] * dx1;

	}


}

void KerPorousSphericalUpdate(ACC<Real>& vP, const ACC<int>& id,
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







#endif /* APPS_LBM_DEM_NOP_POROS_SPHERICAL_KERNEL_INC_ */
