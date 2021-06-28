/*
 * poros_spherical.inc
 *
 *  Created on: May 28, 2021
 *      Author: jpd38567
 */

#ifndef POROS_SPHERICAL_INC_
#define POROS_SPHERICAL_INC_

#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

#include "poros_spherical.h"
#include "poros_spherical_host_device.h"

#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"


void PorosSpherical::KerSolidFracSphere(ACC<int>& id, ACC<Real>& sfp, ACC<Real>& vp,
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

void PorosSpherical::KerSolidVelocitySphere(ACC<Real>& vP, const ACC<int>& id,
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

void PorosSpherical::KerInitializePorousSpherical(ACC<Real>&sfP, ACC<Real>& vP,
		ACC<Real>& xAvg, ACC<int>& id, const int* spacedim, const int* noelem) {

	const int SPACEDIM = (*spacedim);
	const int NELEM = (*noelem);

	for (int iPart = 0; iPart < NELEM; iPart++) {
		sfP(iPart, 0, 0, 0) = 0.0;
		for (int iDir = 0; iDir < SPACEDIM; iDir++) {
			id(iPart * SPACEDIM + iDir, 0, 0, 0) = -1;
			xAvg(iPart * SPACEDIM + iDir, 0, 0, 0) = -1.0;
		}
	}

}


void PorosSpherical::KerPrintPorousData(const ACC<Real>&sfP,const ACC<Real>& vP,
		const ACC<Real>& xAvg, const ACC<int>& id, const ACC<Real>& xf,
		const int* spacedim, const int* noelem) {

	const int SPACEDIM = (*spacedim);
	const int NELEM = (*noelem);

	for (int iPart = 0; iPart < NELEM; iPart++) {
		printf("[%f %f %f]: sfP = %f, idP = %d vAvf=[%f %f %f] xAvg =[%f %f %f]\n",
				xf(0, 0, 0, 0), xf(1, 0, 0, 0), xf(2, 0, 0, 0),
				sfP(iPart, 0, 0, 0), id(iPart, 0, 0, 0), vP(iPart*SPACEDIM, 0, 0, 0),
				vP(iPart*SPACEDIM + 1, 0, 0, 0), vP(iPart*SPACEDIM + 2, 0, 0, 0),
				xAvg(iPart*SPACEDIM, 0, 0, 0),
				xAvg(iPart*SPACEDIM + 1, 0, 0, 0), xAvg(iPart*SPACEDIM + 2, 0, 0, 0));

		}



}


#endif /* APPS_LBM_DEM_POROS_SPHERICAL_INC_ */
