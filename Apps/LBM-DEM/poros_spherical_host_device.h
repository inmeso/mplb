/*
 * poros_spherical_device.h
 *
 *  Created on: May 29, 2021
 *      Author: jpd38567
 */

#ifndef POROS_SPHERICAL_HOST_DEVICE_H_
#define POROS_SPHERICAL_HOST_DEVICE_H_

#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

#include <cmath>
#include "poros_spherical.h"

inline OPS_FUN_PREFIX Real PorosSpherical::CalculateSolidFractionSpheres(const Real* xfl,
		const Real Ravg, const Real* xPos, const Real Rp, Real* xAvg) {


	Real d, x, sf, h1, h2, V1, V2, Vol, norm;
	Real z1, z2, zcm, n[3];
	d = (xPos[0] - xfl[0]) * (xPos[0] - xfl[0]) +
			(xPos[1] - xfl[1]) * (xPos[1] - xfl[1]) +
			(xPos[2] - xfl[2]) * (xPos[2] - xfl[2]);
	d=sqrt(d);
	if (d < Ravg + Rp) {
		if (d<=fabs(Rp-Ravg)) {
			for (int iDim = 0; iDim < 3; iDim++)
				xAvg[iDim] = xfl[iDim];
			return 1.0;
		}
		x = (d * d - Rp * Rp + Ravg * Ravg)/(2.0 * d);
		h1 = Ravg - x;
		h2 = Rp - d + x;
		V1 = (1.0/3.0) * PI * h1 * h1 * (3.0 * Ravg - h1);
		V2 = (1.0/3.0) * PI * h2 * h2 * (3.0 * Rp - h2);
		Vol = (4.0/3.0) * PI * Ravg * Ravg * Ravg;

		z1 = 3.0 * (2.0 * Ravg - h1) * (2.0 * Ravg - h1) / (4.0 * (3.0 * Ravg - h1));
		z2 = 3.0 * (2.0 * Rp - h2) * (2.0 * Rp - h2) / (4.0 * (3.0 * Rp - h2));
		//Transform z2 into z1
		z2 = d - z2;
		zcm = (z2 * V2 + z1 * V1) / (V1 + V2);

		//Transofmration of zcm into global coordinate system
		n[0] = xPos[0] - xfl[0];
		n[1] = xPos[1] - xfl[1];
		n[2] = xPos[2] - xfl[2];
		norm = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		norm = sqrt(norm);
		for (int iDir = 0; iDir < 3; iDir++)
			n[iDir] /= norm;


		for (int iDir = 0; iDir < 3; iDir++)
			xAvg[iDir] = xfl[iDir] + zcm * n[iDir];

		return (V1 + V2) / Vol;
	}
	xAvg[0] = -1.0;
	xAvg[1] = -1.0;
	xAvg[2] = -1.0;
	return 0.0;
}



#endif
