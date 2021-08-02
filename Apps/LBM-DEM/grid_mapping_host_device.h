#ifndef GRID_MAPPING_HOST_DEVICE_H
#define GRID_MAPPING_HOST_DEVICE_H

#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

#include <cmath>
#include "type.h"

static inline OPS_FUN_PREFIX Real ComputeCircularSolidFractionGrid(const Real* xfl,
		const Real* xPart, const Real Rp, Real* xavg, const Real dx, const int noGrid) {

	int Nx = noGrid + 2;
	int Ny = noGrid + 2;
	int Nz = noGrid + 2;
	int inter = 0;

	Real dxp, dyp, dzp, dd;
	Real xmin, ymin, zmin;
	Real xtrial, ytrial, ztrial;
	Real dxc = dx/static_cast<Real>(Nx-1);

	xmin = xfl[0] - 0.5 * dx;
	ymin = xfl[1] - 0.5 * dx;
	zmin = xfl[2] - 0.5 * dx;

	for (int iDir = 0; iDir  < 3; iDir++)
		xavg[iDir] = 0.0;

	for (int ip = 0; ip < Nx; ip++) {
		xtrial = xmin + static_cast<Real>(ip) * dxc;
		dxp = (xtrial - xPart[0]) * (xtrial - xPart[0]);
		for (int jp = 0; jp < Ny; jp++) {
			ytrial = ymin + static_cast<Real>(jp) * dxc;
			dyp = (ytrial - xPart[1]) * (ytrial - xPart[1]);
			for (int kp = 0; kp < Nz; kp++) {
				ztrial = zmin + static_cast<Real>(kp) * dxc;
				dzp = (ztrial - xPart[2]) * (ztrial - xPart[2]);
				dd = dxp + dyp + dzp;
				if (dd < Rp * Rp ) {
					inter+= 1;
					xavg[0] += xtrial;
					xavg[1] += ytrial;
					xavg[2] += ztrial;
				}
			}
		}
	}

	for (int iDir = 0; iDir < 3; iDir++)
		xavg[iDir] /= static_cast<Real>(inter);

	return static_cast<Real>(inter) / static_cast<Real>(Nx * Ny * Nz);

}

#endif
