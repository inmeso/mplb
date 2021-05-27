/*
 * block_particles_wrapper.inc
 *
 *  Created on: May 14, 2021
 *      Author: jpd38567
 */

#ifndef BLOCK_PARTICLES_WRAPPER_INC_
#define BLOCK_PARTICLES_WRAPPER_INC_

#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "block_particles_wrapper.h"

void KerCarBound(const ACC<Real>& xf, Real* xb) {//,int* spacedim) {

	//int size = (* spacedim);
	int size;
	size = 3;
#ifdef OPS_2D
	size = 2;
#endif

	for (int idx = 0; idx < size; idx++) {
#ifdef OPS_2D
		xb[idx] = xf(idx, 0, 0);
#endif

#ifdef OPS_3D
		xb[idx] = xf(idx, 0 , 0, 0);
#endif

	}



}






#endif /* APPS_LBM_DEM_BLOCK_PARTICLES_WRAPPER_INC_ */
