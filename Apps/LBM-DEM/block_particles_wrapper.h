/*
 * block_particles_wrapper.h
 *
 *  Created on: May 14, 2021
 *      Author: jpd38567
 */

#ifndef BLOCK_PARTICLES_WRAPPER_H_
#define BLOCK_PARTICLES_WRAPPER_H_

#include "block.h"
#include "scheme.h"
#include "flowfield.h"
#include "type.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif

void updateParticleList(Block* currentBlock, Real* xb, int* range, int spacedim);

void KerCarBound(const ACC<Real>& xf, Real* xb);




#endif /* APPS_LBM_DEM_BLOCK_PARTICLES_WRAPPER_H_ */
