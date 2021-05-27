/*
 * Block_particles_wrapper.cpp
 *
 *  Created on: May 14, 2021
 *      Author: jpd38567
 */

#include "block_particles_wrapper.h"
#include "block_particles_wrapper.inc"

void updateParticleList(Block* currentBlock, Real* xb, int* iterRng, int spaceDim) {

	int blockIndex{currentBlock->ID()};

	ops_par_loop(KerCarBound, "KerCarBound", currentBlock->Get(), spaceDim, iterRng,
				 ops_arg_dat(g_CoordinateXYZ()[blockIndex], spaceDim,
							 LOCALSTENCIL, "double", OPS_READ),
				 ops_arg_gbl(xb, spaceDim, "double", OPS_READ));//,

}
