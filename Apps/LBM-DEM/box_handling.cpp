/*
 * box_handling.cpp
 *
 *  Created on: Sep 18, 2020
 *      Author: jpd38567
 */
#include "box_handling.h"
#include <limits>
#include <stdlib.h>
#include <stdio.h>


Real*   x_bounds;
Real*    xBoundLocal;
int* Nf;


void  allocateBlockVariables(int blockNum) {

  x_bounds= new Real[2 * SPACEDIM * blockNum ];
  xBoundLocal = new Real [2 * SPACEDIM * blockNum];
  Nf = new int[2 * SPACEDIM * blockNum];
}

void deleteBlockVariables() {

	delete[] x_bounds;
	delete[] xBoundLocal;
	delete[] Nf;

}

void FindingBoxBound() {

	int Nstep_cr;
	int start[SPACEDIM], end[SPACEDIM], range[SPACEDIM], disp[SPACEDIM];
	int* iterRng = (int *) malloc(2 * SPACEDIM *sizeof(int));
	int* BlockIterRngWhole = IterRngWhole();
	int* BlockIterLocal = (int *) malloc (2* SPACEDIM *sizeof(int));
	Real dx = dxLBM();
	Real xb[SPACEDIM];



	std::numeric_limits<Real>::max();
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		xBoundLocal[2 * SPACEDIM * blockIndex] = std::numeric_limits<Real>::max();
		xBoundLocal[2 * SPACEDIM * blockIndex + 1] = -1.0 * std::numeric_limits<Real>::max();
		xBoundLocal[2 * SPACEDIM * blockIndex + 2] = std::numeric_limits<Real>::max();
		xBoundLocal[2 * SPACEDIM * blockIndex + 3] = -1.0 * std::numeric_limits<Real>::max();
		if (SPACEDIM == 3) {
			xBoundLocal[2 * SPACEDIM * blockIndex + 4] = std::numeric_limits<Real>::max(); //10000000.0;
			xBoundLocal[2 * SPACEDIM * blockIndex + 5] = -std::numeric_limits<Real>::min();
		}
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		BlockIterLocal[0] = BlockIterRngWhole[2 * SPACEDIM * blockIndex];
		BlockIterLocal[1] = BlockIterRngWhole[2 * SPACEDIM * blockIndex + 1];
		BlockIterLocal[2] = BlockIterRngWhole[2 * SPACEDIM * blockIndex + 2];
		BlockIterLocal[3] = BlockIterRngWhole[2 * SPACEDIM * blockIndex + 3];
		if (SPACEDIM == 3) {
				BlockIterLocal[4] = BlockIterRngWhole[2 * SPACEDIM * blockIndex + 4];
				BlockIterLocal[5] = BlockIterRngWhole[2 * SPACEDIM * blockIndex + 5];
		}

		ops_get_abs_owned_range(g_Block[blockIndex], BlockIterLocal, start, end, disp);

		for (int idir = 0; idir < SPACEDIM; idir++) {
			Nf[2 * SPACEDIM * blockIndex + 2 * idir] = start[idir];
			Nf[2 * SPACEDIM * blockIndex + 2 * idir + 1] = end[idir];
		}

//		ops_printf("Nf [%d %d %d %d %d]\n",Nf[0], Nf[1], Nf[2], Nf[3], Nf[4], Nf[5]);
//		printf("Start %d %d %d\n", start[0], start[1], start[2]);
//		printf("End %d %d %d\n", end[0], end[1], end[2]);

		//1st point finding (xmin, ymin, zmin)


		iterRng[0] = Nf[0];
		iterRng[1] = Nf[0] + 1;
		iterRng[2] = Nf[2];
		iterRng[3] = Nf[2] + 1;
		if (SPACEDIM==3) {
			iterRng[4] = Nf[4];
			iterRng[5] = Nf[4] + 1;
		}
		for (int idir = 0; idir < SPACEDIM; idir++)
			xb[idir] = -10000.0;

		ops_par_loop(KerCartBounds,"KerCartBounds", g_Block[blockIndex], SPACEDIM, iterRng,
					 ops_arg_gbl(xb,SPACEDIM, "double", OPS_READ),
					 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "double", OPS_READ));


		xBoundLocal[2 * SPACEDIM * blockIndex] = xb[0] - 0.5 * dx;
		xBoundLocal[2 * SPACEDIM * blockIndex + 2] = xb[1] - 0.5 * dx;

		if (SPACEDIM==3)
			xBoundLocal[2 * SPACEDIM * blockIndex + 4] = xb[2] - 0.5 * dx;

		//2nd point (xmax,ymax, zmax)
		iterRng[0] = Nf[1] - 1;
		iterRng[1] = Nf[1];
		iterRng[2] = Nf[3] - 1;
		iterRng[3] = Nf[3];
		if (SPACEDIM==3) {
			iterRng[4] = Nf[5] - 1;
			iterRng[5] = Nf[5];
		}

		for (int idir = 0; idir < SPACEDIM; idir++)
			xb[idir] = -10000.0;

		ops_par_loop(KerCartBounds,"KerCartBounds", g_Block[blockIndex], SPACEDIM, iterRng,
					 ops_arg_gbl(xb,SPACEDIM, "double", OPS_READ),
				     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "double", OPS_READ));

		xBoundLocal[2 * SPACEDIM * blockIndex + 1] = xb[0] + 0.5 * dx;
		xBoundLocal[2 * SPACEDIM * blockIndex + 3] = xb[1] + 0.5 * dx;
		if (SPACEDIM == 3)
			xBoundLocal[2 * SPACEDIM * blockIndex + 5] = xb[2] + 0.5 * dx;

//		printf("proc  %d Nf [%d %d] [%d %d] [%d %d]\n", ops_get_proc(), Nf[0], Nf[1], Nf[2], Nf[3], Nf[4], Nf[5]);
//		printf("proc %d xbound_local x=[%f %f] y=[%f %f] z=[%f %f]\n", ops_get_proc(), xBoundLocal[2 * SPACEDIM * blockIndex],
//				   xBoundLocal[2 * SPACEDIM * blockIndex + 1], xBoundLocal[2 * SPACEDIM * blockIndex + 2], xBoundLocal[2 * SPACEDIM * blockIndex + 3],
//				   xBoundLocal[2 * SPACEDIM * blockIndex + 4], xBoundLocal[2 * SPACEDIM * blockIndex + 5]);

#ifdef OPS_MPI
		Real minMaxTmp[SPACEDIM];
		Real minMaxLoc[SPACEDIM];

		//x-
		minMaxLoc[0] = xBoundLocal[2 * LATTDIM * blockIndex + 0];
		minMaxLoc[1] = xBoundLocal[2 * LATTDIM * blockIndex + 2];
		if (SPACEDIM == 3)
			minMaxLoc[2] = xBoundLocal[2 * LATTDIM * blockIndex + 4];

		MPI_Allreduce(minMaxLoc, minMaxTmp, SPACEDIM, MPI_DOUBLE, MPI_MIN, OPS_MPI_GLOBAL);

		x_bounds[2 * LATTDIM * blockIndex] = minMaxTmp[0];
		x_bounds[2 * LATTDIM * blockIndex + 2] = minMaxTmp[1];
		x_bounds[2 * LATTDIM * blockIndex + 4] = minMaxTmp[2];

		//Finding x_bounds maximum values
		minMaxLoc[0] = xBoundLocal[2 * LATTDIM * blockIndex + 1];
		minMaxLoc[1] = xBoundLocal[2 * LATTDIM * blockIndex + 3];
		if (SPACEDIM == 3)
			minMaxLoc[2] = xBoundLocal[2 * LATTDIM * blockIndex + 5];

		MPI_Allreduce(minMaxLoc, minMaxTmp, SPACEDIM, MPI_DOUBLE, MPI_MAX, OPS_MPI_GLOBAL);

		x_bounds[2 * LATTDIM * blockIndex + 1] = minMaxTmp[0];
		x_bounds[2 * LATTDIM * blockIndex + 3] = minMaxTmp[1];
		if (SPACEDIM == 3)
			x_bounds[2 * LATTDIM * blockIndex + 5] = minMaxTmp[2];
#else
		for (int iIndex = 0; iIndex < 2 * SPACEDIM; iIndex++) {
			x_bounds[2 * LATTDIM * blockIndex + iIndex] = xBoundLocal[2 * LATTDIM * blockIndex + iIndex];
		}
#endif

//		printf("Block: %d Rank %d x_bounds[%f %f] [%f %f] [%f %f]\n",blockIndex, ops_get_proc(), x_bounds[2 * LATTDIM * blockIndex],
//				x_bounds[2 * LATTDIM * blockIndex + 1], x_bounds[2 * LATTDIM * blockIndex + 2], x_bounds[2 * LATTDIM * blockIndex +3],
//				x_bounds[2 * LATTDIM * blockIndex +4], x_bounds[2 * LATTDIM * blockIndex +5]);
	}

	free(iterRng);
	free(BlockIterLocal);


}

