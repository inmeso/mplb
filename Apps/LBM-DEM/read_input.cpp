/*
 * read_input.cpp
 *
 *  Created on: Sep 24, 2020
 *      Author: jpd38567
 */



#include "read_input.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
//#include "flowfield.h"

void ReadBoundaryConditions(configData &param) {

	int spacedim;
	int BCDir[param.spacedim];
	int bcsize = 2 * param.spacedim;
	int BCType[bcsize];
	int size = param.spacedim * param.spacedim * 2;
	double uWall[size];

	for (int iW = 0; iW < size; iW++)
		uWall[iW] = 0;


	if (ops_get_proc()==0) {

		printf("Reading boundary conditions\n");
		FILE *cfilex;
		cfilex = fopen("input_boundaries.txt","r");
		if (cfilex == NULL) {
			printf("FILE CANNOT BE OPENED\n");
			exit (EXIT_FAILURE);
		}
		fscanf(cfilex, "%d", &spacedim);

		if (spacedim != param.spacedim) {
			printf("MP-LBM ERROR: Dimensional space in boundary conditions not consistent\n");
			printf("Actual SPACEDIM: %d and Boundary condition spacedim: %d", param.spacedim, spacedim);
			exit(EXIT_FAILURE);
		}
		for (int iDim = 0; iDim < spacedim; iDim++) {
		//Read first direction walls
			printf("Reading walls in %d direction\n",iDim);
			fscanf(cfilex, "%d", &BCDir[iDim]);
			fscanf(cfilex, "%d", &BCType[2 * iDim]);
#ifdef OPS_3D
			fscanf(cfilex, "%lf %lf %lf", &uWall[2 * spacedim * iDim], &uWall[2 * spacedim * iDim +1], &uWall[2 * spacedim * iDim + 2]);
#endif

#ifdef OPS_2D
			fscanf(cfilex, "%lf %lf",  &uWall[2 * spacedim * iDim], &uWall[2 * spacedim * iDim +1] );
#endif

			//2nd wall
			fscanf(cfilex, "%d", &BCType[2 * iDim + 1]);

#ifdef OPS_3D
			fscanf(cfilex, "%lf %lf %lf", &uWall[2 * spacedim * iDim + 3], &uWall[2 * spacedim * iDim + 4], &uWall[2 * spacedim * iDim + 5]);
#endif

#ifdef OPS_2D
				fscanf(cfilex, "%lf %lf",  &uWall[2 * spacedim * iDim + 2], &uWall[2 * spacedim * iDim + 3] );
#endif
		}

		fclose(cfilex);
	}

	spacedim = param.spacedim;
	MPI_Bcast(BCDir, spacedim, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(BCType, bcsize, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(uWall, size, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);

	param.updateBoundaries(BCDir, BCType, uWall, spacedim);





}


void ReadInputFile(configData &param, long int timestep) {

	char tmp[101], tmpStor[101];
	int blockSize, spacedim, size, muiInterfaceFlag, iters, exportStep;
	int *grid;
	Real *gridStart;

	std::vector<int> gridTmp;
	std::vector<Real> startTmp;
	double Lx, tauRef, RadTmp, gamma;

	if (ops_get_proc()==0) {
		printf("Rank %d: Ready to read input parameters\n", ops_get_proc());
		FILE *cfilex;



		cfilex = fopen("input_params.txt","r");
		if (cfilex == NULL) {
			printf("FILE CANNOT BE OPENED\n");
			exit (EXIT_FAILURE);
		}

		fscanf(cfilex,"%100s",tmp);

		printf("I read the filename\n");

		fscanf(cfilex, "%d", &blockSize);
		fscanf(cfilex, "%d", &spacedim);

		int temp;
		Real temp1;
		for (int iBlock = 0; iBlock < blockSize; iBlock++) {
			for (int iDir = 0; iDir < spacedim; iDir++)  {
				fscanf(cfilex, "%d", &temp);
				gridTmp.push_back(temp);
			}
			for (int iDir  = 0; iDir < spacedim; iDir++) {
				fscanf(cfilex, "%lf", &temp1);
				startTmp.push_back(temp1);
			}
		}

		fscanf(cfilex, "%lf", &Lx);
		fscanf(cfilex, "%d", &muiInterfaceFlag);
		fscanf(cfilex, "%lf", &tauRef);
		fscanf(cfilex, "%lf", &RadTmp);
		fscanf(cfilex, "%d", &iters);
		fscanf(cfilex, "%d", &exportStep);
		fscanf(cfilex, "%lf", &gamma);

		fscanf(cfilex, "%100s", tmpStor);

		fclose(cfilex);
		//Broadcasting the options

	}

	MPI_Bcast(tmp, 100, MPI_CHAR, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&blockSize, 1 , MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&spacedim, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&Lx, 1, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&muiInterfaceFlag, 1, MPI_INT,0,  OPS_MPI_GLOBAL);
	MPI_Bcast(&tauRef, 1, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&RadTmp, 1, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&iters, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&exportStep, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(&gamma, 1, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(tmpStor, 100, MPI_CHAR, 0, OPS_MPI_GLOBAL);

//	printf("I passed broadcastingStep\n");

	//Ready to add to the class object
	//printf("Rank %d: Spacedim: %d blocks: %d\n", ops_get_proc(), spacedim, blockSize);
	size = spacedim * blockSize;
	grid = new int[size];
	gridStart = new Real[size];

	if (ops_get_proc()==0) {
		for (int iLoc = 0; iLoc < size; iLoc++)
			grid[iLoc] = gridTmp[iLoc];

		for (int iLoc = 0; iLoc < size; iLoc++)
			gridStart[iLoc] = startTmp[iLoc];

		for (int iLoc = 0; iLoc < size; iLoc++)
			printf("%d grid =%d starts at %f\n", iLoc, grid[iLoc], gridStart[iLoc]);

	}
	MPI_Bcast(grid, size, MPI_INT, 0, OPS_MPI_GLOBAL);
	MPI_Bcast(gridStart, size, MPI_DOUBLE, 0, OPS_MPI_GLOBAL);




	param.casename = tmp;

	//Sanity checks and data extraction
	if (blockSize < 1) {
		ops_printf("MP-LBM ERROR: Number of blocks smaller than zero\n");
		exit(EXIT_FAILURE);
	}
	else
		param.blockNum = blockSize;

	if (spacedim < 2 || spacedim > 3) {
		ops_printf("MP-LBM ERROR: Current dimension %d not supported\n", spacedim);
		exit(EXIT_FAILURE);
	}
	else
		param.spacedim = spacedim;

	if (tauRef < 0.0) {
		ops_printf("MP-LBM ERROR: Relaxation time smaller than zero\n");
		exit(EXIT_FAILURE);
	}
	else
		param.tauRef = tauRef;

	if (Lx <= 0.0) {
		ops_printf("MP-LBM ERROR: The size of the simulation domain is not finite\n");
		exit(EXIT_FAILURE);
	}
	else
		param.dx = Lx;

	if (RadTmp < 0.0) {
		ops_printf("MP-LBM ERROR: Non positive radius\n");
		exit(EXIT_FAILURE);
	}
	else
		param.Rmax = RadTmp;

	param.storageOption = tmpStor;




	param.updateGrid(grid, size);

	//printf("Rank %d: UpdatedGrid() passed\n", ops_get_proc());
	param.updateStart(gridStart);
	//printf("Rank %d: UpdateStart() passed\n", ops_get_proc());

	if (muiInterfaceFlag < 0 && muiInterfaceFlag > 1) {
		ops_printf("MP-LBM ERROR: The flag for mui should be 0 or 1\n");
		exit(EXIT_FAILURE);
	}
	else
		param.muiFlag = muiInterfaceFlag;

	if (timestep < 0 ) {
		ops_printf("MP-LBM ERROR: The restart step is smaller than zero\n");
		exit(EXIT_FAILURE);
	}
	else
		param.restartStep = timestep;

	if (iters < 0) {
		ops_printf("MP-LBM ERROR: Non positive number of iterations\n");
		exit(EXIT_FAILURE);
	}
	else
		param.iters = static_cast<SizeType>(iters);

	if (exportStep < 0) {
		ops_printf("MP-LBM ERROR: Non positive export steps\n");
		exit(EXIT_FAILURE);
	}
	else
		param.exportSteps = exportStep;

	if (gamma < 0) {
		ops_printf("MP-LBM ERROR: Non positive gamma. Gamma set to zero\n");
		param.gamma = 0.0;
	}
	else
		param.gamma = gamma;

	ReadBoundaryConditions(param);

	param.boundaries.PrintBoundaries();
	delete[] grid;
	delete[] gridStart;
	//We

}

configData::configData() {

	spacedim = 2;
	tauRef = 0.001;
	blockNum = 1;
	dx = 0.00001;
	muiFlag = 0;
	restartStep = 0;
	Rmax = 1.0;
	gamma = 0.0;
}

void configData::updateGrid(int* grid, int gridSize) {

	int size = spacedim * blockNum;
	if (gridSize != size) {
		ops_printf("MP-LBM ERROR: The size of imported not equal to internal size\n");
		exit(EXIT_FAILURE);
	}


	for (int iDir = 0; iDir < size; iDir++) {
		if (grid[iDir] > 0)
			gridBlock.push_back(grid[iDir]);
		else {
			ops_printf("MP-LBM ERROR: Non positive value for block %d in %d direction\n", iDir/(spacedim * blockNum), iDir - spacedim * blockNum);
			exit(EXIT_FAILURE);
		}
	}

	dx /= (static_cast<double>(gridBlock[0]) - 1.0);

}


void configData::printData() {
	printf("Simulation case: %s\n", casename.c_str());
	printf("Rank %d: The dimensional space is %d and %d blocks employed\n",ops_get_proc(), spacedim, blockNum);
	printf("Rank %d: A uniform grid of dx = %f is employed\n",ops_get_proc(), dx);
	for (int iDir = 0; iDir < blockNum; iDir++)  {
		printf("Rank %d [",ops_get_proc());
		for (int ix =0; ix < spacedim; ix++)
			printf("%d, ", gridBlock[iDir * spacedim + ix]);
		printf("] ");

		printf("which starts at[");
		for (int ix = 0; ix <spacedim; ix++)
			printf("%f, ", gridStart[iDir * spacedim + ix]);
		printf("] ");
	}

	printf("Rank %d: Mui interface: %d and restart timestep %d\n",ops_get_proc(), muiFlag, restartStep);

}

void configData::updateStart(Real* start) {

	if (gridBlock.size()==0) {
		ops_printf("MP-LBM Error: Block grid size must be defined first\n");
		exit(EXIT_FAILURE);
	}


	for (int iDir = 0; iDir < gridBlock.size(); iDir++)
		gridStart.push_back(start[iDir]);



}


void configData::updateBoundaries(int *dir, int *BType, double* vels, int spacedim) {

	boundaries.ImportDirection(dir, spacedim);
	boundaries.ImportBoundaryData(BType, vels);
}

BoundaryData::BoundaryData() {

}

void BoundaryData::ImportDirection(int * dir, int spacedim) {

	for (int idir = 0; idir < spacedim; idir++)
		DirBc.push_back(dir[idir]);

}

void BoundaryData::ImportBoundaryData(int *BType, double* vels) {

	int spacedim;
	int Direction;
	if (DirBc.empty()) {
		ops_printf("MP-LBM: Function import boundary type must be called after import direction\n");
		exit(EXIT_FAILURE);
	}
	else
		spacedim = DirBc.size();

	for (int iDir = 0; iDir < spacedim; iDir++) {
		Direction = DirBc[iDir];
		switch(Direction) {
			case 0 : {
				BoundarySurfaces.push_back(BoundarySurface_Left);
				BoundarySurfaces.push_back(BoundarySurface_Right);
				break;
			}
			case 1 : {
				BoundarySurfaces.push_back(BoundarySurface_Bottom);
				BoundarySurfaces.push_back(BoundarySurface_Top);
				break;
			}
			case 2 : {
				BoundarySurfaces.push_back(BoundarySurface_Front);
				BoundarySurfaces.push_back(BoundarySurface_Back);
				break;
			}
			default : {
				ops_printf("MP-LBM: Undefined case\n");
				exit(EXIT_FAILURE);
				break;
			}
		}

	}

	std::vector<double> uVels;

	for (int iDir = 0; iDir < spacedim; iDir++) {
		BoundaryScheme temp = static_cast<BoundaryScheme> (BType[2 * iDir]);
		BoundaryType.push_back(temp);

		for (int jDir = 0; jDir < spacedim ; jDir++) {
			uVels.push_back(vels[2 * iDir * spacedim + jDir]);
		}

		velocities.push_back(uVels);


		//Second wall
		uVels.clear();
		temp = static_cast<BoundaryScheme>(BType[2 * iDir + 1]);
		BoundaryType.push_back(temp);

		for (int jDir = 0; jDir < spacedim; jDir++)
			uVels.push_back(vels[(2 * iDir + 1) * spacedim + jDir]);

		velocities.push_back(uVels);
		uVels.clear();

	}

}


void BoundaryData::PrintBoundaries() {

	int rank = ops_num_procs() - 1;


	if (rank == ops_get_proc()) {


		if (DirBc.empty()) {
			printf("MP-LBM: BoundaryData not defined properly");
			return;
		}
#ifdef OPS_3D
		printf("--------------------------------------------------------------------------------------------------------------------------------\n");
		printf("Printing boundaries at rank %d\n", ops_get_proc());

		printf("Boundary 1: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[0], BoundarySurfaces[0], BoundaryType[0],
				velocities[0][0], velocities[0][1],velocities[0][2]);

		printf("Boundary 2: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[0], BoundarySurfaces[1], BoundaryType[1],
				velocities[1][0], velocities[1][1],velocities[1][2]);

		printf("Boundary 3: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[1], BoundarySurfaces[2], BoundaryType[2],
				velocities[2][0], velocities[2][1],velocities[2][2]);

		printf("Boundary 4: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[1], BoundarySurfaces[3], BoundaryType[3],
				velocities[3][0], velocities[3][1], velocities[3][2]);

		printf("Boundary 5: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[2], BoundarySurfaces[4], BoundaryType[4],
				velocities[4][0], velocities[4][1], velocities[4][2]);

		printf("Boundary 6: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f %f]\n", DirBc[2], BoundarySurfaces[5], BoundaryType[5],
				velocities[5][0], velocities[5][1], velocities[5][2]);
#endif

#ifdef OPS_2D
		printf("Boundary 1: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f]\n", DirBc[0], BoundarySurfaces[0], BoundaryType[0],
				velocities[0][0], velocities[0][1]);


		printf("Boundary 2: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f]\n", DirBc[0], BoundarySurfaces[1], BoundaryType[1],
				velocities[1][0], velocities[1][1]);

		printf("Boundary 3: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f]\n", DirBc[1], BoundarySurfaces[2], BoundaryType[2],
				velocities[2][0], velocities[2][1]);

		printf("Boundary 3: Dir %d wall location is %d boundary scheme %d and Uw=[%f %f]\n", DirBc[1], BoundarySurfaces[3], BoundaryType[3],
				velocities[3][0], velocities[3][1]);
#endif

		printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	}
}
