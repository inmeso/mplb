/*
 * read_data.h
 *
 *  Created on: Sep 24, 2020
 *      Author: jpd38567
 */

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

#include <string>
#include <vector>
#include "type.h"

class BoundaryData {
	public:
		std::vector<int> DirBc;  //Direction of boundary walls
		std::vector<BoundaryScheme> BoundaryType; //Type of boundary
		std::vector<std::vector<double>> velocities;
		std::vector<BoundarySurface> BoundarySurfaces;
		BoundaryData();
		void ImportDirection(int * dir, int spacedim);
		void ImportBoundaryData(int *BType, double* vels);
		void PrintBoundaries();
};

class configData {
	public:
		std::string casename;
		std::string storageOption;
		SizeType blockNum;
		int spacedim;
		int muiFlag;
		Real dx;
		std::vector<SizeType> gridBlock;
		std::vector<Real> gridStart;
		SizeType restartStep;
		double tauRef;
		double Rmax;
		double gamma;
		SizeType iters;
		int exportSteps;
		BoundaryData boundaries;
		configData();
		void updateGrid(int* grid, int gridSize);
		void printData();
		void updateBoundaries(int *dir, int *BType, double* vels, int spacedim);
		void updateStart(Real* gridStart);
};



void ReadInputFile(configData  &data,long int timestep = 0);
#endif /* APPS_LBM_DEM_CLASS_INTERFACE_READ_DATA_H_ */
