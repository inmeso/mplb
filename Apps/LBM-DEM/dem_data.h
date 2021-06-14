/*
 * interaction_data.h
 *
 *  Created on: Jun 4, 2021
 *      Author: jpd38567
 */

#ifndef DEM_DATA_H_
#define DEM_DATA_H_

#include "type.h"


enum ParticleShapeDiscriptor{spherical, quadratic, mesh};

struct InteractionData {
	SizeType nSteps;
	SizeType nStart;
	bool restartFlag;
	bool muiFlag;
	Real dtDEM;
	Real convergenceRate;
	SizeType Nf; //number of DEM iterations
	SizeType Npl; //Need to give definition
	SizeType checkPeriod;
	SizeType maxIters;
	SizeType checkPeriodStS;
	ParticleShapeDiscriptor particleShape;

};



#endif /* APPS_LBM_DEM_INTERACTION_DATA_H_ */
