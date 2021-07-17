#ifndef FORCE_FSI_HOST_DEVICE_H_
#define FORCE_FSI_HOST_DEVICE_H_
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

#include "type.h"
enum ForceType{noForce = 0, forceConstant=1, userDefinedForce = 2};

static inline OPS_FUN_PREFIX void EvaluateNoForce(Real* forceOut, int spaceDim) {

	for (int iDir = 0; iDir < spaceDim; iDir++)
		forceOut[iDir] = 0.0;

}

static inline OPS_FUN_PREFIX void EvaluateForceConstant(Real* forceOut,const Real* forceIn, int spaceDim) {

	for (int iDir = 0; iDir < spaceDim; iDir++)
		forceOut[iDir] = forceIn[iDir];

}

static inline OPS_FUN_PREFIX void EvaluateForceUser(Real* forceOut,const Real* forceIn, Real* xPos,
		int size, int spaceDim) {


	//TO BE DEFINED BY THE USER.


}

static inline OPS_FUN_PREFIX void calculateForce(Real* forceOut,const Real* forceIn, Real* xPos,
		ForceType forceModel, int size, int spaceDim) {

	switch (forceModel) {
		case noForce:
			EvaluateNoForce(forceOut, spaceDim);
			break;
		case forceConstant:
			EvaluateForceConstant(forceOut, forceIn, spaceDim);
			break;
		case userDefinedForce:
			EvaluateForceUser(forceOut, forceIn, xPos, size, spaceDim);
			break;
		default:
			EvaluateNoForce(forceOut, spaceDim);
	}

}



#endif /* APPS_LBM_DEM_NOP_FORCE_FSI_HOST_DEVICE_H_ */
