/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*!
 * @brief   Wrap functions for DEM-LBM iteration scheme
 * @author  Chrysovalantis Tsigginos
 * @details Define functions for DEM-LBM iteration scheme
 * */

#include "iterate_dem.h"

SizeType MAXITER;
int CHECKPERIOD;

void IterateDEMLBM(const SizeType maxIter, const SizeType checkPeriod, const SizeType timeStart, int savingFlag) {

	Real dtD; //TimeStep for force Calculation
	CHECKPERIOD = checkPeriod;
	MAXITER = maxIter;
	SizeType initPoint;
	if (Nfl <= 1)
		dtD = TimeStep();
	else
		dtD  = DemTimeStep();
	int flag;
	ops_printf("MP-LBM: Number of iterations: %d\n", maxIter);
	ops_printf("Npl = %d Nfl= %d\n", Npl, Nfl);


	double ct0I, ct1I, et0I, et1I;
	double et0MUI, et1MUI;
	double timeMUI = 0.0;

	ops_timers(&ct0I, &et0I);


	//exit(EXIT_FAILURE);
	initPoint = timeStart;
	SizeType iter = timeStart;
	for (SizeType iter1 = timeStart; iter1 < MAXITER; iter1++) {

		++iter;
		if (iter1 == initPoint)
			flag = 1;
		else
			flag = 0;
		//ops_printf("Iteration %d\n", iter);

		if ( (iter% static_cast<SizeType>(Npl)) == 0) {
			if (muiInterfaceFlag == 1) {
				ops_timers(&ct0I, &et0MUI);

				ExtractParticleData(iter);

				ops_timers(&ct1I, &et1MUI);
				timeMUI += et1MUI - et0MUI;
			}
			InitializeDragForce();
			for (int jfl = 1; jfl < Nfl + 1; jfl++) {
				//Here we need to introduce DEM-LBM type scheme-since there are many
#ifdef OPS_3D
				StreamCollisionPSM3D(flag);
#endif

#ifdef OPS_2D
				StreamCollisionPSM2D();
#endif
			}
			CalculateDragForce(dtD, iter);
			if (muiInterfaceFlag == 1) {

				ops_timers(&ct0I, &et0MUI);

				SendParticleData(iter);
				//ForgetData(iter);

				ops_timers(&ct1I, &et1MUI);
				timeMUI += et1MUI - et0MUI;
			}
		}

		if (savingFlag == 1) {
		if (iter % checkPeriod == 0) {
#ifdef OPS_3D
			UpdateMacroVars3D();
				if (forceFlag == 1)
					UpdateMacroVarsForceFluidSolid();

				CorrectFluidVelocities();


				CalcResidualError3D();
				DispResidualError3D(iter, CHECKPERIOD * TimeStep());
				WriteFlowfieldToHdf5(iter);
				WriteDistributionsToHdf5(iter);
				WriteNodePropertyToHdf5(iter);
				WritePSMVariablesToHdf5(iter);
#endif

#ifdef OPS_2D
				UpdateMacroVars();
				CalcResidualError();
				DispResidualError(iter, CHECKPERIOD * TimeStep());
				WriteFlowfieldToHdf5(iter);
				WriteDistributionsToHdf5(iter);
				WriteNodePropertyToHdf5(iter);
				WritePSMVariablesToHdf5(iter);
#endif
			}
		}

	}

	ops_timers(&ct1I, &et1I);
	ops_printf("--------------------------------------------------------------------------\n");
	ops_printf("Total running time: %12.9f\n", et1I-et0I);
	ops_printf("MUI Running time: %12.9f\n", timeMUI);
	ops_printf("--------------------------------------------------------------------------\n");




#ifdef OPS_3D
	UpdateMacroVars3D();
	if (forceFlag == 1)
		UpdateMacroVarsForceFluidSolid();

	CorrectFluidVelocities();
	WriteFlowfieldToHdf5(MAXITER);
	WriteDistributionsToHdf5(MAXITER);
	WriteNodePropertyToHdf5(MAXITER);
	WritePSMVariablesToHdf5(MAXITER);

#endif

	if (periodicFlag==1)
		DestroyPeriodic();

	DestroyModel();
	DestroyFlowfield();
	DestroyParticleParams();
}


void IterateDEMLBMSS(const Real convergenceCriteria, const int checkPeriod, const int maxIter, int savingFlag) {

	CHECKPERIOD = checkPeriod;
	int iter =0;
	Real residualError = 10000;
	Real dt = TimeStep();
	int flag;
	ops_printf("maxIter = %d iter1 = %d CHECKPERIOD = %d\n", maxIter, iter,	CHECKPERIOD);

	while (residualError >= convergenceCriteria) {
		if (iter == 0)
			flag = 1;
		else
			flag = 0;
		InitializeDragForce();
#ifdef OPS_3D
		StreamCollisionPSM3D(flag);
		if ((iter % CHECKPERIOD) == 0) {
			UpdateMacroVars3D(); //Add correction term

			if (forceFlag == 1)
				UpdateMacroVarsForceFluidSolid();


			CorrectFluidVelocities();

			CalcResidualError3D();
			residualError = GetMaximumResidual(checkPeriod);

			DispResidualError3D(iter,static_cast<double>( CHECKPERIOD * TimeStep()));
		}
#endif

#ifdef OPS_2D
		StreamCollisionPSM2D();
		if ((iter1 % CHECKPERIOD) == 0) {
			UpdateMacroVars();

			if (forceFlag == 1)
				UpdateMacroVarsForceFluidSolid();


			CorrectFluidVelocities();

			CalcResidualError();
			residualError = GetMaximumResidual(checkPeriod);
			DispResidualError(iter, CHECKPERIOD * TimeStep());
		}
#endif
		CalculateDragForce(dt, iter);
		ops_printf("Residual = %e,  Convergence target = %e\n", residualError, convergenceCriteria);

		iter += 1;


		if (iter > maxIter)
			break;
	}
	if (savingFlag == 1) {
		WriteFlowfieldToHdf5(0);
		WriteDistributionsToHdf5(0);
		WriteNodePropertyToHdf5(0);
		WritePSMVariablesToHdf5(0);
	}
}
