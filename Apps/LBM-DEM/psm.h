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
 * @brief   Functions for the implementation of the PSM scheme
 * @author  Chrysovalantis Tsigginos
 * @details Wrap and Kernel functions for the implementation of the PSM scheme
 */

#ifndef PSM_H_
#define PSM_H_


#include "flowfield.h"
#include "scheme.h"

#include "type.h"

#include "evolution.h"

#include "boundary.h"
#include "dem_handle3d.h"
#include "model.h"

#include "evolution3d.h"
#include "periodic.h"

extern int NELEM; //maximum  number of particles that intersect with a computational cell
extern Real force[3]; //Body acceleration applied to the fluid
extern int forceFlag; //Flag for calculation of forces in PSM collision

extern Field<Real> g_Fd;
extern Field<Real> g_sfm; //Solid fraction of particle p for the computational cell
extern Field<int> g_id; // Local particle id that interact with computational cell i
extern Field<Real> g_vp; // Particle velocities at computational cell i.

//Wrapper functions
void SetInitialMacroVars(Real* initialMacros);
void DefinePSMVariables(SizeType timeStep = 0, Real gamma = 0.0);
void SetupBodyForces();


//Functions for mapping particle in PSM
void InitializePSMLists();

//Functions for velocity calculation on PSM scheme (2nd order)
void UpdateMacroVarsForceFluidSolid();

//2D functions
void EstimateSolidFractionMP2D();
void UpdateSolidVelocity2D();
void SolidFluidCollision2D();
void EstimateDragForce2ndOrder2D();
void StreamCollisionPSM2D();
//Auxiliary 2D particle mapping function
Real  ComputeCellSF(Real* xf,Real* xp,const Real dx,const int intOrder);
Real EvaluateIntegral(Real x, Real ym, Real Rp, Real dx,int n_int);
void GaussIntegrationTable(double *t, double* W, const int n_int);
double ReturnRoot(double xmin,double xmax, int target, double prod, double Rp, double ym);
double sign(double x);
void EstimateVelocity2D(Real* vf,const Real* vel_p,const Real* xp,Real *xfl);

//3D functions
void EstimateSolidFractionMP3D();
void UpdateSolidVelocity3D();
void SolidFluidCollision3D();
void EstimateDragForce2ndOrder3D();
void StreamCollisionPSM3D(int flag);

void WritePSMVariablesToHdf5(const long TimeStep);
//Auxiliary 3D functions
Real  ComputeSolidFraction3D(const Real* xfl,const Real *xp, Real dx,int nmc);
void EstimateVelocity3D(Real* vf,const Real* vel_p,const Real* xp, Real *xfl);
Real CalculateSolidFractionSpheres(Real* xf,Real Rc,Real* xp,Real Rp);
Real CalcBodyForce2ndOrder(const int xiIndex,const Real rho,const Real* u, const Real*  accel);

//Kernel functions
void KerSetInitialMacroVars(ACC<Real>& macros, const ACC<Real>& coords, const Real* initialMacros);
void KerInitialize(ACC<Real>& sfp, ACC<Real>& vp, ACC<int>& id, ACC<Real>& Fd);
//Particle mapping kernel functions
void KerSetInitialMacroVars(ACC<Real>& macros, const ACC<Real>& coords, const Real* initialMacros);
void KerCalcMacroVarsForceFluidSolid(ACC<Real>& macroVars, const ACC<int>& nodeType, const ACC<Real>& sfp, const Real* dt, const Real* force);

//2D Kernel functions
void KerEvaluateSolidFractionMP2D(ACC<Real>& sfp, ACC<Real>& vp, ACC<int>& id, const ACC<Real>& xf, const int* idp, const Real* xp,
								  const Real* velP, const Real* dx, const int* intOrder);
void KerSolVelocityUpdateMP2D(ACC<Real>& vp, const ACC<Real>& xf, const ACC<int>& id, const int* idp, const Real* xp, const Real* up);
void KerPRE2D(ACC<Real>& fcopy, ACC<Real>& Fd,const ACC<Real>& f, const ACC<Real>& macros,const ACC<int>& nodeType,
					 const ACC<int>& id, const ACC<Real>& sfp, const Real* tau, const Real* dt, const int* forceFlag, const Real* force);
void KerDragForcePRE2D(const ACC<int>&id, const ACC<Real>& sfp, const ACC<Real>&xf, const ACC<Real>&Fd, const int* idp, const Real* xpos,
					 Real* FDp, const Real* dt, const Real* tau);

//3D Kernel functions


void KerEvaluateSolidFractionMP3D(ACC<int>& id, ACC<Real>& sfp, ACC<Real>& vp, const ACC<Real>& xf, const int* idp, const Real* xp,
								  const Real* velP, const Real* dx);
void KerSolVelocityUpdateMP3D(ACC<Real>& vp, const ACC<int>& id, const ACC<Real>& xf, const int* idp, const Real* xp, const Real* up);

void  KerPRE3D(ACC<Real>& fcopy, ACC<Real>& Fd,const ACC<Real>& f, const ACC<Real>& macros,const ACC<int>& nodeType,
		 	 	      const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& vp, const Real* tau, const Real* dt, const int* forceFlag,
		 	          const Real* force);

void KerDragPRE3D(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf, const ACC<Real>& Fd,
					   Real* FDp, const int* idp, const Real* xp, const Real* dt, const Real* tau);

void KerStreamPeriodic3D(ACC<Real> & f, const ACC<Real>& fStage, const ACC<int>& nodeType, const ACC<int>& geometry);
#endif /* PSM_H_ */
