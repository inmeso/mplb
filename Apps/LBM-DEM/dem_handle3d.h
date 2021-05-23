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
 * @brief   Functions for particle handling in the LBM code
 * @author  Chrysovalantis Tsigginos
 * @details General functions for direct DEM-LBM simulations
 * 			User must define functions for fluid-particle interaction scheme
 */

#ifndef DEM_HANDLE3D_H_
#define DEM_HANDLE3D_H_

#include <cmath>
#include <string>

#include "dem_particles.h"
#include "type.h"
#include "flowfield.h"
#include "mui.h"
#include "memory_handle.h"
#include "iterate_dem.h"
#include "dem_particles.h"
#include "periodic.h"
#include "mui_interface.h"

extern Real cutoff_distance;  //cut-off distance for building solid fraction list
extern int flag_vel;		 //For violent build of solid fraction list


extern Real dtDEM; //DEM timestep
extern int restartFlag, restartTimeStep; //RestartFlag and TimeStep for
extern int Nfl, Npl; //Number of iterations for LBM and DEM solver before data exchange.

extern Real alpha; // timestep ratio of DEM to LBM

extern int muiInterfaceFlag; //Flag for activation of mui interface


//DEM timestep functions
void SetupDEMLBMcoupling();
Real DemTimeStep();


//Particle and FSI allocation functions
void AllocateParticles(); //Allocate memory for particle lists
void AllocateDomainVariables(); //Allocate domain variables


void SetupDEMLBMcoupling();
void SetupParticles(); //Function for particle initialization
void ReadParticleDetails();
void InitializeDEMLBM(Real convergenceRate, int maxIters, int checkPeriod,SizeType &maxStep, int SavingFlag = 1);
void 	UpdateParticlMappingDragForce(SizeType timeSetup);
void InitializeDEMLBMRestart(SizeType &maxStep, SizeType timeSetup,SizeType RunStep = 1000, int SavingFlag = 1);

//Particle handling auxiliary functions
int CheckDistance(); //Decide for particle re-mapping based on distances
void UpdateOldParticleLocation(); //Particle locations during particle mapping
void InitializeParticleOldLocation(); //Set initial particle locatios
void ParticleMapping(); //Create particle stencil
void CalculateDragForce(Real dt,SizeType iter);
void InitializeDragForce();
void DestroyParticleParams();
//MUI Functions
void DefineMUIInterface(Real Rmax); //Initialization of MUI library
void SetupMUICommunication(int maxIter); //Setup MUI communication
void DestroyMUIInterface();


void ExtractParticles3D(int time); // Particle extraction for three dimensional simulations
void ExtractParticleData(SizeType time); //Extract particle Data through the MUI library

void SendParticleData3D(int t); //Send particle for three dimensional simulations
void SendParticleData(SizeType timestep); //Send data to DEM code

void ForgetData(int t); //Forget data at time t
//Auxiliary MUI functions

bool SphereParallepidedIntersection(Real cx, Real cy, Real cz, Real radius, Real xbMin,
									Real  xbMax, Real ybMin, Real ybMax, Real zbMin, Real zbMax);
#endif
