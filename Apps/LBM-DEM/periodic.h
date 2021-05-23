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
 * @brief   Functions for the implementation of periodic conditions at DEM-LBM simulations
 * @author  Chrysovalantis Tsigginos
 * @details Wrap and Kernel functions for the implementation of periodic B.C
 */

#ifndef PERIODIC_H_
#define PERIODIC_H_

#include "type.h"
#include "flowfield.h"
#include "psm.h"
#include "memory_handle.h"
#include "evolution3d.h"
#include "evolution.h"


extern int periodicFlag; //0: No Periodic B.C in any  direction 1: Periodic B.C at least in one direction
extern int periodic[3]; //Periodic Flag for each direction 1: Periodic B.C 0: No Periodic B.C

void DefinePeriodicHaloTransfer();
void PeriodicHaloTransfer();
void InitializePeriodic();
void DestroyPeriodic();
void ImplementBoundary3DPeriodic();

void InitializePeriodicDragForce();

//Two dimension functions
void SetPeriodicSize2D(Real x1,Real y1);
void SetCutOff2D(Real xcut, Real ycut);

//Three dimension functions
void SetPeriodicSize3D(Real x1,Real y1, Real z1);
void SetupVertexType(VertexType* vertexType);
void SetCutOff3D(Real xcut, Real ycut, Real zcut);

void DefinePeriodicBoundaries(BoundaryScheme* boundaryType);
void DefinePeriodicBoundariesRestart(std::vector<BoundaryScheme>& types,std::vector<int>& dir);

bool decide(Real xp, Real xlo, Real xhi);

void PeriodicPartition();
void PeriodicExchange();
void ForwardComm();
void ReverseComm();
//Packing functions
void packForward2D(int blockIndex, int nsend, int iswap, int* sendList, Real* buff, int* periodic, int* pbc);
void packForward3D(int blockIndex, int nsend, int iswap, int* sendList,Real* buff,int* periodic,int * pbc);

void unpackComm2D(int blockIndex, int nrecv, int first, Real *buf);
void unpackComm3D(int blockIndex, int nrecv, int first, Real *buf);

void packReverse2D(int blockIndex, int nParticles, int nfirst, Real* buff);
void packReverse3D(int blockIndex, int nParticles, int nfirst, Real* buff);
void unPackReverse2D(int blockIndex, int nTot, int* sendList, Real* buff);
void unPackReverse3D(int blockIndex, int nTot, int* sendList, Real* buff);


#endif /* APPS_LBM_DEM_PERIODIC_H_ */
