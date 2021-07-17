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

/*! @brief   Implementing functions related to the flow field
 * @author  Jianping Meng
 * @details Implementing functions related to create the flow
 * field (allocate memory), set up the geometry and the boundary
 * property, and deallocate the memory.
 */
#ifndef FLOWFIELD_H
#define FLOWFIELD_H
#include <map>
#include <string>
#include "type.h"
#include "block.h"
#include "field.h"


/*!
 * This module is set for defining blocks and variables defined on a block
 * including distribution functions, macroscopic variables, node properties,
 * and relevant parameters.
 * The responsibilities including:
 * 1. Create all variables from files or annually written subroutines
 * 2. Initialise the required macroscopic variables and thereby the
 *    distribution functions.
 * 3. Provide some tools for accessing variables.
 */
const BlockGroup& g_Block();
RealField& g_f();
RealField& g_fStage();
RealFieldGroup& g_MacroVars();
RealFieldGroup& g_MacroVarsCopy();

RealFieldGroup& g_MacroBodyforce();

RealField& g_CoordinateXYZ();
IntFieldGroup& g_NodeType();
IntField& g_GeometryProperty();
Real TimeStep();
const Real* pTimeStep();
Real GetDx();
void SetGridSize(Real meshSize);
const std::string& CaseName();
Real TotalMeshSize();
const std::map<std::string,ops_halo_group>& HaloGroups();
void SetTimeStep(Real dt);

/*!
 * the residual error for steady flows
 * for each macroscopic variable, there are two values: the absolute
 * and relative
 * for each component of a vector, two values are allocated
 */
std::map<int, Real>& g_ResidualError();
std::map<int, ops_reduction>& g_ResidualErrorHandle();

// TODO This function is temporary, will be removed in the near future
//void AllocateMemory();
void WriteFlowfieldToHdf5(const SizeType timeStep);
void WriteDistributionsToHdf5(const SizeType timeStep);
void WriteNodePropertyToHdf5(const SizeType timeStep);
void Partition();
void PrepareFlowField();
// caseName: case name
// spaceDim: 2D or 3D application
void DefineCase(const std::string& caseName, const int spaceDim,
                const bool transient=false);
Real GetMaximumResidual(const SizeType checkPeriod);
// blockNum: total number if blocks.
// blockSize: array of integers specifying the block blocksize.
// meshSize: The size of mesh i.e. dx (At present dx = dy = dz).
// startPos: Starting position of each block.
// void DefineBlocks(const SizeType blockNum,
//                   const std::vector<int>& blockSize, const Real meshSize,
//                   const std::vector<Real>& startPos);

void DefineBlocks(const std::vector<int>& blockIds,
                  const std::vector<std::string>& blockNames,
                  const std::vector<int>& blockSizes, const Real meshSize,
                  const std::map<int, std::vector<Real>>& startPos);
bool IsTransient();

void CalcResidualError();
void DispResidualError(const int iter, const SizeType checkPeriod);
void CopyDistribution(RealField& fDest, RealField& fSrc);
void CopyBlockEnvelopDistribution(Field<Real>& fDest, Field<Real>& fSrc);
void NormaliseF(Real* ratio);
void CopyCurrentMacroVar();
void SetBulkandHaloNodesType(const Block& block, int compoId);
void SetBoundaryNodeType();
void SetBlockGeometryProperty(const Block& block);
void AssignCoordinates(const Block& block,
                       const std::vector<std::vector<Real>>& blockCoordinates);
void UpdateMacroscopicBodyForce(const Real time);
void SetInitialMacrosVars();
void DefineBlockConnection(const std::vector<int>& fromBlock,
                           const std::vector<BoundarySurface>& fromSurface,
                           const std::vector<int>& toBlock,
                           const std::vector<BoundarySurface>& toSurface,
                           const std::vector<VertexType>& connectionType);
void TransferHalos();
#endif
