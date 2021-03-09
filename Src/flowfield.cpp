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

#include "flowfield.h"
#include <type_traits>
#include "block.h"
#include "field.h"
#include "model.h"
#include "boundary.h"
#include "scheme.h"
std::string CASENAME;
bool TRANSIENT{false};
/*!
 * SPACEDIM=2 for 2D 3 for three 3D
 */
#ifdef OPS_3D
int SPACEDIM{3};
#endif  //  OPS_3D
#ifdef OPS_2D
int SPACEDIM{2};
#endif // ops_2D
BlockGroup BLOCKS;
RealField f{"f"};
RealField fStage{"fStage"};
RealField MacroVars{"MacroVars"};
RealField MacroVarsCopy{"MacroVars_Copy"};
Real* g_ResidualError{nullptr};
ops_reduction* g_ResidualErrorHandle{nullptr};
RealField MacroBodyforce{"MacroBodyForce"};

const BlockGroup& g_Block() { return BLOCKS; };
RealField& g_f() { return f; };
RealField& g_fStage() { return fStage; };
RealField& g_MacroVars() { return MacroVars; };
RealField& g_MacroVarsCopy() { return MacroVarsCopy; };
RealField& g_MacroBodyforce() { return MacroBodyforce; };
/**
 * DT: time step
 */
Real DT{1};

/**
 * TAUREF: the reference relaxation time
 * In appropriate non-dimensional system, it is the Knudsen number
 * It must be a constant during the run time
 */
Real* TAUREF{nullptr};
RealField CoordinateXYZ{"CoordinateXYZ"};
RealField& g_CoordinateXYZ() { return CoordinateXYZ; };
std::map<SizeType, std::vector<std::vector<Real>>> COORDINATES;

IntField NodeType{"NodeType"};
IntField GeometryProperty{"GeometryProperty"};
IntField& g_NodeType() { return NodeType; };
IntField& g_GeometryProperty() { return GeometryProperty; };

/*!
 * Formal collection of halo relations required by the OPS library
 */
std::map<std::string,ops_halo_group> HALOGROUPS;

void DefineCase(const std::string& caseName, const int spaceDim,
                const bool transient) {
    if (SPACEDIM != spaceDim) {
        ops_printf("Error! The SPACEDIM here is inconsistent with the\n");
        assert(SPACEDIM == spaceDim);
    }
    SetCaseName(caseName);
    ops_decl_const("SPACEDIM", 1, "int", &SPACEDIM);
    TRANSIENT = transient;
}

bool IsTransient() { return TRANSIENT; }

//This version is for the distribution function
//Other version to follow.
//TODO this function only suitable for single halo
#ifdef OPS_3D
void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              ops_dat dat, const int haloDepth) {
    if (dat == nullptr) {
        ops_printf("Data must be allocated before defining halo relations!");
        assert(dat == nullptr);
    }
    // max halo depths for the dat in the positive direction
    int d_p[3] = {haloDepth, haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[3] = {-haloDepth, -haloDepth, -haloDepth};
    // The domain size in the Block 0
    const Block& block{BLOCKS.begin()->second};
    int nx{(int)block.Size().at(0)};
    int ny{(int)block.Size().at(1)};
    int nz{(int)block.Size().at(2)};
    int dir[] = {1, 2, 3};
    for (auto& pair : haloPair) {
        if (0 == pair.first) {
            // left and right pair
            int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1],
                               nz + d_p[2] - d_m[2]};
            int base_from[] = {0, d_m[1], d_m[2]};
            int base_to[] = {nx, d_m[1], d_m[2]};
            ops_halo leftToRight = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[0] = nx + d_m[0];
            base_to[0] = d_m[0];
            ops_halo rightToLeft = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{leftToRight, rightToLeft};
            HALOGROUPS.emplace(pair.second, ops_decl_halo_group(2, group));
        }

        if (1 == pair.first) {
            // top and bottom pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth,
                               nz + d_p[2] - d_m[2]};
            int base_from[] = {d_m[0], 0, d_m[2]};
            int base_to[] = {d_m[0], ny, d_m[2]};
            ops_halo botToTop = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            base_from[1] = ny + d_m[1];
            base_to[1] = d_m[1];
            ops_halo topToBot = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            ops_halo group[]{botToTop, topToBot};
            HALOGROUPS.emplace(pair.second, ops_decl_halo_group(2, group));
        }

        if (2 == pair.first) {
            // front and back pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], ny + d_p[1] - d_m[1],
                               haloDepth};
            int base_from[] = {d_m[0], d_m[1], 0};
            int base_to[] = {d_m[0], d_m[1], nz};
            ops_halo backToFront = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[2] = nz + d_m[2];
            base_to[2] = d_m[2];
            ops_halo frontToBack = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{backToFront, frontToBack};
            HALOGROUPS.emplace(pair.second, ops_decl_halo_group(2, group));
        }
    }
}

void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
    DefinePeriodicHaloPair3D(haloPair, f[blockId], f.HaloDepth());
}

void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              RealField& data) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
    DefinePeriodicHaloPair3D(haloPair, data[blockId], data.HaloDepth());
}

void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              IntField& data) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
    DefinePeriodicHaloPair3D(haloPair, data[blockId], data.HaloDepth());
}

#endif //OPS_3D

void Partition() {
    ops_partition((char*)"LBM Solver");
    PrepareFlowField();
}

/*
 * We need a name to specify which file to input
 * To be decided: a single filename or an array of filenames
 */

void WriteFlowfieldToHdf5(const SizeType timeStep) {
    MacroVars.WriteToHDF5(CASENAME, timeStep);
    CoordinateXYZ.WriteToHDF5(CASENAME, timeStep);
    MacroBodyforce.WriteToHDF5(CASENAME, timeStep);
}

void WriteDistributionsToHdf5(const SizeType timeStep) {
    f.WriteToHDF5(CASENAME, timeStep);

}

void WriteNodePropertyToHdf5(const SizeType timeStep) {
    GeometryProperty.WriteToHDF5(CASENAME, timeStep);
    NodeType.WriteToHDF5(CASENAME, timeStep);
}

const std::string& CaseName() { return CASENAME; }
void SetCaseName(const std::string& caseName) { CASENAME = caseName; }
void setCaseName(const char* caseName) {
    std::string tmp(caseName);
    CASENAME = tmp;
}
int SpaceDim() { return SPACEDIM; }

void DestroyFlowfield() {

    FreeArrayMemory(TAUREF);
    // if steady flow
    // FreeArrayMemory(g_MacroVarsCopy);
    FreeArrayMemory(g_ResidualErrorHandle);
    FreeArrayMemory(g_ResidualError);
    // end if steady flow
    // delete[] halos;
}

Real TotalMeshSize() { return 0; }

Real TimeStep() { return DT; }
const Real* pTimeStep() { return &DT; }
void SetTimeStep(Real dt) { DT = dt; }
const Real* TauRef() { return TAUREF; }

void SetTauRef(const std::vector<Real>& tauRef) {
    const int tauNum = SizeofTau();
    if (tauRef.size() == tauNum) {
        if (nullptr == TAUREF) {
            TAUREF = new Real[tauNum];
        }
        for (int idx = 0; idx < tauNum; idx++) {
            TAUREF[idx] = tauRef[idx];
            // ops_printf("\n tau is %f \n",TAUREF[idx]);
        }
    } else {
        ops_printf("Error! %i taus are required but there are %i!\n", tauNum,
                   tauRef.size());
        assert(tauRef.size() == tauNum);
    }
}

Real GetMaximumResidual(const SizeType checkPeriod) {
    Real maxResError{0};
    Real relResErrorMacroVar{0};
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        relResErrorMacroVar = g_ResidualError[2 * macroVarIdx] /
                              g_ResidualError[2 * macroVarIdx + 1] /
                              (checkPeriod * TimeStep());

        if (maxResError <= relResErrorMacroVar) {
            maxResError = relResErrorMacroVar;
        }
    }
    return maxResError;
}

const std::map<std::string,ops_halo_group>& HaloGroups() { return HALOGROUPS; }

void TransferHalos() {
    if (HALOGROUPS.size() > 0) {
        for (auto haloGroup : HALOGROUPS){
             ops_halo_transfer(haloGroup.second);
        }
    }
}

void TransferHalos(const std::string key) {
    ops_halo_transfer(HALOGROUPS[key]);
}

void TransferHalos(const std::vector<std::string> keys) {
    if (keys.size() > 0) {
        for (auto key : keys) {
            ops_halo_transfer(HALOGROUPS[key]);
        }
    }
}
void DefineBlocks(const std::vector<SizeType>& blockIds,
                  const std::vector<std::string>& blockNames,
                  const std::vector<int>& blockSizes) {
    const SizeType blockNum{blockIds.size()};
    if (blockNum != blockNames.size()) {
        ops_printf(
            "Error! The size of blockIds %i is inconsistent with the size of "
            "blockNames %i!\n",
            blockNum, blockNames.size());
        assert(blockNum == blockNames.size());
    }
    if ((SPACEDIM * blockNum) != blockSizes.size()) {
        ops_printf(
            "Error! The size of blockIds %i is inconsistent with the size of "
            "blockSize %i!\n",
            blockNum, blockSizes.size());
        assert(blockNum == blockSizes.size());
    }
    for (SizeType i = 0; i < blockNum; i++) {
        const SizeType blockId{blockIds.at(i)};
        std::vector<int> blockSize(SPACEDIM);
        for (int j = 0; j < SPACEDIM; j++) {
            blockSize.at(j) = blockSizes.at(i * SPACEDIM + j);
        }
        Block block(blockId, blockNames[i], blockSize);
        BLOCKS.emplace(blockId, block);
    }
}

void DefineBlocks(const std::vector<SizeType>& blockIds,
                  const std::vector<std::string>& blockNames,
                  const std::vector<int>& blockSizes, const Real meshSize,
                  const std::map<SizeType, std::vector<Real>>& startPos) {
    DefineBlocks(blockIds, blockNames, blockSizes);
    const SizeType blockNum{BLOCKS.size()};
    SizeType numBlockStartPos{startPos.size()};
    if (numBlockStartPos == (blockNum)) {
        for (auto& idStartPos : startPos) {
            std::vector<std::vector<Real>> blockCoordinates(SPACEDIM);
            const SizeType id{idStartPos.first};
            const std::vector<Real> blockStartPos{idStartPos.second};
            for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
                const int numOfGridPoints{BLOCKS.at(id).Size().at(coordIndex)};
                blockCoordinates.at(coordIndex).resize(numOfGridPoints);
                for (int nodeIndex = 0; nodeIndex < numOfGridPoints;
                     nodeIndex++) {
                    blockCoordinates.at(coordIndex).at(nodeIndex) =
                        blockStartPos.at(coordIndex) + nodeIndex * meshSize;
                }
            }
            COORDINATES.emplace(id, blockCoordinates);
        }
    } else {
        ops_printf(
            "Error! We expect starting points for %i blocks, but only"
            "received %i blocks!\n",
            blockNum, numBlockStartPos);
        assert(numBlockStartPos == blockNum);
    }
    NodeType.CreateFieldFromScratch(BLOCKS);
    CoordinateXYZ.SetDataDim(SPACEDIM);
    CoordinateXYZ.CreateFieldFromScratch(BLOCKS);
    GeometryProperty.CreateFieldFromScratch(BLOCKS);
}





void PrepareFlowField() {
    ops_printf("The coordinates are assigned!\n");
    for (auto& idBlock: BLOCKS) {
        Block& block{idBlock.second};
        const SizeType blockId{idBlock.first};
        SetBlockGeometryProperty(block);
        ops_printf("The geometry property for Block %i is set!\n", blockId);
        for (int compoId = 0; compoId < NUMCOMPONENTS; compoId++) {
            SetBulkandHaloNodesType(block, compoId);
            ops_printf(
                "The bulk and halo node property are set for Component %i at "
                "Block %i\n",
                compoId, blockId);
        }
        AssignCoordinates(block, COORDINATES.at(blockId));
    }
    SetBoundaryNodeType();
}

void DispResidualError3D(const int iter, const SizeType checkPeriod) {
    ops_printf("##########Residual Error at %i time step##########\n", iter);
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        Real residualError = g_ResidualError[2 * macroVarIdx] /
                             g_ResidualError[2 * macroVarIdx + 1] /
                             (checkPeriod * TimeStep());
        ops_printf("Residual of %s = %.17g\n",
                   MacroVarName()[macroVarIdx].c_str(), residualError);
    }
}