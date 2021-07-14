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
#include  "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
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
RealFieldGroup MacroVars;
RealFieldGroup MacroVarsCopy;
std::map<int,Real> ResidualError;
std::map<int, ops_reduction> ResidualErrorHandle;
std::map<int, Real>& g_ResidualError() { return ResidualError; };
std::map<int, ops_reduction>& g_ResidualErrorHandle() {
    return ResidualErrorHandle;
};

RealFieldGroup MacroBodyforce;
const BlockGroup& g_Block() { return BLOCKS; };
RealField& g_f() { return f; };
RealField& g_fStage() { return fStage; };
RealFieldGroup& g_MacroVars() { return MacroVars; };
RealFieldGroup& g_MacroVarsCopy() { return MacroVarsCopy; };
RealFieldGroup& g_MacroBodyforce() { return MacroBodyforce; };
/**
 * DT: time step
 */
Real DT{1};

RealField CoordinateXYZ{"CoordinateXYZ"};
RealField& g_CoordinateXYZ() { return CoordinateXYZ; };
std::map<SizeType, std::vector<std::vector<Real>>> COORDINATES;
IntFieldGroup NodeType;
IntField GeometryProperty{"GeometryProperty"};
IntFieldGroup& g_NodeType() { return NodeType; };
IntField& g_GeometryProperty() { return GeometryProperty; };

void DefineCase(const std::string& caseName, const int spaceDim,
                const bool transient) {
    if (SPACEDIM != spaceDim) {
        ops_printf("Error! The SPACEDIM here is inconsistent with the\n");
        assert(SPACEDIM == spaceDim);
    }
    CASENAME = caseName;
    TRANSIENT = transient;
}

bool IsTransient() { return TRANSIENT; }

void Partition() {
    ops_partition((char*)"LBM Solver");
    PrepareFlowField();
}

/*
 * We need a name to specify which file to input
 * To be decided: a single filename or an array of filenames
 */

void WriteFlowfieldToHdf5(const SizeType timeStep) {
    for (const auto& macroVar : MacroVars) {
        macroVar.second.WriteToHDF5(CASENAME, timeStep);
    }
    CoordinateXYZ.WriteToHDF5(CASENAME, timeStep);
    for (const auto& force : MacroBodyforce) {
        force.second.WriteToHDF5(CASENAME, timeStep);
    }
}

void WriteDistributionsToHdf5(const SizeType timeStep) {
    f.WriteToHDF5(CASENAME, timeStep);
}

void WriteNodePropertyToHdf5(const SizeType timeStep) {
    GeometryProperty.WriteToHDF5(CASENAME, timeStep);
    for (const auto& pair : NodeType) {
        pair.second.WriteToHDF5(CASENAME, timeStep);
    }
}

const std::string& CaseName() { return CASENAME; }


Real TotalMeshSize() { return 0; }

Real TimeStep() { return DT; }
const Real* pTimeStep() { return &DT; }
void SetTimeStep(Real dt) { DT = dt; }

Real GetMaximumResidual(const SizeType checkPeriod) {
    Real maxResError{0};
    Real relResErrorMacroVar{0};
    for (const auto& error: ResidualError) {
        relResErrorMacroVar = error.second / (checkPeriod * TimeStep());
        if (maxResError <= relResErrorMacroVar) {
            maxResError = relResErrorMacroVar;
        }
    }
    return maxResError;
}

void TransferHalos() {
    fStage.TransferHalos();
}

void DefineBlocks(const std::vector<int>& blockIds,
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
    for (int i = 0; i < blockNum; i++) {
        const int blockId{blockIds.at(i)};
        std::vector<int> blockSize(SPACEDIM);
        for (int j = 0; j < SPACEDIM; j++) {
            blockSize.at(j) = blockSizes.at(i * SPACEDIM + j);
        }
        Block block(blockId, blockNames[i], blockSize);
        BLOCKS.emplace(blockId, block);
    }
}

void DefineBlocks(const std::vector<int>& blockIds,
                  const std::vector<std::string>& blockNames,
                  const std::vector<int>& blockSizes, const Real meshSize,
                  const std::map<int, std::vector<Real>>& startPos) {
    DefineBlocks(blockIds, blockNames, blockSizes);
    const SizeType blockNum{BLOCKS.size()};
    SizeType numBlockStartPos{startPos.size()};
    if (numBlockStartPos == (blockNum)) {
        for (const auto& idStartPos : startPos) {
            std::vector<std::vector<Real>> blockCoordinates(SPACEDIM);
            const int id{idStartPos.first};
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
    //NodeType.CreateFieldFromScratch(BLOCKS);
    CoordinateXYZ.SetDataDim(SPACEDIM);
    CoordinateXYZ.CreateFieldFromScratch(BLOCKS);
    GeometryProperty.CreateFieldFromScratch(BLOCKS);
}
void PrepareFlowField() {
    ops_printf("The coordinates are assigned!\n");
    for (const auto& idBlock: BLOCKS) {
        const Block& block{idBlock.second};
        const int blockId{idBlock.first};
        SetBlockGeometryProperty(block);
        ops_printf("The geometry property for Block %i is set!\n", blockId);
        for (const auto& idCompo:g_Components()) {
            SetBulkandHaloNodesType(block, idCompo.first);
            ops_printf(
                "The bulk and halo node property are set for Component %i at "
                "Block %i\n",
                idCompo.first, blockId);
        }
        AssignCoordinates(block, COORDINATES.at(blockId));
    }
    SetBoundaryNodeType();
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
}

void DispResidualError(const int iter, const SizeType checkPeriod) {
    ops_printf("##########Residual Error at %i time step##########\n", iter);
    for (auto& compo : g_Components()) {
        for (auto& macroVar : compo.second.macroVars) {
            Real residualError = ResidualError.at(macroVar.second.id) /
                                 (checkPeriod * TimeStep());
            ops_printf("Residual of %s = %.17g\n", macroVar.second.name.c_str(),
                       residualError);
        }
    }
}

void DefineBlockConnection(const std::vector<int>& fromBlock,
                           const std::vector<BoundarySurface>& fromSurface,
                           const std::vector<int>& toBlock,
                           const std::vector<BoundarySurface>& toSurface,
                           const std::vector<VertexType>& connectionType) {
    const int fromBlockSize{static_cast<int>(fromBlock.size())};
    const int fromSurfaceSize{static_cast<int>(fromSurface.size())};
    const int toBlockSize{static_cast<int>(toBlock.size())};
    const int toSurfaceSize{static_cast<int>(toSurface.size())};
    const int connectionTypeSize(static_cast<int>(connectionType.size()));
    if ((fromBlockSize != fromSurfaceSize) || (fromBlockSize != toBlockSize) ||
        (fromBlockSize != toSurfaceSize) ||
        (fromSurfaceSize != toSurfaceSize) ||
        (fromSurfaceSize != toBlockSize) || (toBlockSize != toSurfaceSize) ||
        connectionTypeSize != fromBlockSize) {
        ops_printf("Please input consistent halo pairs!\n");
        assert(false);
    }
    for (int idx = 0; idx < fromBlockSize; idx++) {
        Neighbor neighbor;
        neighbor.blockId = toBlock.at(idx);
        neighbor.surface = toSurface.at(idx);
        neighbor.type = connectionType.at(idx);
        BLOCKS.at(fromBlock.at(idx)).AddNeighbor(fromSurface.at(idx), neighbor);
    }
}
