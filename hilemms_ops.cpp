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

#include "hilemms.h"
#include "hilemms_ops_kernel.h"


Real* VERTEXCOORDINATES{nullptr};
int NUMVERTICES{0};

// Structure to hold the values whenever user specifies a boundary condition.
struct BlockBoundary {
    int blockIndex;
    int componentID;
    std::vector<Real> givenVars;
    BoundarySurface boundarySurface;
    VertexTypes boundaryType;
};

// Vector to assemble all boundary conditions so as to use
// in TreatDomainBoundary().
std::vector<BlockBoundary> blockBoundaryConditions;

// This routine finds the index range of boundary surface
int* BoundarySurfaceRange(const int blockId, BoundarySurface surface) {
    int* boundarySurfaceRange;
    switch (surface) {
        case BoundarySurface_Left:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngImin());
            break;

        case BoundarySurface_Right:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngImax());
            break;

        case BoundarySurface_Top:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngJmax());
            break;

        case BoundarySurface_Bottom:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngJmin());
            break;

        case BoundarySurface_Front:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngKmax());
            break;

        case BoundarySurface_Back:
            boundarySurfaceRange = BlockIterRng(blockId, IterRngKmin());
            break;

        default:
            ops_printf("\n Surface entered for the BC is incorrect.");
    }

    return boundarySurfaceRange;
}

void SetBulkandHaloNodesType(int blockIndex, int compoId) {
    int nodeType = (int)Vertex_Fluid;
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    // specify halo points
    nodeType = (int)Vertex_ImmersedSolid;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngImin());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngImax());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));

        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    }
    FreeArrayMemory(haloIterRng);
}

void AssignCoordinates(int blockIndex,
                       const std::vector<std::vector<Real>>& coordinates) {
#ifdef OPS_2D
    if (SPACEDIM == 2) {
        int* range = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerSetCoordinates, "KerSetCoordinates",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_gbl(coordinates[0].data(),
                                 BlockSize(blockIndex)[0], "double", OPS_READ),
                     ops_arg_gbl(coordinates[1].data(),
                                 BlockSize(blockIndex)[1], "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE));
    }
#endif

#ifdef OPS_3D
    if (SPACEDIM == 3) {
        int* range = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerSetCoordinates3D, "KerSetCoordinates3D",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_gbl(coordinates[0].data(),
                                 BlockSize(blockIndex)[0], "double", OPS_READ),
                     ops_arg_gbl(coordinates[1].data(),
                                 BlockSize(blockIndex)[1], "double", OPS_READ),
                     ops_arg_gbl(coordinates[2].data(),
                                 BlockSize(blockIndex)[2], "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE));
    }
#endif
}

void CalculateBlockCoordinates(const int blockIndex, Real* blockStartPos,
                               Real meshSize) {
    std::vector<std::vector<Real>> coordinates(SPACEDIM);
    for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
        int numOfGridPoints =
            BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex];
        coordinates[coordIndex].resize(numOfGridPoints);
        for (int nodeIndex = 0; nodeIndex < numOfGridPoints; nodeIndex++) {
            coordinates[coordIndex][nodeIndex] =
                blockStartPos[coordIndex] + nodeIndex * meshSize;
        }
    }
    AssignCoordinates(blockIndex, coordinates);
}

// This subroutine is for internal use only.
VertexTypes BoundTypeToVertexType(BoundaryType type) {
    VertexTypes vtType;

    switch (type) {
        case BoundaryType_KineticDiffuseWall:
            vtType = Vertex_KineticDiffuseWall;
            break;

        case BoundaryType_KineticSpelluarWall:
            vtType = Vertex_KineticSpelluarWall;
            break;

        case BoundaryType_SlipWall:
            vtType = Vertex_SlipWall;
            break;

        case BoundaryType_VelocityInlet:
            vtType = Vertex_VelocityInlet;
            break;

        case BoundaryType_VelocityOutlet:
            vtType = Vertex_VelocityOutlet;
            break;

        case BoundaryType_ExtrapolPressure1ST:
            vtType = Vertex_ExtrapolPressure1ST;
            break;

        case BoundaryType_ExtrapolPressure2ND:
            vtType = Vertex_ExtrapolPressure2ND;
            break;

        case BoundaryType_Periodic:
            vtType = Vertex_Periodic;
            break;

        case BoundaryType_Uniform:
            vtType = Vertex_Uniform;
            break;

        case BoundaryType_BounceBackWall:
            vtType = Vertex_BounceBackWall;
            break;

        case BoundaryType_FreeFlux:
            vtType = Vertex_FreeFlux;
            break;

        case BoundaryType_ZouHeVelocity:
            vtType = Vertex_ZouHeVelocity;
            break;

        case BoundaryType_EQMDiffuseRefl:
            vtType = Vertex_EQMDiffuseRefl;
            break;

        default:
            ops_printf("\n Current Boundary type is not yet implemented");
    }

    return vtType;
}

//TODO It may be better to add one more init module
void DefineProblemDomain(const int blockNum, const std::vector<int> blockSize,
                         const Real meshSize,
                         const std::vector<Real> startPos) {
    SetBlockNum(blockNum);
    SetBlockSize(blockSize);
    DefineVariables();
// TODO We need to define the halo relation and
#ifdef OPS_3D
    DefineHaloTransfer3D();
#endif  // OPS_3D

#ifdef OPS_2D
    DefineHaloTransfer();
#endif

    ops_partition((char*)"LBM Solver");
    ops_printf("%i blocks are parted and all field variabls allocated!\n",
               BlockNum());
    int numBlockStartPos;
    numBlockStartPos = startPos.size();

    if (numBlockStartPos == blockNum * SPACEDIM) {
        for (int blockIndex = 0; blockIndex < blockNum; blockIndex++) {
            // One block will have 3 values as starting position in x, y, z
            // direction respectively.
            Real* blockStartPosition = new Real[SPACEDIM];

            for (int spaceDim = 0; spaceDim < SPACEDIM; spaceDim++) {
                blockStartPosition[spaceDim] =
                    startPos[blockIndex * SPACEDIM + spaceDim];
            }

            CalculateBlockCoordinates(blockIndex, blockStartPosition, meshSize);
            delete[] blockStartPosition;
        }
    } else {
        ops_printf(
            "Error! Expected %i coordinates of three starting points %i, but "
            "received only =%i \n",
            SPACEDIM * blockNum, numBlockStartPos);
        assert(numBlockStartPos == blockNum * SPACEDIM);
    }
    ops_printf("The coordinates are assigned!\n");
    for (int blockId = 0; blockId < blockNum; blockId++) {
        SetBlockGeometryProperty(blockId);
        ops_printf("The geometry property for Block %i is set!\n", blockId);
        for (int compoId = 0; compoId < NUMCOMPONENTS; compoId++) {
            SetBulkandHaloNodesType(blockId, compoId);
            ops_printf(
                "The bulk and halo node property are set for Component %i at "
                "Block %i\n",
                compoId, blockId);
        }
    }

    for (int bcIdx = 0; bcIdx < blockBoundaryConditions.size(); bcIdx++) {
        const int compoId{blockBoundaryConditions[bcIdx].componentID};
        const int blockIndex{blockBoundaryConditions[bcIdx].blockIndex};
        int* bcRange = BoundarySurfaceRange(
            blockBoundaryConditions[bcIdx].blockIndex,
            blockBoundaryConditions[bcIdx].boundarySurface);
        const int vtType{(int)blockBoundaryConditions[bcIdx].boundaryType};
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, bcRange,
                     ops_arg_gbl(&vtType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    }
}

// Check whether this needs to be defines using OPS Kernel.
Real GetMaximumResidualError(const Real checkPeriod) {
    Real maxResError = 1E-15;
    Real relResErrorMacroVar;
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

void Iterate(const int steps, const int checkPointPeriod) {
    const SchemeType scheme = Scheme();
    ops_printf("Starting the iteration...\n");
    switch (scheme) {
        case Scheme_StreamCollision: {
            for (int iter = 0; iter < steps; iter++) {
#ifdef OPS_3D
                StreamCollision3D();  // Stream-Collision scheme
                // TimeMarching();//Finite difference scheme + cutting cell
                if ((iter % checkPointPeriod) == 0 && iter != 0) {
                    UpdateMacroVars3D();
                    CalcResidualError3D();
                    DispResidualError3D(iter, checkPointPeriod * TimeStep());
                    WriteFlowfieldToHdf5(iter);
                    WriteDistributionsToHdf5(iter);
                    WriteNodePropertyToHdf5(iter);
                }
#endif  // end of OPS_3D
#ifdef OPS_2D
                StreamCollision();  // Stream-Collision scheme
                // TimeMarching();//Finite difference scheme + cutting cell
                if ((iter % checkPointPeriod) == 0 && iter != 0) {
                    UpdateMacroVars();
                    CalcResidualError();
                    DispResidualError(iter, checkPointPeriod * TimeStep());
                    WriteFlowfieldToHdf5(iter);
                    WriteDistributionsToHdf5(iter);
                    WriteNodePropertyToHdf5(iter);
                }
#endif  // end of OPS_2D
            }
        } break;
        default:
            break;
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
    DestroyFlowfield();
}

void Iterate(const Real convergenceCriteria, const int checkPointPeriod) {
    const SchemeType scheme = Scheme();
    ops_printf("Starting the iteration...\n");
    switch (scheme) {
        case Scheme_StreamCollision: {
            int iter{0};
            Real residualError{1};
            do {
#ifdef OPS_3D
                StreamCollision3D();  // Stream-Collision scheme
                if ((iter % checkPointPeriod) == 0) {
                    UpdateMacroVars3D();
                    CalcResidualError3D();
                    residualError =
                        GetMaximumResidualError(checkPointPeriod * TimeStep());
                    DispResidualError3D(iter, checkPointPeriod * TimeStep());
                    WriteFlowfieldToHdf5(iter);
                    WriteDistributionsToHdf5(iter);
                    WriteNodePropertyToHdf5(iter);
                }
#endif  // end of OPS_3D

#ifdef OPS_2D
                StreamCollision();  // Stream-Collision scheme
                // TimeMarching();//Finite difference scheme + cutting cell
                if ((iter % checkPointPeriod) == 0 && iter != 0) {
                    UpdateMacroVars();
                    CalcResidualError();
                    residualError =
                        GetMaximumResidualError(checkPointPeriod * TimeStep());
                    DispResidualError(iter, checkPointPeriod * TimeStep());
                    WriteFlowfieldToHdf5(iter);
                    WriteDistributionsToHdf5(iter);
                    WriteNodePropertyToHdf5(iter);
                }

#endif  // end of OPS_2D
                iter = iter + 1;
            } while (residualError >= convergenceCriteria);
        } break;
        default:
            break;
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
    DestroyFlowfield();
}

void AllocateVertices(const int vertexNum) {
    if (vertexNum == NUMVERTICES) {
        if (nullptr == VERTEXCOORDINATES) {
            VERTEXCOORDINATES = new Real[SPACEDIM * vertexNum];
        }
    }
}

void AddEmbeddedBody(int vertexNum, Real* vertexCoords) {
    NUMVERTICES = vertexNum;
    AllocateVertices(vertexNum);

    int numberVertexCoords;
    numberVertexCoords = sizeof(vertexCoords) / sizeof(vertexCoords[0]);

    if (numberVertexCoords == vertexNum * SPACEDIM) {
#ifdef OPS_2D
        for (int i = 0; i < SPACEDIM * vertexNum; i = i + 2) {
            VERTEXCOORDINATES[i] = vertexCoords[i];          // x_coordinate
            VERTEXCOORDINATES[i + 1] = vertexCoords[i + 1];  // y_coordinate
        }
#endif  // OPS_2D

#ifdef OPS_3D
        for (int i = 0; i < SPACEDIM * vertexNum; i = i + 3) {
            VERTEXCOORDINATES[i] = vertexCoords[i];          // x_coordinate
            VERTEXCOORDINATES[i + 1] = vertexCoords[i + 1];  // y_coordinate
            VERTEXCOORDINATES[i + 2] = vertexCoords[i + 2];  // z_coordinate
        }
#endif  // OPS_3D
    } else {
        ops_printf(
            " For %i dimensional problem, number of vertices should be %i "
            "but received only %i \n",
            SPACEDIM, vertexNum * SPACEDIM, numberVertexCoords);
    }

    ops_decl_const("NUMVERTICES", 1, "int", &NUMVERTICES);
    ops_decl_const("VERTEXCOORDINATES", SPACEDIM * vertexNum, "int",
                   VERTEXCOORDINATES);
}

void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         BoundaryType boundaryType,
                         const std::vector<VariableTypes>& macroVarTypes,
                         const std::vector<Real>& macroVarValues) {
    // Type conversion from BoundaryType to VertexType.
    VertexTypes vtType;
    vtType = BoundTypeToVertexType(boundaryType);
    // Set the number of halo point required by the boundary condition
    // So far all boundary conditions are implemented in a way that requires no
    // halos so we leave it as the initial value 1
    // The only difference is the periodic boundary condition which needs same
    // halos as required by the numerical scheme, which is set by the scheme
    // module
    // If necessary, uncomment the sentence below and give a corret number
    // SetBoundaryHaloNum(1);

    // Here we adopt the assumtion that a boundary is defined by [\rho,u,v,w,T]
    // in 3D or  [\rho,u,v,T] in 2D. For a kernel function for dealing with
    // a boundary condition, these parameters shall be passed in a fixed order
    // as shown.
    const int numMacroVarTypes{(int)macroVarTypes.size()};
    const int numMacroVarValues{(int)macroVarValues.size()};
    // TODO The logic may need rethink
    std::vector<Real> macroVarsAtBoundary;
    if (2 == SPACEDIM) {
        macroVarsAtBoundary.resize(4);
    }
    if (3 == SPACEDIM) {
        macroVarsAtBoundary.resize(5);
    }

    if (numMacroVarTypes == numMacroVarValues) {
        for (int i = 0; i < numMacroVarValues; i++) {
            int varPos{0};
            switch (macroVarTypes[i]) {
                case Variable_Rho:
                    varPos = 0;
                    break;
                case Variable_U:
                    varPos = 1;
                    break;
                case Variable_V:
                    varPos = 2;
                    break;
                case Variable_W:
                    if (3 == SPACEDIM) {
                        varPos = 3;
                    } else {
                        varPos = -1;
                        ops_printf(
                            "Error! The velocity component w is defined/used "
                            "for %iD problem.\n",
                            SPACEDIM);
                        assert(3 == SPACEDIM);
                    }
                    break;
                case Variable_U_Force:
                    varPos = 1;
                    break;
                case Variable_V_Force:
                    varPos = 2;
                    break;
                case Variable_W_Force:
                    if (3 == SPACEDIM) {
                        varPos = 3;
                    } else {
                        varPos = -1;
                        ops_printf(
                            "Error! The velocity component w is defined/used "
                            "for %iD problem.\n",
                            SPACEDIM);
                        assert(3 == SPACEDIM);
                    }
                    break;
                case Variable_T:
                    if (3 == SPACEDIM) {
                        varPos = 4;
                    } else {
                        varPos = 2;
                    }
                    break;
                default:
                    break;
            }
            macroVarsAtBoundary[varPos] = macroVarValues[i];
        }
        BlockBoundary domainBoundaryCondition;
        domainBoundaryCondition.blockIndex = blockIndex;
        domainBoundaryCondition.componentID = componentID;
        domainBoundaryCondition.givenVars = macroVarsAtBoundary;
        domainBoundaryCondition.boundarySurface = boundarySurface;
        domainBoundaryCondition.boundaryType = vtType;
        blockBoundaryConditions.push_back(domainBoundaryCondition);
        ops_printf(
            "The boundary condition %i is adopted for Component %i at Surface "
            "%i, Block %i.\n",
            domainBoundaryCondition.boundaryType,
            domainBoundaryCondition.componentID,
            domainBoundaryCondition.boundarySurface,
            domainBoundaryCondition.blockIndex);
    } else {
        ops_printf("Error! Expected %i values for BC but received only %i \n",
                   numMacroVarTypes, numMacroVarValues);
        assert(numMacroVarTypes == numMacroVarValues);
    }
}

void ImplementBoundaryConditions() {
    int totalNumBoundCond;
    totalNumBoundCond = blockBoundaryConditions.size();
    if (totalNumBoundCond != 0) {
        for (int i = 0; i < totalNumBoundCond; i++) {
            int* rangeBoundaryCondition;
            rangeBoundaryCondition = BoundarySurfaceRange(
                blockBoundaryConditions[i].blockIndex, blockBoundaryConditions[i].boundarySurface);
#ifdef OPS_2D

            TreatDomainBoundary(blockBoundaryConditions[i].blockIndex,
                                blockBoundaryConditions[i].componentID,
                                blockBoundaryConditions[i].givenVars.data(),
                                rangeBoundaryCondition,
                                blockBoundaryConditions[i].boundaryType);
#endif  // End of OPS_2D

#ifdef OPS_3D
            TreatBlockBoundary3D(blockBoundaryConditions[i].blockIndex,
                                 blockBoundaryConditions[i].componentID,
                                 blockBoundaryConditions[i].givenVars.data(),
                                 rangeBoundaryCondition,
                                 blockBoundaryConditions[i].boundaryType);
#endif  // End of OPS_3D
        }

    }

    else {
        ops_printf("\n No Boundary condition has been defined.");
    }
}

void InitialiseNodeMacroVars(Real* nodeMacroVars, const Real* nodeCoordinates) {
#ifdef OPS_2D
    Real x{nodeCoordinates[0]};
    Real y{nodeCoordinates[1]};
    nodeMacroVars[0] = 1;        // rho
    nodeMacroVars[1] = 0;        // u
    nodeMacroVars[2] = 0;        // v
#endif

#ifdef OPS_3D
    Real x{nodeCoordinates[0]};
    Real y{nodeCoordinates[1]};
    Real z{nodeCoordinates[2]};  // for 3D problems
    nodeMacroVars[0] = 1;        // rho
    nodeMacroVars[1] = 0;        // u
    nodeMacroVars[2] = 0;        // v
    nodeMacroVars[3] = 0;        // w
#endif
}

void DefineInitialCondition() {
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        void KerSetInitialMacroVars(Real * macroVars, const Real* coordinates,
                                    const int* idx);
        int* iterRng = BlockIterRng(blockIdx, IterRngWhole());
        ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                     g_Block[blockIdx], SPACEDIM, iterRng,
                     ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW),
                     ops_arg_dat(g_CoordinateXYZ[blockIdx], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_idx());
    }
    ops_printf("Macroscopic variables are initialised!\n");
#ifdef OPS_3D
    InitialiseSolution3D();
#endif
#ifdef OPS_2D
    InitialiseSolution();
#endif
    ops_printf("Distribution functions are initialised\n");
}

void  SetBlockGeometryProperty(int blockIndex) {
    int geometryProperty = (int)VG_Fluid;
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    // specify halo points
    geometryProperty = VG_ImmersedSolid;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));

    iterRange = BlockIterRng(blockIndex, IterRngImax());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    // specify domain
    geometryProperty = VG_JP;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_JM;
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IP;
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IM;
    iterRange = BlockIterRng(blockIndex, IterRngImax());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        geometryProperty = VG_KP;
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        geometryProperty = VG_KM;
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    const int nx = BlockSize(blockIndex)[0];
    const int ny = BlockSize(blockIndex)[1];
    // 2D Domain corner points four types
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    if (3 == SPACEDIM) {
        const int nz = BlockSize(blockIndex)[2];
        // 3D Domain edges 12 types
        int iminjmin[]{0, 1, 0, 1, 0, nz};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int iminkmin[]{0, 1, 0, ny, 0, 1};
        geometryProperty = VG_IPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        geometryProperty = VG_IPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        geometryProperty = VG_IMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        geometryProperty = VG_IMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int jminkmin[]{0, nx, 0, 1, 0, 1};
        geometryProperty = VG_JPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_JPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_JMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_JMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // 3D domain corners 8 types
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        geometryProperty = VG_IPJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        geometryProperty = VG_IPJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        geometryProperty = VG_IPJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IPJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        geometryProperty = VG_IMJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_IMJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_IMJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IMJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }
}

// Define varios halo numbers such as HaloNum, HaloDepth and SchemeHaloPt.
void DefineHaloNumber(int Halo_Number, int Halo_Depth, int Scheme_Halo_points,
                      int Num_Bound_Halo_Points) {
    // g_HaloNum   = Halo_Number;
    // g_HaloDepth = Halo_Depth;
    // schemeHaloPt = Scheme_Halo_points;

    SetHaloRelationNum(Halo_Number);
    SetHaloDepth(Halo_Depth);
    SetSchemeHaloNum(Scheme_Halo_points);

    // boundaryHaloPt = Num_Bound_Halo_Points;
    SetBoundaryHaloNum(Num_Bound_Halo_Points);
}

#ifdef OPS_2D
// mark all solid points inside the circle to be ImmersedSolid
void MarkPtsInsideCircleAsSolid(int blockIndex, Real diameter,
                             std::vector<Real> circlePos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* circlePosition = &circlePos[0];
    ops_par_loop(KerSetEmbeddedCircle, "KerSetEmbeddedCircle",
                 g_Block[blockIndex], SPACEDIM, bulkRng,
                 ops_arg_gbl(&diameter, 1, "double", OPS_READ),
                 ops_arg_gbl(circlePosition, SPACEDIM, "Real", OPS_READ),
                 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                             LOCALSTENCIL, "double", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
}

void MarkPtsInsideEllipseAsSolid(int blockIndex, Real semiMajorAxes,
                              Real semiMinorAxes, std::vector<Real> centerPos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* centerPosition = &centerPos[0];
    ops_par_loop(
        KerSetEmbeddedEllipse, "KerSetEmbeddedCircle", g_Block[blockIndex],
        SPACEDIM, bulkRng, ops_arg_gbl(&semiMajorAxes, 1, "double", OPS_READ),
        ops_arg_gbl(&semiMinorAxes, 1, "double", OPS_READ),
        ops_arg_gbl(centerPosition, SPACEDIM, "Real", OPS_READ),
        ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL,
                    "double", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS, LOCALSTENCIL, "int", OPS_WRITE),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_WRITE));
}

// Wrapper function for embedded body.
void HandleImmersedSolid() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
        // wipe off some solid points that cannot be considered
        // as a good surface point
        ops_par_loop(KerSweep, "KerSweep", g_Block[blockIndex], SPACEDIM,
                     bulkRng,
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // sync the Geometry property to reflect the modifed solid property
        ops_par_loop(KerSyncGeometryProperty, "KerSyncGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_RW));

        // set the correct  geometry property e.g., corner types
        // i.e., mark out the surface points
        ops_par_loop(KerSetEmbeddedBodyGeometry, "KerSetEmbeddedBodyGeometry",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 ONEPTLATTICESTENCIL, "int", OPS_RW),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // set the boundary type
        // int nodeType{ surface };
        int nodeType{Vertex_EQMDiffuseRefl};
        ops_par_loop(KerSetEmbeddedBodyBoundary, "KerSetEmbeddedBodyBoundary",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_RW));
    }
}

// Function to provide details of embedded solid body into the fluid.
void AddEmbeddedBody(SolidBodyType type, int blockIndex,
                  std::vector<Real> centerPos, std::vector<Real> controlParas) {
    int numCoordCenterPos;
    numCoordCenterPos = centerPos.size();

    if (numCoordCenterPos == SPACEDIM) {
        switch (type) {
            case SolidBody_circle: {
                MarkPtsInsideCircleAsSolid(blockIndex, controlParas[0], centerPos);
                break;
            }

            case SolidBody_ellipse: {
                Real semiMajorAxes{controlParas[0]};
                Real semiMinorAxes{controlParas[1]};
                MarkPtsInsideEllipseAsSolid(blockIndex, semiMajorAxes,
                                         semiMinorAxes, centerPos);
                break;
            }

            default:
                ops_printf(
                    "\n This solid body is not yet implemented in the "
                    "code");
                break;
        }
    } else {
        ops_printf(
            "\n For %i dimensional problem, number of coordinates should be "
            "%d, however %d were provided.",
            SPACEDIM, SPACEDIM, numCoordCenterPos);
    }
}
#endif