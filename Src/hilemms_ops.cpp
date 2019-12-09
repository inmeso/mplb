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


//TODO It seems a bit fancy to manually tranform two types.
// Do we have better method?
VertexTypes BoundaryTypeToVertexType(BoundaryType type) {
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
        case BoundaryType_EQN:
            vtType = Vertex_NoslipEQN;
            break;

        default:
            ops_printf("\n Current Boundary type is not yet implemented");
    }

    return vtType;
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
//TODO to be moved to the boundary module soon
void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         BoundaryType boundaryType,
                         const std::vector<VariableTypes>& macroVarTypes,
                         const std::vector<Real>& macroVarValues) {
    // Type conversion from BoundaryType to VertexType.
    VertexTypes vtType;
    vtType = BoundaryTypeToVertexType(boundaryType);
    // Set the number of halo point required by the boundary condition
    // So far all boundary conditions are implemented in a way that requires no
    // halos so we leave it as the initial value 1
    // The only difference is the periodic boundary condition which needs same
    // halos as required by the numerical scheme, which is set by the scheme
    // module
    // If necessary, uncomment the sentence below and give a corret number
    // SetBoundaryHaloNum(1);

    // Here we adopt the assumption that a boundary is defined by [\rho,u,v,w,T]
    // in 3D or  [\rho,u,v,T] in 2D. For a kernel function for dealing with
    // a boundary condition, these parameters shall be passed in a fixed order
    // as shown.
    const int numMacroVarTypes{(int)macroVarTypes.size()};
    const int numMacroVarValues{(int)macroVarValues.size()};

    std::vector<Real> macroVarsAtBoundary;
    const int macroVarNumOfCurrentComponent{
        VARIABLECOMPPOS[2 * componentID + 1] -
        VARIABLECOMPPOS[2 * componentID] + 1};
    macroVarsAtBoundary.resize(macroVarNumOfCurrentComponent);

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

void PrepareBoundary() {
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

#ifdef OPS_2D
// mark all solid points inside the circle to be ImmersedSolid
void MarkPtsInsideCircleAsSolid(int blockIndex, Real diameter,
                                std::vector<Real> circlePos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* circlePosition = &circlePos[0];
    ops_par_loop(KerSetEmbeddedCircle, "KerSetEmbeddedCircle",
                 g_Block[blockIndex], SPACEDIM, bulkRng,
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE),
                 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                             LOCALSTENCIL, "double", OPS_READ)
                     ops_arg_gbl(&diameter, 1, "double", OPS_READ),
                 ops_arg_gbl(circlePosition, SPACEDIM, "Real", OPS_READ));
}

void MarkPtsInsideEllipseAsSolid(int blockIndex, Real semiMajorAxes,
                                 Real semiMinorAxes,
                                 std::vector<Real> centerPos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* centerPosition = &centerPos[0];
    ops_par_loop(KerSetEmbeddedEllipse, "KerSetEmbeddedCircle",
                 g_Block[blockIndex], SPACEDIM, bulkRng,
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE),
                 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                             LOCALSTENCIL, "double", OPS_READ),
                 ops_arg_gbl(&semiMajorAxes, 1, "double", OPS_READ),
                 ops_arg_gbl(&semiMinorAxes, 1, "double", OPS_READ),
                 ops_arg_gbl(centerPosition, SPACEDIM, "Real", OPS_READ));
}

// Function to wipe off some solid points that cannot be considered as a good surface point.
void WipeSolidPtsBasedNeigbours() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
        ops_par_loop(KerSweep, "KerSweep", g_Block[blockIndex], SPACEDIM,
                     bulkRng,
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ));
    }
}
void KerSyncGeometryProperty(ACC<int>& geometryProperty,
                             const ACC<int>& nodeType)
    // Function to sync the Geometry property to reflect the modifed solid
    // property
    void UpdateGeometryAfterWiping() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
        ops_par_loop(KerSyncGeometryProperty, "KerSyncGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_RW),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ));
    }
}

// set the correct  geometry property e.g., corner types i.e., mark out the surface points

void MarkSurfacePoints() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
        ops_par_loop(KerSetEmbeddedBodyGeometry, "KerSetEmbeddedBodyGeometry",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 ONEPTLATTICESTENCIL, "int", OPS_READ));
    }
}

// set the boundary type
// int nodeType{ surface };
void SetBoundaryTypeofImmersedBody() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());

        int nodeType{Vertex_EQMDiffuseRefl};
        ops_par_loop(KerSetEmbeddedBodyBoundary, "KerSetEmbeddedBodyBoundary",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_RW),
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ));
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