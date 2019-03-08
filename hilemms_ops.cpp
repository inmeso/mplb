#include "evolution.h"
#include "evolution3d.h"
#include "flowfield.h"
#include "hilemms.h"
#include "hilemms_ops_kernel.h"
#include "model.h"
#include "scheme.h"
#include "type.h"

//#include "setup_comput_domain.h"

int MAXITER;
int CHECKPERIOD;

Real* VERTEXCOORDINATES{nullptr};
int NUMVERTICES{0};

// void TreatDomainBoundary(const int blockIndex, const int componentID,
//                          const Real* givenVars, int* range,
//                          const VertexTypes boundaryType)

// Strcuture to hold the values whenever user specifies a boundary condition.
struct DomainBoundary {
    int blockIndex;
    int componentID;
    Real* givenVars;
    int* range;
    VertexTypes boundaryType;
};

// Vector to assemble all boundary conditions so as to use
// in TreatDomainBoundary().
vector<DomainBoundary> g_domainBoundCond;

void DefineCase(std::string caseName, const int spaceDim) {
    SetCaseName(caseName);
    SPACEDIM = spaceDim;
}


// Function copied from Setup_comput_domain.cpp (3D File)
//void AssignCoordinates(int blockIndex, Real* coordinates[SPACEDIM]) {
void AssignCoordinates(int blockIndex, Real** coordinates) {
#ifdef OPS_2D
    if (SPACEDIM == 2) {
        int* range = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerSetCoordinates, "KerSetCoordinates",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_gbl(coordinates[0], BlockSize(blockIndex)[0],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[1], BlockSize(blockIndex)[1],
                                 "double", OPS_READ),
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
                     ops_arg_gbl(coordinates[0], BlockSize(blockIndex)[0],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[1], BlockSize(blockIndex)[1],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[2], BlockSize(blockIndex)[2],
                                 "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE));
    }
#endif
}

void CalBlockCoordinates(const int blockIndex, Real* blockStartPos,
                         Real meshSize) {
    Real* coordinates[SPACEDIM];

    for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
        int numCellsOneDir =
            BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex];

        coordinates[coordIndex] = new Real[numCellsOneDir + 1];

        coordinates[coordIndex][0] = blockStartPos[coordIndex];

        for (int nodeIndex = 0; nodeIndex < numCellsOneDir; nodeIndex++) {
            coordinates[coordIndex][nodeIndex + 1] =
                coordinates[coordIndex][nodeIndex] + meshSize;
        }
    }

    AssignCoordinates(blockIndex, coordinates);
    //AssignCoordinates(blockIndex, &coordinates[SPACEDIM]);

    // Deleting the memeory allocated to coordinates.
    for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
        if (coordinates[coordIndex] != nullptr)
            delete[] coordinates[coordIndex];
    }
}


// This subroutine is for internal use only.
VertexTypes BoundTypeToVertexType(BoundaryType type) {
    VertexTypes vtType;

    switch (type) {
        case BoundType_KineticDiffuseWall:
            vtType = Vertex_KineticDiffuseWall;
            break;

        case BoundType_KineticSpelluarWall:
            vtType = Vertex_KineticSpelluarWall;
            break;

        case BoundType_SlipWall:
            vtType = Vertex_SlipWall;
            break;

        case BoundType_VelocityInlet:
            vtType = Vertex_VelocityInlet;
            break;

        case BoundType_VelocityOutlet:
            vtType = Vertex_VelocityOutlet;
            break;

        case BoundType_ExtrapolPressure1ST:
            vtType = Vertex_ExtrapolPressure1ST;
            break;

        case BoundType_ExtrapolPressure2ND:
            vtType = Vertex_ExtrapolPressure2ND;
            break;

        case BoundType_Periodic:
            vtType = Vertex_Periodic;
            break;

        case BoundType_Uniform:
            vtType = Vertex_Uniform;
            break;

        case BoundType_BounceBackWall:
            vtType = Vertex_BounceBackWall;
            break;

        case BoundType_FreeFlux:
            vtType = Vertex_FreeFlux;
            break;

        case BoundType_ZouHeVelocity:
            vtType = Vertex_ZouHeVelocity;
            break;

        case BoundType_NoneqExtrapol:
            vtType = Vertex_NoneqExtrapol;
            break;

        case BoundType_EQMDiffuseRefl:
            vtType = Vertex_EQMDiffuseRefl;
            break;

        case BoundType_NonEqExtrapolPressure:
            vtType = Vertex_NonEqExtrapolPressure;
            break;

        default:
            ops_printf("\n Current Boundary type is not yet implemented");
    }

    return vtType;
}

// This routine is for both 3D and 2D.
void SetupGeomPropAndNodeType(int blockIndex, BoundaryType* boundType) {
#ifdef OPS_3D

    VertexTypes faceType[6];

    // Get Vertx type information from boundary Type. This is because MPLB code
    // uses vertex type information.
    for (int i = 0; i < 6; i++) {
        faceType[i] = BoundTypeToVertexType(boundType[i]);
    }

#endif  // end of OPS_3D

#ifdef OPS_2D

    VertexTypes faceType[4];

    // Get Vertx type information from boundary Type. This is because MPLB code
    // uses vertex type information.
    for (int i = 0; i < 4; i++) {
        faceType[i] = BoundTypeToVertexType(boundType[i]);
    }

#endif  // end of OPS_2D

    // This function assocaites various ranges such as imin, imax etc for a
    // block.
    SetupDomainGeometryProperty(blockIndex);

    SetupDomainNodeType(blockIndex, faceType);
}

void DefineProblemDomain(const int blockNum, const std::vector<int> blockSize,
                         const Real meshSize,
                         const std::vector<Real> startPos) {
    SetBlockNum(blockNum);
    SetBlockSize(blockSize);

    DefineVariables();
    //ops_printf("Variables Defined \n");
    ops_partition((char*)"LBMPreProcessor");

    int numBlockStartPos;
    numBlockStartPos = startPos.size();

    if (numBlockStartPos == blockNum * SPACEDIM) {
        for (int blockIndex = 0; blockIndex < blockNum; blockIndex++) {
            Real* blockStartPosition;  // One block will have 3 values as
                                       // statring position in x, y, x direction
                                       // respectively.
            blockStartPosition = new Real[SPACEDIM];

            for (int spaceDim = 0; spaceDim < SPACEDIM; spaceDim++) {
                blockStartPosition[spaceDim] =
                    startPos[blockIndex * SPACEDIM + spaceDim];
            }

            CalBlockCoordinates(blockIndex, blockStartPosition, meshSize);
            delete[] blockStartPosition;
        }
    } else {
        ops_printf(
            "\n Expected number of starting positions are = %i, but received "
            "only =%i \n",
            SPACEDIM * blockNum, numBlockStartPos);
    }
}

void Iterate(SchemeType scheme, const int steps, const int checkPointPeriod) {
    MAXITER = steps;
    CHECKPERIOD = checkPointPeriod;

    for (int iter = 0; iter < MAXITER; iter++) {
#ifdef OPS_3D
        StreamCollision3D();  // Stream-Collision scheme
        // TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % CHECKPERIOD) == 0 && iter != 0) {
            //#ifdef debug
            UpdateMacroVars3D();
            CalcResidualError3D();
            DispResidualError3D(iter, CHECKPERIOD * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
            // if ((densityResidualError + uResidualError+vResidualError) <=
            // 1e-12) break;
            // WriteDistributionsToHdf5(iter);
            // WriteNodePropertyToHdf5(iter);
            //#endif
        }
#endif  // end of OPS_3D

#ifdef OPS_2D
        StreamCollision();  // Stream-Collision scheme
        // TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % CHECKPERIOD) == 0 && iter != 0) {
            //#ifdef debug
            UpdateMacroVars();
            CalcResidualError();
            DispResidualError(iter, CHECKPERIOD * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
            // if ((densityResidualError + uResidualError+vResidualError) <=
            // 1e-12) break;
            // WriteDistributionsToHdf5(iter);
            // WriteNodePropertyToHdf5(iter);
            //#endif
        }
#endif  // end of OPS_2D
    }

    DestroyModel();
    DestroyFlowfield();
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

void Iterate(SchemeType scheme, const Real convergenceCriteria,
             const int checkPointPeriod) {
    CHECKPERIOD = checkPointPeriod;
    int iter = 0;
    Real residualError = 10000;  // initially defining to be a very high value.

    while (residualError >= convergenceCriteria) {
#ifdef OPS_3D
        StreamCollision3D();  // Stream-Collision scheme
        // TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % CHECKPERIOD) == 0 && iter != 0) {
            //#ifdef debug
            UpdateMacroVars3D();
            CalcResidualError3D();
            residualError = GetMaximumResidualError(CHECKPERIOD * TimeStep());
            DispResidualError3D(iter, CHECKPERIOD * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
            // if ((densityResidualError + uResidualError+vResidualError) <=
            // 1e-12) break;
            // WriteDistributionsToHdf5(iter);
            // WriteNodePropertyToHdf5(iter);
            //#endif
        }
        iter = iter + 1;  // Required for checkpoint.
#endif                    // end of OPS_3D

#ifdef OPS_2D
        StreamCollision();  // Stream-Collision scheme

        // TimeMarching();//Finite difference scheme + cutting cell
        if ((iter % CHECKPERIOD) == 0 && iter != 0) {
            //#ifdef debug
            UpdateMacroVars();
            CalcResidualError();
            residualError = GetMaximumResidualError(CHECKPERIOD * TimeStep());
            DispResidualError(iter, CHECKPERIOD * TimeStep());
            WriteFlowfieldToHdf5(iter);
            WriteDistributionsToHdf5(iter);
            WriteNodePropertyToHdf5(iter);
            // if ((densityResidualError + uResidualError+vResidualError) <=
            // 1e-12) break;
            // WriteDistributionsToHdf5(iter);
            // WriteNodePropertyToHdf5(iter);
            //#endif
        }

        iter = iter + 1;  // Required for checkpoint.
#endif                    // end of OPS_2D
    }

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

void AddEmbededBody(int vertexNum, Real* vertexCoords) {
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

// This routine is for internal use and gives the cells on which a
// particular BC is to be applied.
int* RangeBoundCond(const int blockId, BoundarySurface surface) {
    int* rangeBoundaryCond;
    switch (surface) {
        case BoundSurf_Inlet:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngImin());
            break;

        case BoundSurf_Outlet:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngImax());
            break;

        case BoundSurf_Top:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngJmax());
            break;

        case BoundSurf_Bottom:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngJmin());
            break;

        case BoundSurf_Front:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngKmax());
            break;

        case BoundSurf_Back:
            rangeBoundaryCond = BlockIterRng(blockId, IterRngKmin());
            break;

        default:
            ops_printf("\n Surface entered for the BC is incorrect.");
    }

    return rangeBoundaryCond;
}

void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface surface, BoundaryType type,
                         std::vector<VariableTypes> macroVarsComp,
                         std::vector<Real> valuesMacroVarsComp) {
    // Find the range of cells on which a BC will act on.
    int* rangeBoundaryCond;
    rangeBoundaryCond = RangeBoundCond(blockIndex, surface);

    // Type conversion from BoundaryType to VertexType.
    VertexTypes vtType;
    vtType = BoundTypeToVertexType(type);

    // A component can only define BC for the macrovars which were used in
    // the DefineMacroVars. The number of BC can be less than or equal to the
    // number of macrovars for the component defined in DefineMacroVars.
    int numBcComponent;
    numBcComponent = VARIABLECOMPPOS[2 * componentID + 1] -
                     VARIABLECOMPPOS[2 * componentID] + 1;

    // Array which will store all maco vars according to predefined order
    // and pass it to treat domain boundary. This way will involve less
    // change in the MPLB code.
    Real* macroVarsBoundCond{nullptr};

    const int numMacroVarsComp{(int)macroVarsComp.size()};
    const int numValuesMacroVarsComp{(int)valuesMacroVarsComp.size()};

    if (numMacroVarsComp == numValuesMacroVarsComp) {
        macroVarsBoundCond = new Real[numBcComponent];

        for (int i = 0; i < numValuesMacroVarsComp; i++) {
            macroVarsBoundCond[(int)macroVarsComp[i]] = valuesMacroVarsComp[i];
        }

        DomainBoundary domainBoundCondition;
        domainBoundCondition.blockIndex = blockIndex;
        domainBoundCondition.componentID = componentID;
        domainBoundCondition.givenVars = macroVarsBoundCond;
        domainBoundCondition.range = rangeBoundaryCond;
        domainBoundCondition.boundaryType = vtType;

        g_domainBoundCond.push_back(domainBoundCondition);

    } else {
        ops_printf("\n Expected %i values for BC but received only %i ",
                   numMacroVarsComp, numValuesMacroVarsComp);
    }
}

void ImplementBoundaryConditions() {
    int totalNumBoundCond;
    totalNumBoundCond = g_domainBoundCond.size();
    if (totalNumBoundCond != 0) {
        for (int i = 0; i < totalNumBoundCond; i++) {
#ifdef OPS_2D

            TreatDomainBoundary(g_domainBoundCond[i].blockIndex,
                                g_domainBoundCond[i].componentID,
                                g_domainBoundCond[i].givenVars,
                                g_domainBoundCond[i].range,
                                g_domainBoundCond[i].boundaryType);
#endif  // End of OPS_2D

#ifdef OPS_3D
            TreatBlockBoundary3D(g_domainBoundCond[i].blockIndex,
                                 g_domainBoundCond[i].componentID,
                                 g_domainBoundCond[i].givenVars,
                                 g_domainBoundCond[i].range,
                                 g_domainBoundCond[i].boundaryType);
#endif  // End of OPS_3D
        }

    }

    else {
        ops_printf("\n No Boundary condition has been defined.");
    }
}

void DefineIntialCond(const int blockIndex, const int componentId,
                      std::vector<Real> initialMacroValues) {
    int numMacroVarsComp;  // Number of macroscopic variables for a given
                           // component.
    numMacroVarsComp = VARIABLECOMPPOS[2 * componentId + 1] -
                       VARIABLECOMPPOS[2 * componentId] + 1;

    const int numInitialMacroComp{(int)initialMacroValues.size()};

    Real* initMacroVal;
    initMacroVal = &initialMacroValues[0];

    if (numInitialMacroComp == numMacroVarsComp) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());

        ops_par_loop(KerSetInitialMacroVarsHilemms,
                     "KerSetInitialMacroVarsHilemms", g_Block[blockIndex],
                     SPACEDIM, iterRng,
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW),
                     ops_arg_gbl(initMacroVal, NUMMACROVAR, "Real", OPS_READ),
                     ops_arg_gbl(&componentId, 1, "int", OPS_READ));
    } else {
        ops_printf(
            "\n For component = %i expected number of macroscopic values "
            "for initial condition = %i however number of values received "
            "= %i",
            componentId, numMacroVarsComp, numInitialMacroComp);
    }
}

void SetupDomainNodeType(int blockIndex, VertexTypes* faceType) {
    int nodeType = (int)Vertex_Fluid;
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));

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
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngImax());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }
    // specify faces for 2D cases, they are actually edges.
    // 0 VG_IP left, 1 VG_IM right
    // 2 VG_JP bottom, 3 VG_JM top
    // 4 VG_KP back, 5 VG_KM front
    nodeType = (int)faceType[0];
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[1];
    iterRange = BlockIterRng(blockIndex, IterRngImax());
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[2];
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[3];
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));

    if (3 == SPACEDIM) {
        nodeType = (int)faceType[4];
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        nodeType = (int)faceType[5];
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }

#if 0

    const int nx = BlockSize(blockIndex)[0];
    const int ny = BlockSize(blockIndex)[1];
    // 2D Domain corner points 4 types
    // 0 IPJP leftBottom, 1 IPJM leftTop
    // 2 IMJP rightBottom, 3 IMJM rightTop
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        nodeType = (int)cornerType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        nodeType = (int)cornerType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        nodeType = (int)cornerType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        nodeType = (int)cornerType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }

    if (3 == SPACEDIM) {
        const int nz = BlockSize(blockIndex)[2];
        // 3D Domain edges 12 types
        // 0 IPJP leftBottom, 1 IPJM leftTop, 2 IMJP rightBottom, 3 IMJM
        // rightTop 4 IPKP leftBack, 5  IPKM leftFront, 6 IMKP rightBack, 7 IMKM
        // rightFront 8 JPKP bottomBack, 9 JPKM bottomFront, 10 JMKP topBack, 11
        // JMKM topFront
        int iminjmin[]{0, 1, 0, 1, 0, nz};
        nodeType = (int)edgeType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        nodeType = (int)edgeType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        nodeType = (int)edgeType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        nodeType = (int)edgeType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminkmin[]{0, 1, 0, ny, 0, 1};
        nodeType = (int)edgeType[4];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        nodeType = (int)edgeType[5];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        nodeType = (int)edgeType[6];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        nodeType = (int)edgeType[7];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jminkmin[]{0, nx, 0, 1, 0, 1};
        nodeType = (int)edgeType[8];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        nodeType = (int)edgeType[9];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        nodeType = (int)edgeType[10];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        nodeType = (int)edgeType[11];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));

        // 3D domain corners 8 types
        // 0 IPJPKP LeftBottomBack 1 IPJPKM LeftBottomFront
        // 2 IPJMKP LeftTopBack 3 IPJMKM LeftTopFront
        // 4 IMJPKP RightBottomBack 5 IMJPKM RightBottomFront
        // 6 IMJMKP RightTopBack 7 IMJMKM RightTopFront
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        nodeType = (int)cornerType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        nodeType = (int)cornerType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        nodeType = (int)cornerType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny-1, ny, nz - 1, nz};
        nodeType = (int)cornerType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        nodeType = (int)cornerType[4];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        nodeType = (int)cornerType[5];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        nodeType = (int)cornerType[6];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny-1, ny, nz - 1, nz};
        nodeType = (int)cornerType[7];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }
#endif  // end of #ifdef 0. Currently disabling this piece of code.
}

void SetupDomainGeometryProperty(int blockIndex) {
    int geometryProperty = (int)VG_Fluid;
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
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
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
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
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
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
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
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
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
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
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
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
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    // specify domain
    geometryProperty = VG_JP;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_JM;
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IP;
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IM;
    iterRange = BlockIterRng(blockIndex, IterRngImax());
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        geometryProperty = VG_KP;
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        geometryProperty = VG_KM;
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

#if 0
    const int nx = BlockSize(blockIndex)[0];
    const int ny = BlockSize(blockIndex)[1];
    // 2D Domain corner points four types
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
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
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int iminkmin[]{0, 1, 0, ny, 0, 1};
        geometryProperty = VG_IPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        geometryProperty = VG_IPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        geometryProperty = VG_IMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        geometryProperty = VG_IMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int jminkmin[]{0, nx, 0, 1, 0, 1};
        geometryProperty = VG_JPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_JPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_JMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_JMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // 3D domain corners 8 types
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        geometryProperty = VG_IPJPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        geometryProperty = VG_IPJPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        geometryProperty = VG_IPJMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny-1, ny, nz - 1, nz};
        geometryProperty = VG_IPJMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        geometryProperty = VG_IMJPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_IMJPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_IMJMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny-1, ny, nz - 1, nz};
        geometryProperty = VG_IMJMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }
#endif  // end of if 0.
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

#ifdef OPS_2d
// mark all solid points inside the circle to be ImmersedSolid
void SolidPointsInsideCircle(int blockIndex, Real diameter,
                             std::vector<Real> circlePos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* circlePosition = &circlePos[0];
    ops_par_loop(
        KerSetEmbededCircle, "KerSetEmbededCircle", g_Block[blockIndex],
        SPACEDIM, bulkRng, ops_arg_gbl(&diameter, 1, "double", OPS_READ),
        ops_arg_gbl(circlePosition, SPACEDIM, "Real", OPS_READ),
        ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL,
                    "double", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_WRITE));
}

void SolidPointsInsideEllipse(int blockIndex, Real semiMajorAxes,
                              Real semiMinorAxes, std::vector<Real> centerPos) {
    int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
    Real* centerPosition = &centerPos[0];
    ops_par_loop(
        KerSetEmbeddedEllipse, "KerSetEmbededCircle", g_Block[blockIndex],
        SPACEDIM, bulkRng, ops_arg_gbl(&semiMajorAxes, 1, "double", OPS_READ),
        ops_arg_gbl(&semiMinorAxes, 1, "double", OPS_READ),
        ops_arg_gbl(centerPosition, SPACEDIM, "Real", OPS_READ),
        ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL,
                    "double", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_WRITE));
}

// Wrapper function for embedded body.
void HandleImmersedSoild() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* bulkRng = BlockIterRng(blockIndex, IterRngBulk());
        // wipe off some solid points that cannot be consideres
        // as a good surface point
        ops_par_loop(KerSweep, "KerSweep", g_Block[blockIndex], SPACEDIM,
                     bulkRng,
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));

        // sync the Geometry property to reflect the modifed solid property
        ops_par_loop(KerSyncGeometryProperty, "KerSyncGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_RW));

        // set the correct  geometry property e.g., corner types
        // i.e., mark out the surface points
        ops_par_loop(KerSetEmbededBodyGeometry, "KerSetEmbededBodyGeometry",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_dat(g_NodeType[blockIndex], 1, ONEPTLATTICESTENCIL,
                                 "int", OPS_RW),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // set the boundary type
        // int nodeType{ surface };
        int nodeType{Vertex_EQMDiffuseRefl};
        ops_par_loop(KerSetEmbededBodyBoundary, "KerSetEmbededBodyBoundary",
                     g_Block[blockIndex], SPACEDIM, bulkRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_RW));
    }
}

// Function to provide details of embedded solid body into the fluid.
void EmbeddedBody(SolidBodyType type, int blockIndex,
                  std::vector<Real> centerPos, std::vector<Real> controlParas) {
    int numCoordCenterPos;
    numCoordCenterPos = centerPos.size();

    if (numCoordCenterPos == SPACEDIM) {
        switch (type) {
            case SolidBody_circle: {
                SolidPointsInsideCircle(blockIndex, controlParas[0], centerPos);
                break;
            }

            case SolidBody_ellipse: {
                Real semiMajorAxes{controlParas[0]};
                Real semiMinorAxes{controlParas[1]};
                SolidPointsInsideEllipse(blockIndex, semiMajorAxes,
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