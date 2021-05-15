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
//Note: this function can be dropped out for the new vertexType
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