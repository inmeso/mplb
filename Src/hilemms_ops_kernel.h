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

#ifndef HILEMMS_OPS_KERNEL
#define HILEMMS_OPS_KERNEL
#include "hilemms.h"

void KerSetCoordinates(ACC<Real>& coordinates, const int* idx,
                       const Real* coordX, const Real* coordY) {
#ifdef OPS_2D
    coordinates(0, 0, 0) = coordX[idx[0]];
    coordinates(1, 0, 0) = coordY[idx[1]];
#endif
}

void KerSetCoordinates3D(ACC<Real>& coordinates, const int* idx,
                         const Real* coordX, const Real* coordY,
                         const Real* coordZ) {
#ifdef OPS_3D
    coordinates(0, 0, 0, 0) = coordX[idx[0]];
    coordinates(1, 0, 0, 0) = coordY[idx[1]];
    coordinates(2, 0, 0, 0) = coordZ[idx[2]];
#endif
}

void KerSetEmbeddedBodyBoundary(ACC<int>& nodeType,
                                const ACC<int>& geometryProperty,
                                int* surfaceBoundary) {
#ifdef OPS_2D
    VertexGeometryTypes gp = (VertexGeometryTypes)geometryProperty(0, 0);
    VertexTypes vt = (VertexTypes)nodeType(0, 0);
    if (gp != VG_Fluid && gp != VG_ImmersedSolid && Vertex_Fluid != vt) {
        nodeType(0, 0) = *surfaceBoundary;
    }
#endif
}

void KerSetEmbeddedCircle(ACC<int>& nodeType, ACC<int>& geometryProperty,
                          const Real* coordinates, Real* diameter,
                          Real* centerPos) {
#ifdef OPS_2D
    if ((coordinates(0, 0, 0) - centerPos[0]) *
                (coordinates(0, 0, 0) - centerPos[0]) +
            (coordinates(1, 0, 0) - centerPos[1]) *
                (coordinates(1, 0, 0) - centerPos[1]) <=
        (*diameter) * (*diameter) / 4) {
        nodeType(0, 0) = (int)Vertex_ImmersedSolid;
        geometryProperty(0, 0) = (int)VG_ImmersedSolid;
    }
#endif
}

void KerSetEmbeddedEllipse(ACC<int>& nodeType, ACC<int>& geometryProperty,
                           const ACC<Real>& coordinates, Real* semiMajorAxes,
                           Real* semiMinorAxis, Real* centerPos) {
#ifdef OPS_2D
    if ((coordinates(0, 0, 0) - centerPos[0]) / (*semiMajorAxes) *
                (coordinates(0, 0, 0) - centerPos[0]) / (*semiMajorAxes) +
            (coordinates(1, 0, 0) - centerPos[1]) / (*semiMinorAxis) *
                (coordinates(1, 0, 0) - centerPos[1]) / (*semiMinorAxis) <=
        1.0) {
        nodeType(0, 0) = (int)Vertex_ImmersedSolid;
        geometryProperty(0, 0) = (int)VG_ImmersedSolid;
    }
#endif
}

void KerSweep(ACC<int> nodeType, const ACC<int>& geometryProperty) {
#ifdef OPS_2D
    VertexGeometryTypes vg = (VertexGeometryTypes)geometryProperty(0, 0);
    if (VG_ImmersedSolid == vg) {
        VertexGeometryTypes neiborVertexType[8];
        neiborVertexType[0] = (VertexGeometryTypes)geometryProperty(1, 0);
        neiborVertexType[1] = (VertexGeometryTypes)geometryProperty(-1, 0);
        neiborVertexType[2] = (VertexGeometryTypes)geometryProperty(0, 1);
        neiborVertexType[3] = (VertexGeometryTypes)geometryProperty(0, -1);
        neiborVertexType[4] = (VertexGeometryTypes)geometryProperty(1, 1);
        neiborVertexType[5] = (VertexGeometryTypes)geometryProperty(-1, -1);
        neiborVertexType[6] = (VertexGeometryTypes)geometryProperty(-1, 1);
        neiborVertexType[7] = (VertexGeometryTypes)geometryProperty(1, -1);
        int fluidNeighborNum = 0;
        for (int i = 0; i < 8; i++) {
            if (VG_ImmersedSolid != neiborVertexType[i]) {
                fluidNeighborNum++;
            }
        }

        int solidNeighborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VG_ImmersedSolid == neiborVertexType[i]) {
                solidNeighborNumatCoord++;
            }
        }
        if (fluidNeighborNum > 0 && solidNeighborNumatCoord <= 1) {
            nodeType(0, 0) = Vertex_Fluid;
            ops_printf(
                "A solid point is wiped off due to there are %d fluid points "
                "surrounded and only %d solid points at x and y coordinates\n ",
                fluidNeighborNum, solidNeighborNumatCoord);
        }
    }
#endif
}

void KerSyncGeometryProperty(ACC<int>& geometryProperty,
                             const ACC<int>& nodeType) {
#ifdef OPS_2D
    VertexGeometryTypes gp = (VertexGeometryTypes)geometryProperty(0, 0);
    VertexTypes vt = (VertexTypes)nodeType(0, 0);
    if (Vertex_Fluid == vt && gp != VG_Fluid) {
        geometryProperty(0, 0) = (int)VG_Fluid;
    }
#endif
}

void KerSetEmbeddedBodyGeometry(ACC<int>& geometryProperty,
                                const ACC<int>& nodeType) {
#ifdef OPS_2D
    VertexTypes vt = (VertexTypes)nodeType(0, 0);
    if (Vertex_ImmersedSolid == vt) {
        VertexTypes neiborVertexType[8];
        /*
                        6*****2*****4
                        *     *     *
                        *     *     *
                        1***********0
                        *     *     *
                        *     *     *
                        5*****3*****7
        */
        neiborVertexType[0] = (VertexTypes)nodeType(1, 0);
        neiborVertexType[1] = (VertexTypes)nodeType(-1, 0);
        neiborVertexType[2] = (VertexTypes)nodeType(0, 1);
        neiborVertexType[3] = (VertexTypes)nodeType(0, -1);
        neiborVertexType[4] = (VertexTypes)nodeType(1, 1);
        neiborVertexType[5] = (VertexTypes)nodeType(-1, -1);
        neiborVertexType[6] = (VertexTypes)nodeType(-1, 1);
        neiborVertexType[7] = (VertexTypes)nodeType(1, -1);
        int fluidNeiborNum{0};
        for (int i = 0; i < 8; i++) {
            if (Vertex_ImmersedSolid != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (Vertex_ImmersedSolid == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        if (fluidNeiborNum > 0) {
            // outer corner
            if (2 == solidNeiborNumatCoord) {
                if ((Vertex_ImmersedSolid == neiborVertexType[0] &&
                     Vertex_ImmersedSolid == neiborVertexType[1]) ||
                    (Vertex_ImmersedSolid == neiborVertexType[0] &&
                     Vertex_ImmersedSolid == neiborVertexType[1])) {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points,i.e., "
                               "the solid body may be too thin");
                }
                if (Vertex_ImmersedSolid == neiborVertexType[2] &&
                    Vertex_ImmersedSolid == neiborVertexType[1]) {
                    if (Vertex_ImmersedSolid == neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IPJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if (Vertex_ImmersedSolid == neiborVertexType[3] &&
                    Vertex_ImmersedSolid == neiborVertexType[1]) {
                    if (Vertex_ImmersedSolid == neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if (Vertex_ImmersedSolid == neiborVertexType[3] &&
                    Vertex_ImmersedSolid == neiborVertexType[0]) {
                    if (Vertex_ImmersedSolid == neiborVertexType[7]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
                if (Vertex_ImmersedSolid == neiborVertexType[2] &&
                    Vertex_ImmersedSolid == neiborVertexType[0]) {
                    if (Vertex_ImmersedSolid == neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
            }
            // Planlar corner
            if (3 == solidNeiborNumatCoord) {
                if (Vertex_ImmersedSolid != neiborVertexType[0]) {
                    geometryProperty(0, 0) = (int)VG_IP;
                }
                if (Vertex_ImmersedSolid != neiborVertexType[1]) {
                    geometryProperty(0, 0) = (int)VG_IM;
                }
                if (Vertex_ImmersedSolid != neiborVertexType[2]) {
                    geometryProperty(0, 0) = (int)VG_JP;
                }
                if (Vertex_ImmersedSolid != neiborVertexType[3]) {
                    geometryProperty(0, 0) = (int)VG_JM;
                }
            }
            // Inner corner
            if (4 == solidNeiborNumatCoord) {
                if (1 == fluidNeiborNum) {
                    if (Vertex_ImmersedSolid != neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_I;
                    }
                    if (Vertex_ImmersedSolid != neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_I;
                    }
                    if (Vertex_ImmersedSolid != neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_I;
                    }
                    if (Vertex_ImmersedSolid != neiborVertexType[7]) {
                        geometryProperty(0, 0) = (int)VG_IPJM_I;
                    }
                } else {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points");
                }
            }
        }
    }
#endif
}

#endif  // HILEMMS_OPS_KERNEL