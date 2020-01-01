
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

/*! @brief   Head files for importing geometry from HDF5 file
  * @author  Jianping Meng
  * @details Declaring kernel functions related to create computing domain
  */

#ifndef SETUP_COMPUT_DOMAIN_KERNEL_H
#define SETUP_COMPUT_DOMAIN_KERNEL_H
#include "setup_comput_domain.h"
#ifdef OPS_2D
void KerSetCoordinates(const Real* coordX, const Real* coordY, const int* idx,
                       Real* coordinates) {
    coordinates[OPS_ACC_MD3(0, 0, 0)] = coordX[idx[0]];
    coordinates[OPS_ACC_MD3(1, 0, 0)] = coordY[idx[1]];
}
void KerSetEmbeddedCircle(Real* diameter, Real* centerPos,
                         const Real* coordinates, int* nodeType,
                         int* geometryProperty) {
    if ((coordinates[OPS_ACC_MD3(0, 0, 0)] - centerPos[0]) *
                (coordinates[OPS_ACC_MD3(0, 0, 0)] - centerPos[0]) +
            (coordinates[OPS_ACC_MD3(1, 0, 0)] - centerPos[1]) *
                (coordinates[OPS_ACC_MD3(1, 0, 0)] - centerPos[1]) <=
        (*diameter) * (*diameter) / 4) {
        nodeType[OPS_ACC4(0, 0)] = (int)VertexType::ImmersedSolid;
        geometryProperty[OPS_ACC5(0, 0)] = (int)VG_ImmersedSolid;
    }
}
void KerSweep(const int* geometryProperty, int* nodeType) {
    VertexGeometryTypes vg =
        (VertexGeometryTypes)geometryProperty[OPS_ACC0(0, 0)];
    if (VG_ImmersedSolid == vg) {
        VertexGeometryTypes neiborVertexType[8];
        neiborVertexType[0] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(1, 0)];
        neiborVertexType[1] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(-1, 0)];
        neiborVertexType[2] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(0, 1)];
        neiborVertexType[3] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(0, -1)];
        neiborVertexType[4] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(1, 1)];
        neiborVertexType[5] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(-1, -1)];
        neiborVertexType[6] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(-1, 1)];
        neiborVertexType[7] =
            (VertexGeometryTypes)geometryProperty[OPS_ACC0(1, -1)];
        int fluidNeiborNum = 0;
        for (int i = 0; i < 8; i++) {
            if (VG_ImmersedSolid != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        // A point
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VG_ImmersedSolid == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        if (fluidNeiborNum > 0 && solidNeiborNumatCoord <= 1) {
            nodeType[OPS_ACC1(0, 0)] = VertexType::Fluid;
            ops_printf(
                "A solid point is wiped off due to there are %d fluid points "
                "surrounded and only %d solid points at x and y coordinates\n ",
                fluidNeiborNum, solidNeiborNumatCoord);
        }
    }
}

void KerSyncGeometryProperty(const int* nodeType, int* geometryProperty) {
    VertexGeometryTypes gp =
        (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
    VertexType vt = (VertexType)nodeType[OPS_ACC0(0, 0)];
    if (VertexType::Fluid == vt && gp != VG_Fluid) {
        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_Fluid;
    }
}

void KerSetEmbeddedBodyGeometry(const int* nodeType, int* geometryProperty) {
    VertexType vt = (VertexType)nodeType[OPS_ACC1(0, 0)];
    if (VertexType::ImmersedSolid == vt) {
        VertexType neiborVertexType[8];
        /*
                        6*****2*****4
                        *     *     *
                        *     *     *
                        1***********0
                        *     *     *
                        *     *     *
                        5*****3*****7
        */
        neiborVertexType[0] = (VertexType)nodeType[OPS_ACC0(1, 0)];
        neiborVertexType[1] = (VertexType)nodeType[OPS_ACC0(-1, 0)];
        neiborVertexType[2] = (VertexType)nodeType[OPS_ACC0(0, 1)];
        neiborVertexType[3] = (VertexType)nodeType[OPS_ACC0(0, -1)];
        neiborVertexType[4] = (VertexType)nodeType[OPS_ACC0(1, 1)];
        neiborVertexType[5] = (VertexType)nodeType[OPS_ACC0(-1, -1)];
        neiborVertexType[6] = (VertexType)nodeType[OPS_ACC0(-1, 1)];
        neiborVertexType[7] = (VertexType)nodeType[OPS_ACC0(1, -1)];
        int fluidNeiborNum{0};
        for (int i = 0; i < 8; i++) {
            if (VertexType::ImmersedSolid != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        if (fluidNeiborNum > 0) {
            // outer corner
            if (2 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1]) ||
                    (VertexType::ImmersedSolid == neiborVertexType[2] &&
                     VertexType::ImmersedSolid == neiborVertexType[3])) {
                    ops_printf("%s\n", "This may be a dangling point.");
                }
                if (VertexType::ImmersedSolid == neiborVertexType[2] &&
                    VertexType::ImmersedSolid == neiborVertexType[1]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IPJM_O;
                    } else {
                        ops_printf("%s\n", "The solid body may be too thin.");
                    }
                }

                if (VertexType::ImmersedSolid == neiborVertexType[3] &&
                    VertexType::ImmersedSolid == neiborVertexType[1]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[5]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IPJP_O;
                    } else {
                       ops_printf("%s\n", "The solid body may be too thin.");
                    }
                }

                if (VertexType::ImmersedSolid == neiborVertexType[3] &&
                    VertexType::ImmersedSolid == neiborVertexType[0]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IMJP_O;
                    } else {
                       ops_printf("%s\n", "The solid body may be too thin.");
                    }
                }
                if (VertexType::ImmersedSolid == neiborVertexType[2] &&
                    VertexType::ImmersedSolid == neiborVertexType[0]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[4]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IMJM_O;
                    } else {
                       ops_printf("%s\n", "The solid body may be too thin.");
                    }
                }
            }
            // Planlar corner
            if (3 == solidNeiborNumatCoord) {
                if (VertexType::ImmersedSolid != neiborVertexType[0]) {
                    geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]) {
                    geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IM;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]) {
                    geometryProperty[OPS_ACC1(0, 0)] = (int)VG_JP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]) {
                    geometryProperty[OPS_ACC1(0, 0)] = (int)VG_JM;
                }
            }
            // Inner corner
            if (4 == solidNeiborNumatCoord) {
                if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[4]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[5]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[6]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]) {
                        geometryProperty[OPS_ACC1(0, 0)] = (int)VG_IPJM_I;
                    }
                } else {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points");
                }
            }
        }
    }
}
#endif  // OPS_2D
#ifdef OPS_3D
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx,
                         Real* coordinates) {
    coordinates[OPS_ACC_MD4(0, 0, 0, 0)] = coordX[idx[0]];
    coordinates[OPS_ACC_MD4(1, 0, 0, 0)] = coordY[idx[1]];
    coordinates[OPS_ACC_MD4(2, 0, 0, 0)] = coordZ[idx[2]];
}
#endif  // OPS_3D

void KerSetInitialMacroVars(const Real* coordinates, const int* idx,
                            Real* macroVars) {
#ifdef OPS_2D
    Real dist = (coordinates[OPS_ACC_MD0(0, 0, 0)]-1) * (coordinates[OPS_ACC_MD0(0, 0, 0)]-1) +(coordinates[OPS_ACC_MD0(1, 0, 0)]-1) * (coordinates[OPS_ACC_MD0(1, 0, 0)]-1);
    if (dist <=(0.35*0.35)){
        macroVars[OPS_ACC_MD2(0, 0, 0)] = 1;
        macroVars[OPS_ACC_MD2(1, 0, 0)] = 0;
        macroVars[OPS_ACC_MD2(2, 0, 0)] = 0;
    } else {
        macroVars[OPS_ACC_MD2(0, 0, 0)] = 0.1;
        macroVars[OPS_ACC_MD2(1, 0, 0)] = 0;
        macroVars[OPS_ACC_MD2(2, 0, 0)] = 0;
    }
#endif
#ifdef OPS_3D
    macroVars[OPS_ACC_MD2(0, 0, 0, 0)] = 1;
    macroVars[OPS_ACC_MD2(1, 0, 0, 0)] = 0;//sin(coordinates[OPS_ACC_MD0(0, 0, 0, 0)]); //u
    macroVars[OPS_ACC_MD2(2, 0, 0, 0)] = 0;         // v
    macroVars[OPS_ACC_MD2(3, 0, 0, 0)] = 0;         // w
#endif
}

void KerSetEmbeddedBodyBoundary(int* surfaceBoundary,
                               const int* geometryProperty, int* nodeType) {
#ifdef OPS_2D
    VertexGeometryTypes gp =
        (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
    VertexType vt = (VertexType)nodeType[OPS_ACC2(0, 0)];
    if (gp != VG_Fluid && gp != VG_ImmersedSolid && VertexType::Fluid != vt) {
        nodeType[OPS_ACC2(0, 0)] = *surfaceBoundary;
    }
#endif
#ifdef OPS_3D
    VertexGeometryType gp =
        (VertexGeometryType)geometryProperty[OPS_ACC1(0, 0, 0)];
    VertexType vt = (VertexType)nodeType[OPS_ACC2(0, 0, 0)];
    if (gp != VG_Fluid && gp != VG_ImmersedSolid && VertexType::Fluid != vt) {
        nodeType[OPS_ACC2(0, 0, 0)] = *surfaceBoundary;
    }
#endif
}

#endif /* SETUP_COMPUT_DOMAIN_KERNEL_H */
