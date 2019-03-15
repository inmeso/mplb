// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief   Kernel functions for boundary conditions
 * @author  Jianping Meng
 * @details Defining kernel functions for various boundary conditions.
 */
#ifndef BOUNDARY_KERNEL_H
#define BOUNDARY_KERNEL_H
#include <iostream>
#include "boundary.h"
// As we are using update-halo method for the discretisation,
// we need to deal with halo points when treating boundary
/*
 * TODO: for a corner point, there may be a issue for boundary conditions
 * including the free flux boundary since they also rely on the surface
 * direction. For example, the free flux boundary at a VG_IPJP_I point may not
 * know which direction to pose the zero gradient, although this may have
 * negligible effect.
 */

#ifdef OPS_2D
void KerCutCellZeroFlux(const int *nodeType, const int *geometryProperty,
                        Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    if (vt == Vertex_FreeFlux) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
        switch (vg) {
            case VG_IP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 1, 0)];
                }
                break;
            case VG_IM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, -1, 0)];
                }
                break;
            case VG_JP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, 1)];
                }
                break;
            case VG_JM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, -1)];
                }
                break;
            case VG_IPJP_I:                          // inner corner point
                if (vt == nodeType[OPS_ACC0(0, 1)])  // VG_IP
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(1, 0)])  // VG_JP
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1)];
                    }
                }
                break;
            case VG_IPJM_I:                           // inner corner point
                if (vt == nodeType[OPS_ACC0(0, -1)])  // VG_IP
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(1, 0)])  // VG_JM
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1)];
                    }
                }
                break;
            case VG_IMJP_I:                          // inner corner point
                if (vt == nodeType[OPS_ACC0(0, 1)])  // VG_IM
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(-1, 0)])  // VG_JP
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1)];
                    }
                }
                break;
            case VG_IMJM_I:                           // inner corner point
                if (vt == nodeType[OPS_ACC0(0, -1)])  // VG_IM
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(-1, 0)])  // VG_JM
                {
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1)];
                    }
                }
                break;
            default:
                break;
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a free flux boundary "
                   "point: KerCutCellZeroFlux");
#endif
    }
}

void KerCutCellEmbededBoundary(const int *nodeType, const int *geometryProperty,
                               Real *f) {
    /*!
     For the bounce back scheme,We consider zero velocity boundary first.
     To make sure the velocity at boundary is zero, the implementation
     is lattice specific.
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    VertexGeometryTypes vg =
        (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
    if (vt >= Vertex_Boundary) {
        switch (vt) {
            case Vertex_EQMDiffuseRefl: {
                Real u{0};
                Real v{0};
                const Real sqrt3 = sqrt(3);
                switch (vg) {
                    case VG_IP: {
                        const Real f3 = f[OPS_ACC_MD2(3, 0, 0)];
                        const Real f7 = f[OPS_ACC_MD2(7, 0, 0)];
                        const Real f6 = f[OPS_ACC_MD2(6, 0, 0)];
                        const Real rhow =
                            6 * (f3 + f6 + f7) / (u * u - sqrt3 * u + 1);
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            f7 + rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            f3 + 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            f6 + rhow * (u - v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            -(u * u - 2 * (1 + sqrt3 * v + v * v) * rhow) / 18;
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            -((-2 + u * u + 2 * sqrt3 * v - 2 * v * v) * rhow) /
                            18;
                    } break;
                    case VG_IM: {
                        const Real f5 = f[OPS_ACC_MD2(5, 0, 0)];
                        const Real f1 = f[OPS_ACC_MD2(1, 0, 0)];
                        const Real f8 = f[OPS_ACC_MD2(8, 0, 0)];
                        const Real rhow =
                            6 * (f1 + f5 + f8) / (u * u + sqrt3 * u + 1);
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            f5 - rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            f1 - 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            f8 + rhow * (v - u) / (6 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            -(u * u - 2 * (1 + sqrt3 * v + v * v) * rhow) / 18;
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            -((-2 + u * u + 2 * sqrt3 * v - 2 * v * v) * rhow) /
                            18;
                    } break;
                    case VG_JP: {
                        const Real f4 = f[OPS_ACC_MD2(4, 0, 0)];
                        const Real f8 = f[OPS_ACC_MD2(8, 0, 0)];
                        const Real f7 = f[OPS_ACC_MD2(7, 0, 0)];
                        const Real rhow =
                            6 * (f4 + f8 + f7) / (v * v - sqrt3 * v + 1);
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            f4 + 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            f8 + rhow * (v - u) / (6 * sqrt3);
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            f7 + rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            ((2 + 2 * sqrt3 * u + 2 * u * u - v * v) * rhow) /
                            18;
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            -((-2 + 2 * sqrt3 * u - 2 * u * u + v * v) * rhow) /
                            18;
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                    } break;
                    case VG_JM: {
                        const Real f2 = f[OPS_ACC_MD2(2, 0, 0)];
                        const Real f5 = f[OPS_ACC_MD2(5, 0, 0)];
                        const Real f6 = f[OPS_ACC_MD2(6, 0, 0)];
                        const Real rhow =
                            6 * (f2 + f5 + f6) / (v * v + sqrt3 * v + 1);
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            f2 - 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            f6 + rhow * (u - v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            f5 - rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            ((2 + 2 * sqrt3 * u + 2 * u * u - v * v) * rhow) /
                            18;
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            -((-2 + 2 * sqrt3 * u - 2 * u * u + v * v) * rhow) /
                            18;
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                    } break;
                    case VG_IPJP_I: {  // inner corner point
                        const Real f3 = f[OPS_ACC_MD2(3, 0, 0)];
                        const Real f7 = f[OPS_ACC_MD2(7, 0, 0)];
                        const Real f4 = f[OPS_ACC_MD2(4, 0, 0)];
                        const Real rhow =
                            (-36 * (f3 + f4 + f7)) /
                            (-9 + 5 * sqrt3 * u - 3 * u * u + 5 * sqrt3 * v -
                             3 * u * v - 3 * v * v);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            f3 + 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            f7 + rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            f4 + 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v -
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            (1 + u * u + u * (sqrt3 - 3 * v) - sqrt3 * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IPJM_I: {  // inner corner point
                        const Real f2 = f[OPS_ACC_MD2(2, 0, 0)];
                        const Real f3 = f[OPS_ACC_MD2(3, 0, 0)];
                        const Real f6 = f[OPS_ACC_MD2(6, 0, 0)];
                        const Real rhow =
                            (-36 * (f2 + f3 + f6)) /
                            (-9 + 5 * sqrt3 * u - 3 * u * u - 5 * sqrt3 * v +
                             3 * u * v - 3 * v * v);
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            f6 + rhow * (u - v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            f3 + 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            f2 - 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v +
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            (1 - sqrt3 * u + u * u - sqrt3 * v + 3 * u * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IMJP_I: {  // inner corner point
                        const Real f1 = f[OPS_ACC_MD2(1, 0, 0)];
                        const Real f4 = f[OPS_ACC_MD2(4, 0, 0)];
                        const Real f8 = f[OPS_ACC_MD2(8, 0, 0)];
                        const Real rhow =
                            (36 * (f1 + f4 + f8)) /
                            (9 + 5 * sqrt3 * u + 3 * u * u - 5 * sqrt3 * v -
                             3 * u * v + 3 * v * v);
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            f8 + rhow * (v - u) / (6 * sqrt3);
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            f1 - 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            f4 + 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v +
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            (1 - sqrt3 * u + u * u - sqrt3 * v + 3 * u * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IMJM_I: {  // inner corner point
                        const Real f1 = f[OPS_ACC_MD2(1, 0, 0)];
                        const Real f5 = f[OPS_ACC_MD2(5, 0, 0)];
                        const Real f2 = f[OPS_ACC_MD2(2, 0, 0)];
                        const Real rhow =
                            (36 * (f1 + f2 + f5)) /
                            (9 + 5 * sqrt3 * u + 3 * u * u + 5 * sqrt3 * v +
                             3 * u * v + 3 * v * v);
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            f1 - 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            f2 - 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            f5 - rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v -
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            (1 + u * u + u * (sqrt3 - 3 * v) - sqrt3 * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IPJP_O: {  // outter corner point
                        const Real f3 = f[OPS_ACC_MD2(3, 0, 0)];
                        const Real f7 = f[OPS_ACC_MD2(7, 0, 0)];
                        const Real f4 = f[OPS_ACC_MD2(4, 0, 0)];
                        const Real rhow =
                            (-36 * (f3 + f4 + f7)) /
                            (-9 + 5 * sqrt3 * u - 3 * u * u + 5 * sqrt3 * v -
                             3 * u * v - 3 * v * v);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            f3 + 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            f7 + rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            f4 + 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v -
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            (1 + u * u + u * (sqrt3 - 3 * v) - sqrt3 * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IPJM_O: {  // outter corner point
                        const Real f2 = f[OPS_ACC_MD2(2, 0, 0)];
                        const Real f3 = f[OPS_ACC_MD2(3, 0, 0)];
                        const Real f6 = f[OPS_ACC_MD2(6, 0, 0)];
                        const Real rhow =
                            (-36 * (f2 + f3 + f6)) /
                            (-9 + 5 * sqrt3 * u - 3 * u * u - 5 * sqrt3 * v +
                             3 * u * v - 3 * v * v);
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            f6 + rhow * (u - v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(1, 0, 0)] =
                            f3 + 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            f2 - 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v +
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            (1 - sqrt3 * u + u * u - sqrt3 * v + 3 * u * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IMJP_O: {  // outter corner point
                        const Real f1 = f[OPS_ACC_MD2(1, 0, 0)];
                        const Real f4 = f[OPS_ACC_MD2(4, 0, 0)];
                        const Real f8 = f[OPS_ACC_MD2(8, 0, 0)];
                        const Real rhow =
                            (36 * (f1 + f4 + f8)) /
                            (9 + 5 * sqrt3 * u + 3 * u * u - 5 * sqrt3 * v -
                             3 * u * v + 3 * v * v);
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            f8 + rhow * (v - u) / (6 * sqrt3);
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            f1 - 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(2, 0, 0)] =
                            f4 + 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(5, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v +
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            (1 - sqrt3 * u + u * u - sqrt3 * v + 3 * u * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    case VG_IMJM_O: {  // outter corner point
                        const Real f1 = f[OPS_ACC_MD2(1, 0, 0)];
                        const Real f5 = f[OPS_ACC_MD2(5, 0, 0)];
                        const Real f2 = f[OPS_ACC_MD2(2, 0, 0)];
                        const Real rhow =
                            (36 * (f1 + f2 + f5)) /
                            (9 + 5 * sqrt3 * u + 3 * u * u + 5 * sqrt3 * v +
                             3 * u * v + 3 * v * v);
                        f[OPS_ACC_MD2(3, 0, 0)] =
                            f1 - 2 * rhow * u / (3 * sqrt3);
                        f[OPS_ACC_MD2(4, 0, 0)] =
                            f2 - 2 * rhow * v / (3 * sqrt3);
                        f[OPS_ACC_MD2(7, 0, 0)] =
                            f5 - rhow * (u + v) / (6 * sqrt3);
                        f[OPS_ACC_MD2(0, 0, 0)] =
                            2 * rhow * (2 - u * u - v * v) / 9;
                        f[OPS_ACC_MD2(6, 0, 0)] =
                            (1 + u * u + sqrt3 * v + v * v -
                             u * (sqrt3 + 3 * v)) *
                            rhow / 36;
                        f[OPS_ACC_MD2(8, 0, 0)] =
                            (1 + u * u + u * (sqrt3 - 3 * v) - sqrt3 * v +
                             v * v) *
                            rhow / 36;
                    } break;
                    default:
                        break;
                }  // vg
            }      // case Vertex_EQMDiffRefl
            break;
            case Vertex_KineticDiffuseWall: {
                Real u = 0;
                Real v = 0;  // means non-moving wall
                Real wallNormalVector[]{0, 0};
                const Real sqrt2Inverse = 1 / sqrt(2);
                switch (vg) {
                    case VG_IP: {
                        wallNormalVector[0] = 1;
                        wallNormalVector[1] = 0;
                    } break;
                    case VG_IM: {
                        wallNormalVector[0] = -1;
                        wallNormalVector[1] = 0;
                    } break;
                    case VG_JP: {
                        wallNormalVector[0] = 0;
                        wallNormalVector[1] = 1;
                    } break;
                    case VG_JM: {
                        wallNormalVector[0] = 0;
                        wallNormalVector[1] = -1;
                    } break;
                    case VG_IPJP_I:  // inner corner point
                    {
                        wallNormalVector[0] = sqrt2Inverse;
                        wallNormalVector[1] = sqrt2Inverse;
                    } break;
                    case VG_IPJM_I:  // inner corner point
                    {
                        wallNormalVector[0] = sqrt2Inverse;
                        wallNormalVector[1] = -sqrt2Inverse;
                    } break;
                    case VG_IMJP_I:  // inner corner point
                    {
                        wallNormalVector[0] = -sqrt2Inverse;
                        wallNormalVector[1] = sqrt2Inverse;
                    } break;
                    case VG_IMJM_I:  // inner corner point
                    {
                        wallNormalVector[0] = -sqrt2Inverse;
                        wallNormalVector[1] = -sqrt2Inverse;
                    } break;
                    case VG_IPJP_O: {
                        wallNormalVector[0] = sqrt2Inverse;
                        wallNormalVector[1] = sqrt2Inverse;
                    } break;
                    case VG_IPJM_O: {
                        wallNormalVector[0] = sqrt2Inverse;
                        wallNormalVector[1] = -sqrt2Inverse;
                    } break;
                    case VG_IMJP_O: {
                        wallNormalVector[0] = -sqrt2Inverse;
                        wallNormalVector[1] = sqrt2Inverse;
                    } break;
                    case VG_IMJM_O: {
                        wallNormalVector[0] = -sqrt2Inverse;
                        wallNormalVector[1] = -sqrt2Inverse;
                    } break;

                    default:
                        break;
                }
                Real outFlux = 0;  // flow into wall
                Real inFlux = 0;   // flow into fluid bulk
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    const Real cx = XI[xiIndex * LATTDIM];
                    const Real cy = XI[xiIndex * LATTDIM + 1];
                    Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
                                      (CS * cy - v) * wallNormalVector[1];
                    if (cDotNormal < 0) {
                        outFlux +=
                            (-cDotNormal * f[OPS_ACC_MD2(xiIndex, 0, 0)]);
                    }
                    if (cDotNormal > 0) {
                        Real cu = (CS * cx * u + CS * cy * v);
                        inFlux +=
                            (cDotNormal * WEIGHTS[xiIndex] *
                             (1 + cu +
                              0.5 * (cu * cu - (u * u + v * v))));  // i.e., the
                        // equilibrium
                    }
                }
                Real rho = outFlux / inFlux;
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    const Real cx = XI[xiIndex * LATTDIM];
                    const Real cy = (int)XI[xiIndex * LATTDIM + 1];
                    Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
                                      (CS * cy - v) * wallNormalVector[1];
                    if (cDotNormal >= 0) {
                        Real cu = (CS * cx * u + CS * cy * v);
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            (rho * WEIGHTS[xiIndex] *
                             (1 + cu + 0.5 * (cu * cu - (u * u + v * v))));
                    }
                }
            }  // case Vertex_KineticDiffuseWall
            break;
            default:
#ifdef debug
                ops_printf("%s\n",
                           "Warning: KerCutCellImmersedBoundary: there seems a "
                           "boundary condition that has note been implemented");
#endif
                break;
        }
    }
}

void KerCutCellExtrapolPressure1ST(const Real *givenBoundaryVars,
                                   const int *nodeType,
                                   const int *geometryProperty, Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_ExtrapolPressure1ST) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        Real rhoGiven = givenBoundaryVars[0];
        Real rho = 0;
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            const int cx = (int)XI[xiIndex * LATTDIM];
            const int cy = (int)XI[xiIndex * LATTDIM + 1];
            switch (vg) {
                case VG_IP: {
                    if (cx > 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD3(xiIndex, 1, 0)];
                    }
                } break;
                case VG_IM: {
                    if (cx < 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD3(xiIndex, -1, 0)];
                    }
                } break;
                case VG_JP: {
                    if (cy > 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD3(xiIndex, 0, 1)];
                    }
                } break;
                case VG_JM: {
                    if (cy < 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD3(xiIndex, 0, -1)];
                    }

                } break;
                case VG_IPJP_I: {                        // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, 1)])  // VG_IP
                    {
                        if ((cx > 0 && cy >= 0) || (cx >= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 1, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(1, 0)])  // VG_JP
                    {
                        if ((cx > 0 && cy >= 0) || (cx >= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 0, 1)];
                        }
                    }
                } break;
                case VG_IPJM_I: {                         // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, -1)])  // VG_IP
                    {
                        if ((cx > 0 && cy <= 0) || (cx >= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 1, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(1, 0)])  // VG_JM
                    {
                        if ((cx > 0 && cy <= 0) || (cx >= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 0, -1)];
                        }
                    }
                } break;

                case VG_IMJP_I: {                        // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, 1)])  // VG_IM
                    {
                        if ((cx < 0 && cy >= 0) || (cx <= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, -1, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(-1, 0)])  // VG_JP
                    {
                        if ((cx < 0 && cy >= 0) || (cx <= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 0, 1)];
                        }
                    }
                } break;
                case VG_IMJM_I: {                         // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, -1)])  // VG_IM
                    {
                        if ((cx < 0 && cy <= 0) || (cx <= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, -1, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(-1, 0)])  // VG_JM
                    {
                        if ((cx < 0 && cy <= 0) || (cx <= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD3(xiIndex, 0, -1)];
                        }
                    }
                } break;
                default:
                    break;
            }
            rho += f[OPS_ACC_MD3(xiIndex, 0, 0)];
        }
        Real ratio = rhoGiven / rho;
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            f[OPS_ACC_MD3(xiIndex, 0, 0)] *= ratio;
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a extrapol pressure "
                   "boundary point: KerCutCellExtraolPressure1ST");
#endif
    }
}

void KerCutCellExtrapolPressure2ND(const Real *givenBoundaryVars,
                                   const int *nodeType,
                                   const int *geometryProperty, Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_ExtrapolPressure2ND) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        Real rhoGiven = givenBoundaryVars[0];
        Real rho = 0;
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            const int cx = (int)XI[xiIndex * LATTDIM];
            const int cy = (int)XI[xiIndex * LATTDIM + 1];
            switch (vg) {
                case VG_IP: {
                    if (cx > 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            2 * f[OPS_ACC_MD3(xiIndex, 1, 0)] -
                            f[OPS_ACC_MD3(xiIndex, 2, 0)];
                    }
                } break;
                case VG_IM: {
                    if (cx < 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            2 * f[OPS_ACC_MD3(xiIndex, -1, 0)] -
                            f[OPS_ACC_MD3(xiIndex, -2, 0)];
                    }
                } break;
                case VG_JP: {
                    if (cy > 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            2 * f[OPS_ACC_MD3(xiIndex, 0, 1)] -
                            f[OPS_ACC_MD3(xiIndex, 0, 2)];
                    }
                } break;
                case VG_JM: {
                    if (cy < 0) {
                        f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                            2 * f[OPS_ACC_MD3(xiIndex, 0, -1)] -
                            f[OPS_ACC_MD3(xiIndex, 0, -2)];
                    }
                } break;
                case VG_IPJP_I: {                        // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, 1)])  // VG_IP
                    {
                        if ((cx > 0 && cy >= 0) || (cx >= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 1, 0)] -
                                f[OPS_ACC_MD3(xiIndex, 2, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(1, 0)])  // VG_JP
                    {
                        if ((cx > 0 && cy >= 0) || (cx >= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 0, 1)] -
                                f[OPS_ACC_MD3(xiIndex, 0, 2)];
                        }
                    }
                } break;
                case VG_IPJM_I: {                         // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, -1)])  // VG_IP
                    {
                        if ((cx > 0 && cy <= 0) || (cx >= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 1, 0)] -
                                f[OPS_ACC_MD3(xiIndex, 2, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(1, 0)])  // VG_JM
                    {
                        if ((cx > 0 && cy <= 0) || (cx >= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 0, -1)] -
                                f[OPS_ACC_MD3(xiIndex, 0, -2)];
                        }
                    }
                } break;
                case VG_IMJP_I: {                        // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, 1)])  // VG_IM
                    {
                        if ((cx < 0 && cy >= 0) || (cx <= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, -1, 0)] -
                                f[OPS_ACC_MD3(xiIndex, -2, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(-1, 0)])  // VG_JP
                    {
                        if ((cx < 0 && cy >= 0) || (cx <= 0 && cy > 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 0, 1)] -
                                f[OPS_ACC_MD3(xiIndex, 0, 2)];
                        }
                    }
                } break;
                case VG_IMJM_I: {                         // inner corner point
                    if (vt == nodeType[OPS_ACC1(0, -1)])  // VG_IM
                    {
                        if ((cx < 0 && cy <= 0) || (cx <= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, -1, 0)] -
                                f[OPS_ACC_MD3(xiIndex, -2, 0)];
                        }
                    }
                    if (vt == nodeType[OPS_ACC1(-1, 0)])  // VG_JM
                    {
                        if ((cx < 0 && cy <= 0) || (cx <= 0 && cy < 0)) {
                            f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                                2 * f[OPS_ACC_MD3(xiIndex, 0, -1)] -
                                f[OPS_ACC_MD3(xiIndex, 0, -2)];
                        }
                    }
                } break;
                default:
                    break;
            }
            rho += f[OPS_ACC_MD3(xiIndex, 0, 0)];
        }
        Real ratio = rhoGiven / rho;
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            f[OPS_ACC_MD3(xiIndex, 0, 0)] *= ratio;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a extrapol pressure "
                   "boundary point: KerCutCellExtraolPressure2ND");
#endif
    }
}

/*
void KerCutCellKinetic(const Real *givenMacroVars, const int *nodeType,
                       const int *geometryProperty, Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_KineticDiffuseWall) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        Real u = givenMacroVars[1];
        Real v = givenMacroVars[2];
        Real T = givenMacroVars[3];
        Real primaryVector[]{0, 0};
        Real secondVector[]{0, 0};
        int boundaryType{0};  // 0 normal boundary 1 inner corner 2 outter
                              // corner
        switch (vg) {
            case VG_IP: {
                boundaryType = 0;
                primaryVector[0] = 1;
                primaryVector[1] = 0;
            } break;
            case VG_IM: {
                boundaryType = 0;
                primaryVector[0] = -1;
                primaryVector[1] = 0;
            } break;
            case VG_JP: {
                boundaryType = 0;
                primaryVector[0] = 0;
                primaryVector[1] = 1;
            } break;
            case VG_JM: {
                boundaryType = 0;
                primaryVector[0] = 0;
                primaryVector[1] = -1;
            } break;
            case VG_IPJP_I:  // inner corner point
            {
                boundaryType = 1;
                primaryVector[0] = 1;
                primaryVector[1] = 0;
                secondVector[0] = 0;
                secondVector[1] = 1;
            } break;
            case VG_IPJM_I:  // inner corner point
            {
                boundaryType = 1;
                primaryVector[0] = 1;
                primaryVector[1] = 0;
                secondVector[0] = 0;
                secondVector[1] = -1;
            } break;
            case VG_IMJP_I:  // inner corner point
            {
                boundaryType = 1;
                primaryVector[0] = -1;
                primaryVector[1] = 0;
                secondVector[0] = 0;
                secondVector[1] = 1;
            } break;
            case VG_IMJM_I:  // inner corner point
            {
                boundaryType = 1;
                primaryVector[0] = -1;
                primaryVector[1] = 0;
                secondVector[0] = 0;
                secondVector[1] = -1;
            } break;
            default:
                break;
        }
        Real outFlux = 0;  // flow into wall
        Real inFlux = 0;   // flow into fluid bulk
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            const Real relVeloX = CS * XI[xiIndex * LATTDIM] - u;
            const Real relVeloY = CS * XI[xiIndex * LATTDIM + 1] - v;
            const Real speed = relVeloX * relVeloX + relVeloY * relVeloY;
            Real cDotPrimary =
                relVeloX * primaryVector[0] + relVeloY * primaryVector[1];
            bool isInflux = cDotPrimary > 0;
            bool isOutFlux = cDotPrimary < 0;
            //ops_printf("1 isIn=%i isOut= %i speed=
%f\n",isInflux,isOutFlux,speed); if (1 == boundaryType) { Real cDotSecond =
                    relVeloX * secondVector[0] + relVeloY * secondVector[1];
                isInflux = isInflux && (cDotSecond > 0);
                isOutFlux = isOutFlux && (cDotSecond < 0);
                //ops_printf("2 isIn=%i isOut= %i\n",isInflux,isOutFlux);
            }
            if (2 == boundaryType) {
                Real cDotSecond =
                    relVeloX * secondVector[0] + relVeloY * secondVector[1];
                isInflux = isInflux && (cDotSecond > 0);
                isOutFlux = isOutFlux && (cDotSecond < 0);
            }

            if (isOutFlux) {
                outFlux += (speed * f[OPS_ACC_MD3(xiIndex, 0, 0)]);
            }
            if (isInflux) {
                inFlux += (speed*CalcBGKFeq(xiIndex, 1, u, v, T, FEQORDER));
            }
        }
        Real rho = outFlux / inFlux;
        //ops_printf("vg=%i rho=%f outflux=%f
influx=%f\n",vg,rho,outFlux,inFlux); for (int xiIndex = 0; xiIndex < NUMXI;
xiIndex++) { const Real relVeloX = CS * XI[xiIndex * LATTDIM] - u; const Real
relVeloY = CS * XI[xiIndex * LATTDIM + 1] - v; Real cDotPrimary = relVeloX *
primaryVector[0] + relVeloY * primaryVector[1]; bool isInflux = cDotPrimary >=
0; if (0 == boundaryType) { Real cDotSecond = relVeloX * secondVector[0] +
relVeloY * secondVector[1]; isInflux = isInflux && (cDotSecond >= 0);
            }
            if (1 == boundaryType) {
                Real cDotSecond =
                    relVeloX * secondVector[0] + relVeloY * secondVector[1];
                isInflux = isInflux || (cDotSecond > 0);
            }
            if (isInflux) {
                f[OPS_ACC_MD3(xiIndex, 0, 0)] =
                    CalcBGKFeq(xiIndex, rho, u, v, T, FEQORDER);
                ;
            }
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a kinetic boundary "
                   "point: KerCutCellKinetic");
#endif
    }
}
*/

void KerCutCellCorrectedKinetic(const Real *givenMacroVars, const Real *dt,
                                const int *nodeType,
                                const int *geometryProperty, const Real *tau,
                                const Real *feq, Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC2(0, 0)];
    if (vt == Vertex_KineticDiffuseWall) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC3(0, 0)];
        const Real u = givenMacroVars[1];
        const Real v = givenMacroVars[2];
        const Real kn =
            tau[OPS_ACC_MD4(0, 0, 0)];  // only for single componment
        Real wallNormalVector[]{0, 0};
        const Real sqrt2Inverse = 1 / sqrt(2);
        switch (vg) {
            case VG_IP: {
                wallNormalVector[0] = 1;
                wallNormalVector[1] = 0;
            } break;
            case VG_IM: {
                wallNormalVector[0] = -1;
                wallNormalVector[1] = 0;
            } break;
            case VG_JP: {
                wallNormalVector[0] = 0;
                wallNormalVector[1] = 1;
            } break;
            case VG_JM: {
                wallNormalVector[0] = 0;
                wallNormalVector[1] = -1;
            } break;
            case VG_IPJP_I:  // inner corner point
            {
                wallNormalVector[0] = sqrt2Inverse;
                wallNormalVector[1] = sqrt2Inverse;
            } break;
            case VG_IPJM_I:  // inner corner point
            {
                wallNormalVector[0] = sqrt2Inverse;
                wallNormalVector[1] = -sqrt2Inverse;
            } break;
            case VG_IMJP_I:  // inner corner point
            {
                wallNormalVector[0] = -sqrt2Inverse;
                wallNormalVector[1] = sqrt2Inverse;
            } break;
            case VG_IMJM_I:  // inner corner point
            {
                wallNormalVector[0] = -sqrt2Inverse;
                wallNormalVector[1] = -sqrt2Inverse;
            } break;
            default:
                break;
        }
        Real outFlux = 0;  // flow into wall
        Real inFlux = 0;   // flow into fluid bulk
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            const Real cx = XI[xiIndex * LATTDIM];
            const Real cy = XI[xiIndex * LATTDIM + 1];
            Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
                              (CS * cy - v) * wallNormalVector[1];
            if (cDotNormal < 0) {
                outFlux +=
                    (-cDotNormal *
                     ((f[OPS_ACC_MD6(xiIndex, 0, 0)] +
                       (*dt) * feq[OPS_ACC_MD5(xiIndex, 0, 0)] / (2 * kn)) /
                      (1 + (*dt) / (2 * kn))));
            }
            if (cDotNormal > 0) {
                Real cu = (CS * cx * u + CS * cy * v);
                inFlux += (cDotNormal * WEIGHTS[xiIndex] *
                           (1 + cu +
                            0.5 * (cu * cu - (u * u + v * v))));  // i.e., the
                // equilibrium
            }
        }
        Real rho = outFlux / inFlux;
        for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            const Real cx = XI[xiIndex * LATTDIM];
            const Real cy = (int)XI[xiIndex * LATTDIM + 1];
            Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
                              (CS * cy - v) * wallNormalVector[1];
            if (cDotNormal >= 0) {
                Real cu = (CS * cx * u + CS * cy * v);
                f[OPS_ACC_MD6(xiIndex, 0, 0)] =
                    (rho * WEIGHTS[xiIndex] *
                     (1 + cu + 0.5 * (cu * cu - (u * u + v * v))));
            }
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a kinetic boundary "
                   "point: KerCutCellKinetic");
#endif
    }
}

void KerCutCellBounceBack(const int *nodeType, const int *geometryProperty,
                          Real *f) {
    /*!
     We consider zero velocity boundary first
     To make sure the velocity at boundary is zero, the implementation
     is lattice specific.
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    if (vt == Vertex_BounceBackWall) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
        switch (vg) {
            case VG_IP: {
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
                f[OPS_ACC_MD2(1, 0, 0)] = f[OPS_ACC_MD2(3, 0, 0)];
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
            } break;
            case VG_IM: {
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
                f[OPS_ACC_MD2(3, 0, 0)] = f[OPS_ACC_MD2(1, 0, 0)];
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
            } break;
            case VG_JP: {
                f[OPS_ACC_MD2(2, 0, 0)] = f[OPS_ACC_MD2(4, 0, 0)];
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
            } break;
            case VG_JM: {
                f[OPS_ACC_MD2(4, 0, 0)] = f[OPS_ACC_MD2(2, 0, 0)];
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
            } break;
            case VG_IPJP_I:  // inner corner point
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
                break;
            case VG_IPJM_I:  // inner corner point
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
                break;
            case VG_IMJP_I:  // inner corner point
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
                break;
            case VG_IMJM_I:  // inner corner point
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
                break;
            default:
                break;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a bounce-back boundary "
                   "point: KerCutCellBounceBack");
#endif
    }
}

void KerCutCellEQMDiffuseRefl(const Real *givenMacroVars, const int *nodeType,
                              const int *geometryProperty, Real *f,
                              const int *componentId) {
    // This kernel is suitable for a single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    const int equilibriumOrder{2};
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_EQMDiffuseRefl) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        Real u = givenMacroVars[1];
        Real v = givenMacroVars[2];
        // for (int compoIdx = 0; compoIdx < NUMCOMPONENTS; compoIdx++) {
        const int compoIdx{*componentId};
        int numOutgoing{0};
        int numIncoming{0};
        int numParallel{0};
        int *outgoing = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        int *incoming = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        int *parallel = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        Real rhoIncoming{0};
        Real rhoParallel{0};
        Real deltaRho{0};
        for (int xiIdx = COMPOINDEX[2 * compoIdx];
             xiIdx <= COMPOINDEX[2 * compoIdx + 1]; xiIdx++) {
            BndryDvType bdt = FindBdyDvType(vg, &XI[xiIdx * LATTDIM]);
            switch (bdt) {
                case BndryDv_Incoming: {
                    incoming[numIncoming] = xiIdx;
                    rhoIncoming += f[OPS_ACC_MD3(xiIdx, 0, 0)];
                    numIncoming++;
                } break;
                case BndryDv_Outgoing: {
                    outgoing[numOutgoing] = xiIdx;
                    Real cx{CS * XI[xiIdx * LATTDIM]};
                    Real cy{CS * XI[xiIdx * LATTDIM + 1]};
                    deltaRho += (2 * WEIGHTS[xiIdx]) * (cx * u + cy * v);
                    numOutgoing++;
                } break;
                case BndryDv_Parallel: {
                    parallel[numParallel] = xiIdx;
                    rhoParallel +=
                        CalcBGKFeq(xiIdx, 1, u, v, 1, equilibriumOrder);
                    numParallel++;
                } break;
                default:
                    break;
            }
        }
        Real rhoWall = 2 * rhoIncoming / (1 - deltaRho - rhoParallel);
        for (int idx = 0; idx < numParallel; idx++) {
            f[OPS_ACC_MD3(parallel[idx], 0, 0)] =
                CalcBGKFeq(parallel[idx], rhoWall, u, v, 1, equilibriumOrder);
        }
        for (int idx = 0; idx < numOutgoing; idx++) {
            int xiIdx = outgoing[idx];
            Real cx{CS * XI[xiIdx * LATTDIM]};
            Real cy{CS * XI[xiIdx * LATTDIM + 1]};
            f[OPS_ACC_MD3(xiIdx, 0, 0)] =
                f[OPS_ACC_MD3(OPP[xiIdx], 0, 0)] +
                2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
        }
        delete[] outgoing;
        delete[] incoming;
        delete[] parallel;
        //}
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a equilibrium diffuse reflection "
                   "boundary condition point: KerCutCellEQMDiffuseRefl");
#endif
    }
}

void KerCutCellBounceBackNew(const int *nodeType, const int *geometryProperty,
                             Real *f) {
    /*!
     We consider zero velocity boundary first
     To make sure the velocity at boundary is zero, the implementation
     is lattice specific.
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    if (vt == Vertex_BounceBackWall) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
        switch (vg) {
            case VG_IP: {
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
                f[OPS_ACC_MD2(1, 0, 0)] = f[OPS_ACC_MD2(3, 0, 0)];
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
            } break;
            case VG_IM: {
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
                f[OPS_ACC_MD2(3, 0, 0)] = f[OPS_ACC_MD2(1, 0, 0)];
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
            } break;
            case VG_JP: {
                f[OPS_ACC_MD2(2, 0, 0)] = f[OPS_ACC_MD2(4, 0, 0)];
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
                const Real rhow =
                    6 * (f[OPS_ACC_MD2(4, 0, 0)] + f[OPS_ACC_MD2(8, 0, 0)] +
                         f[OPS_ACC_MD2(7, 0, 0)]);
                f[OPS_ACC_MD2(1, 0, 0)] = rhow / 9;
                f[OPS_ACC_MD2(3, 0, 0)] = rhow / 9;
                f[OPS_ACC_MD2(0, 0, 0)] = 4 * rhow / 9;
            } break;
            case VG_JM: {
                f[OPS_ACC_MD2(4, 0, 0)] = f[OPS_ACC_MD2(2, 0, 0)];
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
                const Real rhow =
                    6 * (f[OPS_ACC_MD2(2, 0, 0)] + f[OPS_ACC_MD2(6, 0, 0)] +
                         f[OPS_ACC_MD2(5, 0, 0)]);
                f[OPS_ACC_MD2(1, 0, 0)] = rhow / 9;
                f[OPS_ACC_MD2(3, 0, 0)] = rhow / 9;
                f[OPS_ACC_MD2(0, 0, 0)] = 4 * rhow / 9;
            } break;
            case VG_IPJP_I:  // inner corner point
                f[OPS_ACC_MD2(5, 0, 0)] = f[OPS_ACC_MD2(7, 0, 0)];
                break;
            case VG_IPJM_I:  // inner corner point
                f[OPS_ACC_MD2(8, 0, 0)] = f[OPS_ACC_MD2(6, 0, 0)];
                break;
            case VG_IMJP_I:  // inner corner point
                f[OPS_ACC_MD2(6, 0, 0)] = f[OPS_ACC_MD2(8, 0, 0)];
                break;
            case VG_IMJM_I:  // inner corner point
                f[OPS_ACC_MD2(7, 0, 0)] = f[OPS_ACC_MD2(5, 0, 0)];
                break;
            default:
                break;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a bounce-back boundary "
                   "point: KerCutCellBounceBack");
#endif
    }
}

void KerCutCellPeriodic(const int *nodeType, const int *geometryProperty,
                        Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    if (vt == Vertex_Periodic) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0)];
        switch (vg) {
            case VG_IP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    Real cx = XI[xiIndex * LATTDIM];
                    if (cx > 0) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0)];
                    }
                }
                break;
            case VG_IM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    Real cx = XI[xiIndex * LATTDIM];
                    if (cx < 0) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0)];
                    }
                }
                break;
            case VG_JP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    Real cy = XI[xiIndex * LATTDIM + 1];
                    if (cy > 0) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1)];
                    }
                }
                break;
            case VG_JM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    Real cy = XI[xiIndex * LATTDIM + 1];
                    if (cy < 0) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1)];
                    }
                }
                break;
            // as we are dealing with domain boundary, there can be only inner
            // corner
            case VG_IPJP_I: {                          // corner point
                if (vt == nodeType[OPS_ACC0(0, 1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy > 0 || cx > 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, -1, 0)];
                        }
                    }
                }
                if (vt == nodeType[OPS_ACC0(1, 0)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy > 0 || cx > 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 0, -1)];
                        }
                    }
                }
            } break;
            case VG_IPJM_I: {                           // inner corner point
                if (vt == nodeType[OPS_ACC1(0, -1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy < 0 || cx > 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, -1, 0)];
                        }
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy < 0 || cx > 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 0, 1)];
                        }
                    }
                }
            } break;
            case VG_IMJP_I: {                          // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy > 0 || cx < 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 1, 0)];
                        }
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy > 0 || cx < 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 0, -1)];
                        }
                    }
                }
            } break;
            case VG_IMJM_I: {
                // inner corner point
                if (vt == nodeType[OPS_ACC1(0, -1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy < 0 || cx < 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 1, 0)];
                        }
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        Real cx = XI[xiIndex * LATTDIM];
                        Real cy = XI[xiIndex * LATTDIM + 1];
                        if (cy < 0 || cx < 0) {
                            f[OPS_ACC_MD2(xiIndex, 0, 0)] =
                                f[OPS_ACC_MD2(xiIndex, 0, 1)];
                        }
                    }
                }
            } break;
            default:
                break;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a periodic boundary "
                   "point: KerCutCellPeriodic");
#endif
    }
}

void KerEquibriumVelocity(const Real *givenMacroVars, Real *f) {
    Real rho = 1;
    Real u = givenMacroVars[1];
    Real v = givenMacroVars[2];
    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
        const int cx = (int)XI[xiIndex * LATTDIM];
        const int cy = (int)XI[xiIndex * LATTDIM + 1];
        const Real cu = (CS * cx * u + CS * cy * v);
        f[OPS_ACC_MD1(xiIndex, 0, 0)] =
            WEIGHTS[xiIndex] * rho *
            (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
    }
}

void KerCutCellZouHeVelocity(const Real *givenMacroVars, const int *nodeType,
                             const int *geometryProperty, const Real *macroVars,
                             Real *f) {
    /*!
    Note: This boundary condition requires both stream and collision happenning
    at a boundary point.
    Note: This boundary condition is lattice specific.
    */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_ZouHeVelocity) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        Real rho{0};
        Real u{givenMacroVars[1]};
        Real v{givenMacroVars[2]};
        Real sqrt3 = sqrt(3);
        switch (vg) {
            case VG_IP: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f2 = f[OPS_ACC_MD4(2, 0, 0)];
                Real f3 = f[OPS_ACC_MD4(3, 0, 0)];
                Real f4 = f[OPS_ACC_MD4(4, 0, 0)];
                Real f6 = f[OPS_ACC_MD4(6, 0, 0)];
                Real f7 = f[OPS_ACC_MD4(7, 0, 0)];
                rho = sqrt3 * (f0 + f2 + 2 * f3 + f4 + 2 * f6 + 2 * f7) /
                      (sqrt3 - u);
                f[OPS_ACC_MD4(1, 0, 0)] = (2 * sqrt3 * rho * u + 9 * f3) / 9.0;
                f[OPS_ACC_MD4(5, 0, 0)] =
                    (sqrt3 * rho * u + 3 * sqrt3 * rho * v - 9 * f2 + 9 * f4 +
                     18 * f7) /
                    18.0;
                f[OPS_ACC_MD4(8, 0, 0)] =
                    (sqrt3 * rho * u - 3 * sqrt3 * rho * v + 9 * f2 - 9 * f4 +
                     18 * f6) /
                    18.0;
            } break;
            case VG_IM: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f2 = f[OPS_ACC_MD4(2, 0, 0)];
                Real f4 = f[OPS_ACC_MD4(4, 0, 0)];
                Real f1 = f[OPS_ACC_MD4(1, 0, 0)];
                Real f5 = f[OPS_ACC_MD4(5, 0, 0)];
                Real f8 = f[OPS_ACC_MD4(8, 0, 0)];
                rho = (sqrt3 * f0 + 2 * sqrt3 * f1 + sqrt3 * f2 + sqrt3 * f4 +
                       2 * sqrt3 * f5 + 2 * sqrt3 * f8) /
                      (sqrt3 + u);
                f[OPS_ACC_MD4(3, 0, 0)] = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
                f[OPS_ACC_MD4(6, 0, 0)] =
                    (-(sqrt3 * u * rho) + 3 * sqrt3 * v * rho - 9 * f2 +
                     9 * f4 + 18 * f8) /
                    18.0;
                f[OPS_ACC_MD4(7, 0, 0)] =
                    (-(sqrt3 * u * rho) - 3 * sqrt3 * v * rho + 9 * f2 -
                     9 * f4 + 18 * f5) /
                    18.0;
            } break;
            case VG_JP: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f1 = f[OPS_ACC_MD4(1, 0, 0)];
                Real f3 = f[OPS_ACC_MD4(3, 0, 0)];
                Real f4 = f[OPS_ACC_MD4(4, 0, 0)];
                Real f7 = f[OPS_ACC_MD4(7, 0, 0)];
                Real f8 = f[OPS_ACC_MD4(8, 0, 0)];
                rho = (sqrt3 * f0 + sqrt3 * f1 + sqrt3 * f3 + 2 * sqrt3 * f4 +
                       2 * sqrt3 * f7 + 2 * sqrt3 * f8) /
                      (sqrt3 - v);
                f[OPS_ACC_MD4(2, 0, 0)] = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
                f[OPS_ACC_MD4(5, 0, 0)] =
                    (3 * sqrt3 * u * rho + sqrt3 * v * rho - 9 * f1 + 9 * f3 +
                     18 * f7) /
                    18.0;
                f[OPS_ACC_MD4(6, 0, 0)] =
                    (-3 * sqrt3 * u * rho + sqrt3 * v * rho + 9 * f1 - 9 * f3 +
                     18 * f8) /
                    18.0;
            } break;
            case VG_JM: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f1 = f[OPS_ACC_MD4(1, 0, 0)];
                Real f3 = f[OPS_ACC_MD4(3, 0, 0)];
                Real f2 = f[OPS_ACC_MD4(2, 0, 0)];
                Real f5 = f[OPS_ACC_MD4(5, 0, 0)];
                Real f6 = f[OPS_ACC_MD4(6, 0, 0)];
                rho = (sqrt3 * f0 + sqrt3 * f1 + 2 * sqrt3 * f2 + sqrt3 * f3 +
                       2 * sqrt3 * f5 + 2 * sqrt3 * f6) /
                      (sqrt3 + v);
                f[OPS_ACC_MD4(4, 0, 0)] = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
                f[OPS_ACC_MD4(7, 0, 0)] =
                    (-3 * sqrt3 * u * rho - sqrt3 * v * rho + 9 * f1 - 9 * f3 +
                     18 * f5) /
                    18.0;
                f[OPS_ACC_MD4(8, 0, 0)] =
                    (3 * sqrt3 * u * rho - sqrt3 * v * rho - 9 * f1 + 9 * f3 +
                     18 * f6) /
                    18.0;
            } break;
            case VG_IPJM_I: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f2 = f[OPS_ACC_MD4(2, 0, 0)];
                Real f6 = f[OPS_ACC_MD4(6, 0, 0)];
                Real f3 = f[OPS_ACC_MD4(3, 0, 0)];
                rho = macroVars[OPS_ACC_MD3(0, 1, -1)];
                f[OPS_ACC_MD4(1, 0, 0)] = (2 * sqrt3 * u * rho + 9 * f3) / 9.0;
                f[OPS_ACC_MD4(5, 0, 0)] =
                    (9 * rho - 2 * sqrt3 * u * rho + 3 * sqrt3 * v * rho -
                     9 * f0 - 18 * f2 - 18 * f3 - 18 * f6) /
                    18.0;
                f[OPS_ACC_MD4(7, 0, 0)] =
                    (9 * rho - 3 * sqrt3 * u * rho + 2 * sqrt3 * v * rho -
                     9 * f0 - 18 * f2 - 18 * f3 - 18 * f6) /
                    18.0;
                f[OPS_ACC_MD4(4, 0, 0)] = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
                f[OPS_ACC_MD4(8, 0, 0)] =
                    (sqrt3 * u * rho - sqrt3 * v * rho + 18 * f6) / 18.0;
            } break;
            case VG_IPJP_I: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f4 = f[OPS_ACC_MD4(4, 0, 0)];
                Real f3 = f[OPS_ACC_MD4(3, 0, 0)];
                Real f7 = f[OPS_ACC_MD4(7, 0, 0)];
                rho = macroVars[OPS_ACC_MD3(0, 1, 1)];
                f[OPS_ACC_MD4(1, 0, 0)] = (2 * sqrt3 * u * rho + 9 * f3) / 9.0;
                f[OPS_ACC_MD4(5, 0, 0)] =
                    (sqrt3 * u * rho + sqrt3 * v * rho + 18 * f7) / 18.0;
                f[OPS_ACC_MD4(8, 0, 0)] =
                    (9 * rho - 2 * sqrt3 * u * rho - 3 * sqrt3 * v * rho -
                     9 * f0 - 18 * f3 - 18 * f4 - 18 * f7) /
                    18.0;
                f[OPS_ACC_MD4(2, 0, 0)] = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
                f[OPS_ACC_MD4(6, 0, 0)] =
                    (9 * rho - 3 * sqrt3 * u * rho - 2 * sqrt3 * v * rho -
                     9 * f0 - 18 * f3 - 18 * f4 - 18 * f7) /
                    18.0;
            } break;
            case VG_IMJP_I: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f1 = f[OPS_ACC_MD4(1, 0, 0)];
                Real f8 = f[OPS_ACC_MD4(8, 0, 0)];
                Real f4 = f[OPS_ACC_MD4(4, 0, 0)];
                rho = macroVars[OPS_ACC_MD3(0, -1, 1)];
                f[OPS_ACC_MD4(5, 0, 0)] =
                    (9 * rho + 3 * sqrt3 * u * rho - 2 * sqrt3 * v * rho -
                     9 * f0 - 18 * f1 - 18 * f4 - 18 * f8) /
                    18.0;
                f[OPS_ACC_MD4(2, 0, 0)] = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
                f[OPS_ACC_MD4(6, 0, 0)] =
                    (-(sqrt3 * u * rho) + sqrt3 * v * rho + 18 * f8) / 18.0;
                f[OPS_ACC_MD4(3, 0, 0)] = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
                f[OPS_ACC_MD4(7, 0, 0)] =
                    (9 * rho + 2 * sqrt3 * u * rho - 3 * sqrt3 * v * rho -
                     9 * f0 - 18 * f1 - 18 * f4 - 18 * f8) /
                    18.0;
            } break;
            case VG_IMJM_I: {
                // Knows
                Real f0 = f[OPS_ACC_MD4(0, 0, 0)];
                Real f1 = f[OPS_ACC_MD4(1, 0, 0)];
                Real f2 = f[OPS_ACC_MD4(2, 0, 0)];
                Real f5 = f[OPS_ACC_MD4(5, 0, 0)];
                rho = macroVars[OPS_ACC_MD3(0, -1, -1)];
                f[OPS_ACC_MD4(6, 0, 0)] =
                    (9 * rho + 2 * sqrt3 * u * rho + 3 * sqrt3 * v * rho -
                     9 * f0 - 18 * f1 - 18 * f2 - 18 * f5) /
                    18.0;
                f[OPS_ACC_MD4(3, 0, 0)] = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
                f[OPS_ACC_MD4(7, 0, 0)] =
                    (-(sqrt3 * u * rho) - sqrt3 * v * rho + 18 * f5) / 18.0;
                f[OPS_ACC_MD4(4, 0, 0)] = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
                f[OPS_ACC_MD4(8, 0, 0)] =
                    (9 * rho + 3 * sqrt3 * u * rho + 2 * sqrt3 * v * rho -
                     9 * f0 - 18 * f1 - 18 * f2 - 18 * f5) /
                    18.0;
            } break;
            default:
                break;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a Zou-He velocity boundary "
                   "point: KerCutCellZhouHeVelocity");
#endif
    }
}

void KerCutCellNonEqExtrapol(const Real *givenMacroVars, const int *nodeType,
                             const int *geometryProperty, const Real *macroVars,
                             const Real *feq, Real *f, const int componentID) {
    /*!
     Note: Here we are implementing the version defined in the book "Lattice
     Boltzmann Method and Its Applications in Engineering" by Guo and Shu.
     This version is different from the original verion presented in their
     paper. Note: The implementation can only be used for singe-speed lattice,
     e.g., D2Q9
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_NoneqExtrapol) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        const Real u{givenMacroVars[1]};
        const Real v{givenMacroVars[2]};
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            for (int xiIndex = COMPOINDEX[2 * compoIndex];
                 xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                switch (vg) {
                    case VG_IP: {
                        if (cx > 0) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                    } break;
                    case VG_IM: {
                        if (cx < 0) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                    } break;
                    case VG_JP: {
                        if (cy > 0) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (cy==0) {
                        //                            const Real
                        //                            rho{macroVars[OPS_ACC_MD3(0,
                        //                            0,
                        //                            1)]/*-macroVars[OPS_ACC_MD3(0,
                        //                            0, 2)]*/};
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)]=WEIGHTS[xiIndex] *
                        //                            rho;
                        //                        }
                    } break;
                    case VG_JM: {
                        if (cy < 0) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (cy==0) {
                        //                            const Real
                        //                            rho{macroVars[OPS_ACC_MD3(0,
                        //                            0,
                        //                            -1)]/*-macroVars[OPS_ACC_MD3(0,
                        //                            0, -2)]*/};
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)]=WEIGHTS[xiIndex] *
                        //                            rho;
                        //                        }

                    } break;
                    case VG_IPJP_I: {  // inner corner point
                        if (cy >= 0 && cx >= 0 && (cx != 0 || cy != 0)) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            // const Real rho = 1.1;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                       Real rhotmp = (f[OPS_ACC_MD2(0,
                        //                       0, 0)] +
                        //                                      f[OPS_ACC_MD2(6,
                        //                                      0, 0)] +
                        //                                      f[OPS_ACC_MD2(8,
                        //                                      0, 0)]);
                        //                        const Real rho = 1.1;
                        //                        f[OPS_ACC_MD2(0, 0, 0)] =
                        //                        4*rho/9; f[OPS_ACC_MD2(6, 0,
                        //                        0)] = rho/36; f[OPS_ACC_MD2(8,
                        //                        0, 0)] = rho/36;
                    } break;
                    case VG_IPJM_I: {  // inner corner point·
                        if (cy <= 0 && cx >= 0 && (cx != 0 || cy != 0)) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            // const Real rho = 1.1;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                       Real rhotmp = (f[OPS_ACC_MD2(0,
                        //                       0, 0)] +
                        //                                      f[OPS_ACC_MD2(5,
                        //                                      0, 0)] +
                        //                                      f[OPS_ACC_MD2(7,
                        //                                      0, 0)]);
                        //                        const Real rho = 1.1;
                        //                        f[OPS_ACC_MD2(0, 0, 0)] =
                        //                        4*rho/9; f[OPS_ACC_MD2(5, 0,
                        //                        0)] = rho/36; f[OPS_ACC_MD2(7,
                        //                        0, 0)] = rho/36;
                    } break;
                    case VG_IMJP_I: {  // inner corner point
                        if (cy >= 0 && cx <= 0 && (cx != 0 || cy != 0)) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            // const Real rho = 1;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                       Real rhotmp = (f[OPS_ACC_MD2(0,
                        //                       0, 0)] +
                        //                                      f[OPS_ACC_MD2(5,
                        //                                      0, 0)] +
                        //                                      f[OPS_ACC_MD2(7,
                        //                                      0, 0)]);
                        //                        const Real rho = 1;
                        //                        f[OPS_ACC_MD2(0, 0, 0)] =
                        //                        4*rho/9; f[OPS_ACC_MD2(5, 0,
                        //                        0)] = rho/36; f[OPS_ACC_MD2(7,
                        //                        0, 0)] = rho/36;
                    } break;
                    case VG_IMJM_I: {  // inner corner point
                        if (cy <= 0 && cx <= 0 && (cx != 0 || cy != 0)) {
                            const Real rho{macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            // const Real rho = 1;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rho *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                       Real rhotmp = (f[OPS_ACC_MD2(0,
                        //                       0, 0)] +
                        //                                      f[OPS_ACC_MD2(6,
                        //                                      0, 0)] +
                        //                                      f[OPS_ACC_MD2(8,
                        //                                      0, 0)]);
                        //                        const Real rho=1;
                        //                        f[OPS_ACC_MD2(0, 0, 0)] =
                        //                        4*rho/9; f[OPS_ACC_MD2(6, 0,
                        //                        0)] = rho/36; f[OPS_ACC_MD2(8,
                        //                        0, 0)] = rho/36;
                    } break;
                    default:
                        break;
                }
            }
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a non-equilibrium extraploation "
                   "boundary "
                   "point: KerCutCellNoneqExtrapol");
#endif
    }
}

void KerCutCellNonEqExtrapolPressure(const Real *givenMacroVars,
                                     const int *nodeType,
                                     const int *geometryProperty,
                                     const Real *macroVars, const Real *feq,
                                     Real *f) {
    /*!
     Note: Here we are implementing the version defined in the book "Lattice
     Boltzmann Method and Its Applications in Engineering" by Guo and Shu.
     Thie version is different from the original paper.
     Note: The implementation can only be used for singe-speed lattice, e.g.,
     D2Q9
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0)];
    if (vt == Vertex_NonEqExtrapolPressure) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0)];
        const Real rhoGiven{givenMacroVars[0]};
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            Real rhoLocal = 0;
            for (int xiIndex = COMPOINDEX[2 * compoIndex];
                 xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                switch (vg) {
                    case VG_IP: {
                        if (cx > 0) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            // const Real v =0;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                            // printf("%s%f%s%f\n","u=",u," v=",v);
                        }
                        //                        if (cx == 0) {
                        //                            const Real
                        //                            rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                            1,
                        //                            0)]-macroVars[OPS_ACC_MD3(0,
                        //                            2, 0)]}; const Real
                        //                            u=2*(macroVars[OPS_ACC_MD3(1,
                        //                            1,
                        //                            0)]/macroVars[OPS_ACC_MD3(0,
                        //                            1, 0)])
                        //                            -macroVars[OPS_ACC_MD3(1,
                        //                            2,
                        //                            0)]/macroVars[OPS_ACC_MD3(0,
                        //                            2, 0)]; const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                            WEIGHTS[xiIndex] * rhoNext
                        //                            * (1 + cu + 0.5 * (cu * cu
                        //                            - (u * u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                            //printf("%s%f%s%f\n","u=",u,"
                        //                            v=",v);
                        //                        }
                    } break;
                    case VG_IM: {
                        if (cx < 0) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            // const Real v =0;
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (cx == 0) {
                        //                            const Real
                        //                            rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                            -1,
                        //                            0)]-macroVars[OPS_ACC_MD3(0,
                        //                            -2, 0)]}; const Real
                        //                            u=2*(macroVars[OPS_ACC_MD3(1,
                        //                            -1,
                        //                            0)]/macroVars[OPS_ACC_MD3(0,
                        //                            -1, 0)])
                        //                            -macroVars[OPS_ACC_MD3(1,
                        //                            -2,
                        //                            0)]/macroVars[OPS_ACC_MD3(0,
                        //                            -2, 0)]; const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                            WEIGHTS[xiIndex] *
                        //                            0.5*(rhoNext+rhoGiven) *
                        //                            (1 + cu + 0.5 * (cu * cu -
                        //                            (u * u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                            //printf("%s%f%s%f\n","u=",u,"
                        //                            v=",v);
                        //                        }
                    } break;
                    case VG_JP: {
                        if (cy > 0) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                    } break;
                    case VG_JM: {
                        if (cy < 0) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                    } break;
                    case VG_IPJP_I: {  // inner corner point
                        if (cy >= 0 && cx >= 0 && (cx != 0 || cy != 0)) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (
                        //                            (cx == 0 && cy == 0) ||
                        //                            (cx == 1 && cy == -1) ||
                        //                            (cx == -1 && cy == 1)
                        //                            ) {
                        //                            const Real
                        //                            rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                            1,
                        //                            0)]-macroVars[OPS_ACC_MD3(0,
                        //                            2, 0)]}; const Real
                        //                            u{macroVars[OPS_ACC_MD3(1,
                        //                            1, 0)] /
                        //                                         rhoNext};
                        //                            const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                                WEIGHTS[xiIndex] *
                        //                                rhoGiven * (1 + cu +
                        //                                0.5 * (cu * cu - (u *
                        //                                u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                        }
                        // printf("%s%f\n","f0=",f[OPS_ACC_MD5(0,0,0)]);
                    } break;
                    case VG_IPJM_I: {  // inner corner point·
                        if (cy <= 0 && cx >= 0 && (cx != 0 || cy != 0)) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (
                        //                            (cx == 0 && cy == 0) ||
                        //                            (cx == 1 && cy == 1) ||
                        //                            (cx == -1 && cy == -1)
                        //                            ) {
                        //                            const Real
                        //                            rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                            1,
                        //                            0)]-macroVars[OPS_ACC_MD3(0,
                        //                            2, 0)]}; const Real
                        //                            u{macroVars[OPS_ACC_MD3(1,
                        //                            1, 0)] /
                        //                                         rhoNext};
                        //                            const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                                WEIGHTS[xiIndex] *
                        //                                rhoGiven * (1 + cu +
                        //                                0.5 * (cu * cu - (u *
                        //                                u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                        }

                    } break;
                    case VG_IMJP_I: {  // inner corner point
                        if (cy >= 0 && cx <= 0 && (cx != 0 || cy != 0)) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (
                        //                            (cx == 0 && cy == 0) ||
                        //                            (cx == 1 && cy == 1) ||
                        //                            (cx == -1 && cy == -1)
                        //                            ) {
                        //                            const Real
                        //                            rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                            -1,
                        //                            0)]-macroVars[OPS_ACC_MD3(0,
                        //                            -2, 0)]}; const Real
                        //                            u{macroVars[OPS_ACC_MD3(1,
                        //                            -1, 0)] /
                        //                                         rhoNext};
                        //                            const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                                WEIGHTS[xiIndex] *
                        //                                rhoGiven * (1 + cu +
                        //                                0.5 * (cu * cu - (u *
                        //                                u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                        }
                    } break;
                    case VG_IMJM_I: {  // inner corner point
                        if (cy <= 0 && cx <= 0 && (cx != 0 || cy != 0)) {
                            const Real rhoNext{
                                macroVars[OPS_ACC_MD3(0, cx, cy)]};
                            const Real u{macroVars[OPS_ACC_MD3(1, cx, cy)] /
                                         rhoNext};
                            const Real v{macroVars[OPS_ACC_MD3(2, cx, cy)] /
                                         rhoNext};
                            const Real cu = (CS * cx * u + CS * cy * v);
                            const Real eqPart =
                                WEIGHTS[xiIndex] * rhoGiven *
                                (1 + cu + 0.5 * (cu * cu - (u * u + v * v)));
                            const Real noneqPart =
                                (f[OPS_ACC_MD5(xiIndex, cx, cy)] -
                                 feq[OPS_ACC_MD4(xiIndex, cx, cy)]);
                            f[OPS_ACC_MD5(xiIndex, 0, 0)] = eqPart + noneqPart;
                        }
                        //                        if (
                        //                            (cx == 0 && cy == 0) ||
                        //                            (cx == 1 && cy == -1) ||
                        //                            (cx == -1 && cy == 1)
                        //                            ) {
                        //                           const Real
                        //                           rhoNext{2*macroVars[OPS_ACC_MD3(0,
                        //                           -1,
                        //                           0)]-macroVars[OPS_ACC_MD3(0,
                        //                           -2, 0)]};
                        //                            const Real
                        //                            u{macroVars[OPS_ACC_MD3(1,
                        //                            -1, 0)] /
                        //                                         rhoNext};
                        //                            const Real v{0};
                        //                            const Real cu = (CS * cx *
                        //                            u + CS * cy * v); const
                        //                            Real eqPart =
                        //                                WEIGHTS[xiIndex] *
                        //                                rhoGiven * (1 + cu +
                        //                                0.5 * (cu * cu - (u *
                        //                                u + v * v)));
                        //                            f[OPS_ACC_MD5(xiIndex, 0,
                        //                            0)] = eqPart;
                        //                        }
                    } break;
                    default:
                        break;
                }
                rhoLocal += f[OPS_ACC_MD5(xiIndex, 0, 0)];
            }
            //            Real ratio = rhoGiven / rhoLocal;
            // printf("%s%f\n","rhoLocal=",rhoLocal);
            //    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
            //        f[OPS_ACC_MD5(xiIndex, 0, 0)] *= ratio;
            //     }
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a non-equilibrium extraploation "
                   "boundary "
                   "point: KerCutCellNoneqExtrapol");
#endif
    }
}
#endif
#ifdef OPS_3D
void KerCutCellExtrapolPressure1ST3D(const Real *givenBoundaryVars,
                                     const int *nodeType,
                                     const int *geometryProperty, Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0, 0)];
    if (vt == Vertex_ExtrapolPressure1ST) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0, 0)];
        Real rhoGiven = givenBoundaryVars[0];
        Real rho = 0;
        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
            Real cx{CS * XI[xiIdx * LATTDIM]};
            Real cy{CS * XI[xiIdx * LATTDIM + 1]};
            Real cz{CS * XI[xiIdx * LATTDIM + 2]};
            switch (vg) {
                case VG_IP: {
                    if (cx > 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                    }
                } break;
                case VG_IM: {
                    if (cx < 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                    }
                } break;
                case VG_JP: {
                    if (cy > 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                    }
                } break;
                case VG_JM: {
                    if (cy < 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                    }
                } break;
                case VG_KP: {
                    if (cz > 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                    }
                } break;
                case VG_KM: {
                    if (cz < 0) {
                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                            f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                    }
                } break;
                case VG_IPJP_I: {
                    if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                    }
                } break;
                case VG_IPJM_I: {
                    if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                    }
                } break;
                case VG_IMJP_I: {
                    if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                    }
                } break;
                case VG_IMJM_I: {
                    if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                    }
                } break;
                case VG_IPKP_I: {
                    if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IPKM_I: {
                    if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_IMKP_I: {
                    if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IMKM_I: {
                    if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_JPKP_I: {
                    if ((cy >= 0 && cz > 0) || (cy > 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_JPKM_I: {
                    if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_JMKP_I: {
                    if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_JMKM_I: {
                    if ((cy <= 0 && cz < 0) || (cy < 0 && cz == 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_IPJPKP_I: {
                    if ((cx >= 0 && cy >= 0 && cz >= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IPJPKM_I: {
                    if ((cx >= 0 && cy >= 0 && cz <= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_IPJMKP_I: {
                    if ((cx >= 0 && cy <= 0 && cz >= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IPJMKM_I: {
                    if ((cx >= 0 && cy <= 0 && cz <= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(1, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_IMJPKP_I: {
                    if ((cx <= 0 && cy >= 0 && cz >= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IMJPKM_I: {
                    if ((cx <= 0 && cy >= 0 && cz <= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, 1, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                case VG_IMJMKP_I: {
                    if ((cx <= 0 && cy <= 0 && cz >= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, 1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 1)];
                        }
                    }
                } break;
                case VG_IMJMKM_I: {
                    if ((cx <= 0 && cy <= 0 && cz <= 0) &&
                        (cx != 0 || cy != 0 || cz != 0)) {
                        if (vt == nodeType[OPS_ACC1(0, -1, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, -1, 0, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, 0, -1)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, -1, 0)];
                        }
                        if (vt == nodeType[OPS_ACC1(-1, -1, 0)]) {
                            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                                f[OPS_ACC_MD3(xiIdx, 0, 0, -1)];
                        }
                    }
                } break;
                default:
                    break;
            }
            rho += f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
        }
        Real ratio = rhoGiven / rho;
        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] *= ratio;
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a extrapol pressure "
                   "boundary point: KerCutCellExtraolPressure1ST3D");
#endif
    }
}

void KerCutCellEQMDiffuseRefl3D(const Real *givenMacroVars, const int *nodeType,
                                const int *geometryProperty, Real *f,
                                const int *componentId) {
    // This kernel is suitable for any single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    const int equilibriumOrder{2};
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0, 0)];
    if (vt == Vertex_EQMDiffuseRefl) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC2(0, 0, 0)];
        Real u = givenMacroVars[1];
        Real v = givenMacroVars[2];
        Real w = givenMacroVars[3];
        // loop to classify types of discrete velocity i.e., incoming, outgoing
        // and parallel
        // for (int compoIdx = 0; compoIdx < NUMCOMPONENTS; compoIdx++) {
        const int compoIdx{*componentId};
        // compoIdx = *componentId;

        int numOutgoing{0};
        int numIncoming{0};
        int numParallel{0};
        int *outgoing = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        int *incoming = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        int *parallel = new int[COMPOINDEX[2 * compoIdx + 1] -
                                COMPOINDEX[2 * compoIdx] + 1];
        Real rhoIncoming{0};
        Real rhoParallel{0};
        Real deltaRho{0};
        for (int xiIdx = COMPOINDEX[2 * compoIdx];
             xiIdx <= COMPOINDEX[2 * compoIdx + 1]; xiIdx++) {
            Real cx{CS * XI[xiIdx * LATTDIM]};
            Real cy{CS * XI[xiIdx * LATTDIM + 1]};
            Real cz{CS * XI[xiIdx * LATTDIM + 2]};
            BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
            switch (bdt) {
                case BndryDv_Incoming: {
                    incoming[numIncoming] = xiIdx;
                    rhoIncoming += f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                    numIncoming++;
                } break;
                case BndryDv_Outgoing: {
                    outgoing[numOutgoing] = xiIdx;
                    deltaRho +=
                        (2 * WEIGHTS[xiIdx]) * (cx * u + cy * v + cz * w);
                    numOutgoing++;
                } break;
                case BndryDv_Parallel: {
                    parallel[numParallel] = xiIdx;
                    rhoParallel +=
                        CalcBGKFeq(xiIdx, 1, u, v, w, 1, equilibriumOrder);
                    numParallel++;
                } break;
                default:
                    break;
            }
        }
        Real rhoWall = 2 * rhoIncoming / (1 - deltaRho - rhoParallel);
        for (int idx = 0; idx < numParallel; idx++) {
            f[OPS_ACC_MD3(parallel[idx], 0, 0, 0)] = CalcBGKFeq(
                parallel[idx], rhoWall, u, v, w, 1, equilibriumOrder);
        }
        for (int idx = 0; idx < numOutgoing; idx++) {
            int xiIdx = outgoing[idx];
            Real cx{CS * XI[xiIdx * LATTDIM]};
            Real cy{CS * XI[xiIdx * LATTDIM + 1]};
            Real cz{CS * XI[xiIdx * LATTDIM + 2]};
            f[OPS_ACC_MD3(xiIdx, 0, 0, 0)] =
                f[OPS_ACC_MD3(OPP[xiIdx], 0, 0, 0)] +
                2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v + cz * w);
        }
        delete[] outgoing;
        delete[] incoming;
        delete[] parallel;
        //}
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a equilibrium diffuse reflection "
                   "boundary condition point: KerCutCellEQMDiffuseRefl3D");
#endif
    }
}

void KerCutCellNoslipEQN3D(const Real *givenMacroVars, const int *nodeType,
                           Real *f, const int *compoId) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0, 0)];

    if (vt == Vertex_NoslipEQN) {
        Real u = givenMacroVars[1];
        Real v = givenMacroVars[2];
        Real w = givenMacroVars[3];
        Real rhoIntermidate{0};
        Real uIntermidate{0};
        Real vIntermidate{0};
        Real wIntermidate{0};

        for (int xiIdx = COMPOINDEX[2 * (*compoId)];
             xiIdx <= COMPOINDEX[2 * (*compoId) + 1]; xiIdx++) {
            Real cx{CS * XI[xiIdx * LATTDIM]};
            Real cy{CS * XI[xiIdx * LATTDIM + 1]};
            Real cz{CS * XI[xiIdx * LATTDIM + 2]};
            rhoIntermidate += f[OPS_ACC_MD2(xiIdx, 0, 0, 0)];
            uIntermidate += (cx * f[OPS_ACC_MD2(xiIdx, 0, 0, 0)]);
            vIntermidate += (cy * f[OPS_ACC_MD2(xiIdx, 0, 0, 0)]);
            wIntermidate += (cz * f[OPS_ACC_MD2(xiIdx, 0, 0, 0)]);
        }

        for (int xiIdx = COMPOINDEX[2 * (*compoId)];
             xiIdx <= COMPOINDEX[2 * (*compoId) + 1]; xiIdx++) {
            f[OPS_ACC_MD2(xiIdx, 0, 0, 0)] =
                
                       
                                  
                CalcBGKFeq(xiIdx, rhoIntermidate, u, v, w, 1, 2);
        }

    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a equilibrium diffuse reflection "
                   "boundary condition point: KerCutCellEQMDiffuseRefl3D");
#endif
    }
}

void KerCutCellPeriodic3D(const int *nodeType, const int *geometryProperty,
                          Real *f) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0, 0)];
    if (vt == Vertex_Periodic) {
        VertexGeometryTypes vg =
            (VertexGeometryTypes)geometryProperty[OPS_ACC1(0, 0, 0)];
        switch (vg) {
            case VG_IP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                }
                break;
            case VG_IM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                }
                break;
            case VG_JP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                }
                break;
            case VG_JM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                }
                break;
            case VG_KP:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                }
                break;
            case VG_KM:
                for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                    f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                }
                break;
            // as we are dealing with domain boundary, there can be only inner
            // corner
            case VG_IPJP_I: {                             // corner point
                if (vt == nodeType[OPS_ACC0(0, 1, 0)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(1, 0, 0)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
            } break;
            case VG_IPJM_I: {                              // inner corner point
                if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
            } break;
            case VG_IMJP_I: {                             // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
            } break;
            case VG_IMJM_I: {
                // inner corner point
                if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
            } break;

            case VG_IPKP_I: {                             // corner point
                if (vt == nodeType[OPS_ACC0(0, 0, 1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(1, 0, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_IPKM_I: {                              // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_IMKP_I: {                             // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_IMKM_I: {
                // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_JPKP_I: {                             // corner point
                if (vt == nodeType[OPS_ACC0(0, 0, 1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC0(0, 1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_JPKM_I: {                              // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }

                if (vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_JMKP_I: {                             // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_JMKM_I: {
                // inner corner point
                if (vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_IPJPKP_I: {
                if (vt == nodeType[OPS_ACC1(0, 1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_IPJPKM_I: {
                if (vt == nodeType[OPS_ACC1(0, 1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_IPJMKP_I: {
                if (vt == nodeType[OPS_ACC1(0, -1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_IPJMKM_I: {
                if (vt == nodeType[OPS_ACC1(0, -1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, -1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_IMJPKP_I: {
                if (vt == nodeType[OPS_ACC1(0, 1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }
            } break;
            case VG_IMJPKM_I: {
                if (vt == nodeType[OPS_ACC1(0, 1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, -1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            case VG_IMJMKP_I: {
                if (vt == nodeType[OPS_ACC1(0, -1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, 1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KP
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, -1)];
                    }
                }

            } break;
            case VG_IMJMKM_I: {
                if (vt == nodeType[OPS_ACC1(0, -1, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_IM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 1, 0, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, 0, -1)]) {  // VG_JM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 1, 0)];
                    }
                }
                if (vt == nodeType[OPS_ACC1(-1, 0, 0)] &&
                    vt == nodeType[OPS_ACC1(0, -1, 0)]) {  // VG_KM
                    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
                        f[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                            f[OPS_ACC_MD2(xiIndex, 0, 0, 1)];
                    }
                }
            } break;
            default:
                break;
        }
    } else {
#ifdef debug
        ops_printf("%s\n",
                   "Warning: this node is not a periodic boundary "
                   "point: KerCutCellPeriodic");
#endif
    }
}
#endif  // OPS_3D
#endif  /* BOUNDARY_KERNEL_H */
