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

/*! @brief  Kernel functions for boundary conditions
 * @author  Jianping Meng
 * @details Defining kernel functions for various boundary conditions.
 */

// As we are using update-halo method for the discretisation,
// we need to deal with halo points when treating boundary
#ifndef BOUNDARY_KERNEL_INC
#define BOUNDARY_KERNEL_INC
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "model.h"
#include "model_host_device.h"
#include "boundary.h"
#include "flowfield_host_device.h"
#include "boundary_host_device.h"

//TODO: It looks that g_GeometryProperty is not necessary in domain
// (block) boundary conditions

#ifdef OPS_2D// Boundary conditions for two-dimensional problems
// Need to be modified for ImmersedSolid
// void KerCutCellEmbeddedBoundary(const ACC<int> &nodeType,
//                                 const ACC<int> &geometryProperty,
//                                 ACC<Real> &f) {
// #ifdef OPS_2D
//     /*!
//      For the bounce back scheme,We consider zero velocity boundary first.
//      To make sure the velocity at boundary is zero, the implementation
//      is lattice specific.
//      */
//     VertexType vt = (VertexType)nodeType(0, 0);
//     VertexGeometryTypes vg = (VertexGeometryTypes)geometryProperty(0, 0);
//     // TODO to be changed for embeded boundary
//     if (vt == VertexType::ImmersedBoundary) {
//         switch (vt) {
//             case Vertex_EQMDiffuseRefl: {
//                 Real u{0};
//                 Real v{0};
//                 const Real sqrt3 = sqrt(3);
//                 switch (vg) {
//                     case VG_IP: {
//                         const Real f3 = f(3, 0, 0);
//                         const Real f7 = f(7, 0, 0);
//                         const Real f6 = f(6, 0, 0);
//                         const Real rhow =
//                             6 * (f3 + f6 + f7) / (u * u - sqrt3 * u + 1);
//                         f(5, 0, 0) = f7 + rhow * (u + v) / (6 * sqrt3);
//                         f(1, 0, 0) = f3 + 2 * rhow * u / (3 * sqrt3);
//                         f(8, 0, 0) = f6 + rhow * (u - v) / (6 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(2, 0, 0) =
//                             -(u * u - 2 * (1 + sqrt3 * v + v * v) * rhow) / 18;
//                         f(4, 0, 0) =
//                             -((-2 + u * u + 2 * sqrt3 * v - 2 * v * v) * rhow) /
//                             18;
//                     } break;
//                     case VG_IM: {
//                         const Real f5 = f(5, 0, 0);
//                         const Real f1 = f(1, 0, 0);
//                         const Real f8 = f(8, 0, 0);
//                         const Real rhow =
//                             6 * (f1 + f5 + f8) / (u * u + sqrt3 * u + 1);
//                         f(7, 0, 0) = f5 - rhow * (u + v) / (6 * sqrt3);
//                         f(3, 0, 0) = f1 - 2 * rhow * u / (3 * sqrt3);
//                         f(6, 0, 0) = f8 + rhow * (v - u) / (6 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(2, 0, 0) =
//                             -(u * u - 2 * (1 + sqrt3 * v + v * v) * rhow) / 18;
//                         f(4, 0, 0) =
//                             -((-2 + u * u + 2 * sqrt3 * v - 2 * v * v) * rhow) /
//                             18;
//                     } break;
//                     case VG_JP: {
//                         const Real f4 = f(4, 0, 0);
//                         const Real f8 = f(8, 0, 0);
//                         const Real f7 = f(7, 0, 0);
//                         const Real rhow =
//                             6 * (f4 + f8 + f7) / (v * v - sqrt3 * v + 1);
//                         f(2, 0, 0) = f4 + 2 * rhow * v / (3 * sqrt3);
//                         f(6, 0, 0) = f8 + rhow * (v - u) / (6 * sqrt3);
//                         f(5, 0, 0) = f7 + rhow * (u + v) / (6 * sqrt3);
//                         f(1, 0, 0) =
//                             ((2 + 2 * sqrt3 * u + 2 * u * u - v * v) * rhow) /
//                             18;
//                         f(3, 0, 0) =
//                             -((-2 + 2 * sqrt3 * u - 2 * u * u + v * v) * rhow) /
//                             18;
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                     } break;
//                     case VG_JM: {
//                         const Real f2 = f(2, 0, 0);
//                         const Real f5 = f(5, 0, 0);
//                         const Real f6 = f(6, 0, 0);
//                         const Real rhow =
//                             6 * (f2 + f5 + f6) / (v * v + sqrt3 * v + 1);
//                         f(4, 0, 0) = f2 - 2 * rhow * v / (3 * sqrt3);
//                         f(8, 0, 0) = f6 + rhow * (u - v) / (6 * sqrt3);
//                         f(7, 0, 0) = f5 - rhow * (u + v) / (6 * sqrt3);
//                         f(1, 0, 0) =
//                             ((2 + 2 * sqrt3 * u + 2 * u * u - v * v) * rhow) /
//                             18;
//                         f(3, 0, 0) =
//                             -((-2 + 2 * sqrt3 * u - 2 * u * u + v * v) * rhow) /
//                             18;
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                     } break;
//                     case VG_IPJP_I: {
//                         const Real f3 = f(3, 0, 0);
//                         const Real f7 = f(7, 0, 0);
//                         const Real f4 = f(4, 0, 0);
//                         const Real rhow =
//                             (-36 * (f3 + f4 + f7)) /
//                             (-9 + 5 * sqrt3 * u - 3 * u * u + 5 * sqrt3 * v -
//                              3 * u * v - 3 * v * v);
//                         f(1, 0, 0) = f3 + 2 * rhow * u / (3 * sqrt3);
//                         f(5, 0, 0) = f7 + rhow * (u + v) / (6 * sqrt3);
//                         f(2, 0, 0) = f4 + 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(6, 0, 0) = (1 + u * u + sqrt3 * v + v * v -
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(8, 0, 0) = (1 + u * u + u * (sqrt3 - 3 * v) -
//                                       sqrt3 * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IPJM_I: {
//                         const Real f2 = f(2, 0, 0);
//                         const Real f3 = f(3, 0, 0);
//                         const Real f6 = f(6, 0, 0);
//                         const Real rhow =
//                             (-36 * (f2 + f3 + f6)) /
//                             (-9 + 5 * sqrt3 * u - 3 * u * u - 5 * sqrt3 * v +
//                              3 * u * v - 3 * v * v);
//                         f(8, 0, 0) = f6 + rhow * (u - v) / (6 * sqrt3);
//                         f(1, 0, 0) = f3 + 2 * rhow * u / (3 * sqrt3);
//                         f(4, 0, 0) = f2 - 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(5, 0, 0) = (1 + u * u + sqrt3 * v + v * v +
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(7, 0, 0) = (1 - sqrt3 * u + u * u - sqrt3 * v +
//                                       3 * u * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IMJP_I: {
//                         const Real f1 = f(1, 0, 0);
//                         const Real f4 = f(4, 0, 0);
//                         const Real f8 = f(8, 0, 0);
//                         const Real rhow =
//                             (36 * (f1 + f4 + f8)) /
//                             (9 + 5 * sqrt3 * u + 3 * u * u - 5 * sqrt3 * v -
//                              3 * u * v + 3 * v * v);
//                         f(6, 0, 0) = f8 + rhow * (v - u) / (6 * sqrt3);
//                         f(3, 0, 0) = f1 - 2 * rhow * u / (3 * sqrt3);
//                         f(2, 0, 0) = f4 + 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(5, 0, 0) = (1 + u * u + sqrt3 * v + v * v +
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(7, 0, 0) = (1 - sqrt3 * u + u * u - sqrt3 * v +
//                                       3 * u * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IMJM_I: {
//                         const Real f1 = f(1, 0, 0);
//                         const Real f5 = f(5, 0, 0);
//                         const Real f2 = f(2, 0, 0);
//                         const Real rhow =
//                             (36 * (f1 + f2 + f5)) /
//                             (9 + 5 * sqrt3 * u + 3 * u * u + 5 * sqrt3 * v +
//                              3 * u * v + 3 * v * v);
//                         f(3, 0, 0) = f1 - 2 * rhow * u / (3 * sqrt3);
//                         f(4, 0, 0) = f2 - 2 * rhow * v / (3 * sqrt3);
//                         f(7, 0, 0) = f5 - rhow * (u + v) / (6 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(6, 0, 0) = (1 + u * u + sqrt3 * v + v * v -
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(8, 0, 0) = (1 + u * u + u * (sqrt3 - 3 * v) -
//                                       sqrt3 * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IPJP_O: {  // outter corner point
//                         const Real f3 = f(3, 0, 0);
//                         const Real f7 = f(7, 0, 0);
//                         const Real f4 = f(4, 0, 0);
//                         const Real rhow =
//                             (-36 * (f3 + f4 + f7)) /
//                             (-9 + 5 * sqrt3 * u - 3 * u * u + 5 * sqrt3 * v -
//                              3 * u * v - 3 * v * v);
//                         f(1, 0, 0) = f3 + 2 * rhow * u / (3 * sqrt3);
//                         f(5, 0, 0) = f7 + rhow * (u + v) / (6 * sqrt3);
//                         f(2, 0, 0) = f4 + 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(6, 0, 0) = (1 + u * u + sqrt3 * v + v * v -
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(8, 0, 0) = (1 + u * u + u * (sqrt3 - 3 * v) -
//                                       sqrt3 * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IPJM_O: {  // outter corner point
//                         const Real f2 = f(2, 0, 0);
//                         const Real f3 = f(3, 0, 0);
//                         const Real f6 = f(6, 0, 0);
//                         const Real rhow =
//                             (-36 * (f2 + f3 + f6)) /
//                             (-9 + 5 * sqrt3 * u - 3 * u * u - 5 * sqrt3 * v +
//                              3 * u * v - 3 * v * v);
//                         f(8, 0, 0) = f6 + rhow * (u - v) / (6 * sqrt3);
//                         f(1, 0, 0) = f3 + 2 * rhow * u / (3 * sqrt3);
//                         f(4, 0, 0) = f2 - 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(5, 0, 0) = (1 + u * u + sqrt3 * v + v * v +
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(7, 0, 0) = (1 - sqrt3 * u + u * u - sqrt3 * v +
//                                       3 * u * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IMJP_O: {  // outter corner point
//                         const Real f1 = f(1, 0, 0);
//                         const Real f4 = f(4, 0, 0);
//                         const Real f8 = f(8, 0, 0);
//                         const Real rhow =
//                             (36 * (f1 + f4 + f8)) /
//                             (9 + 5 * sqrt3 * u + 3 * u * u - 5 * sqrt3 * v -
//                              3 * u * v + 3 * v * v);
//                         f(6, 0, 0) = f8 + rhow * (v - u) / (6 * sqrt3);
//                         f(3, 0, 0) = f1 - 2 * rhow * u / (3 * sqrt3);
//                         f(2, 0, 0) = f4 + 2 * rhow * v / (3 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(5, 0, 0) = (1 + u * u + sqrt3 * v + v * v +
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(7, 0, 0) = (1 - sqrt3 * u + u * u - sqrt3 * v +
//                                       3 * u * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     case VG_IMJM_O: {  // outter corner point
//                         const Real f1 = f(1, 0, 0);
//                         const Real f5 = f(5, 0, 0);
//                         const Real f2 = f(2, 0, 0);
//                         const Real rhow =
//                             (36 * (f1 + f2 + f5)) /
//                             (9 + 5 * sqrt3 * u + 3 * u * u + 5 * sqrt3 * v +
//                              3 * u * v + 3 * v * v);
//                         f(3, 0, 0) = f1 - 2 * rhow * u / (3 * sqrt3);
//                         f(4, 0, 0) = f2 - 2 * rhow * v / (3 * sqrt3);
//                         f(7, 0, 0) = f5 - rhow * (u + v) / (6 * sqrt3);
//                         f(0, 0, 0) = 2 * rhow * (2 - u * u - v * v) / 9;
//                         f(6, 0, 0) = (1 + u * u + sqrt3 * v + v * v -
//                                       u * (sqrt3 + 3 * v)) *
//                                      rhow / 36;
//                         f(8, 0, 0) = (1 + u * u + u * (sqrt3 - 3 * v) -
//                                       sqrt3 * v + v * v) *
//                                      rhow / 36;
//                     } break;
//                     default:
//                         break;
//                 }  // vg
//             }      // case Vertex_EQMDiffRefl
//             break;
//             case Vertex_KineticDiffuseWall: {
//                 Real u = 0;
//                 Real v = 0;  // means non-moving wall
//                 Real wallNormalVector[]{0, 0};
//                 const Real sqrt2Inverse = 1 / sqrt(2);
//                 switch (vg) {
//                     case VG_IP: {
//                         wallNormalVector[0] = 1;
//                         wallNormalVector[1] = 0;
//                     } break;
//                     case VG_IM: {
//                         wallNormalVector[0] = -1;
//                         wallNormalVector[1] = 0;
//                     } break;
//                     case VG_JP: {
//                         wallNormalVector[0] = 0;
//                         wallNormalVector[1] = 1;
//                     } break;
//                     case VG_JM: {
//                         wallNormalVector[0] = 0;
//                         wallNormalVector[1] = -1;
//                     } break;
//                     case VG_IPJP_I: {
//                         wallNormalVector[0] = sqrt2Inverse;
//                         wallNormalVector[1] = sqrt2Inverse;
//                     } break;
//                     case VG_IPJM_I: {
//                         wallNormalVector[0] = sqrt2Inverse;
//                         wallNormalVector[1] = -sqrt2Inverse;
//                     } break;
//                     case VG_IMJP_I: {
//                         wallNormalVector[0] = -sqrt2Inverse;
//                         wallNormalVector[1] = sqrt2Inverse;
//                     } break;
//                     case VG_IMJM_I: {
//                         wallNormalVector[0] = -sqrt2Inverse;
//                         wallNormalVector[1] = -sqrt2Inverse;
//                     } break;
//                     case VG_IPJP_O: {
//                         wallNormalVector[0] = sqrt2Inverse;
//                         wallNormalVector[1] = sqrt2Inverse;
//                     } break;
//                     case VG_IPJM_O: {
//                         wallNormalVector[0] = sqrt2Inverse;
//                         wallNormalVector[1] = -sqrt2Inverse;
//                     } break;
//                     case VG_IMJP_O: {
//                         wallNormalVector[0] = -sqrt2Inverse;
//                         wallNormalVector[1] = sqrt2Inverse;
//                     } break;
//                     case VG_IMJM_O: {
//                         wallNormalVector[0] = -sqrt2Inverse;
//                         wallNormalVector[1] = -sqrt2Inverse;
//                     } break;

//                     default:
//                         break;
//                 }
//                 Real outFlux = 0;  // flow into wall
//                 Real inFlux = 0;   // flow into fluid bulk
//                 for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
//                     const Real cx = XI[xiIndex * LATTDIM];
//                     const Real cy = XI[xiIndex * LATTDIM + 1];
//                     Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
//                                       (CS * cy - v) * wallNormalVector[1];
//                     if (cDotNormal < 0) {
//                         outFlux += (-cDotNormal * f(xiIndex, 0, 0));
//                     }
//                     if (cDotNormal > 0) {
//                         Real cu = (CS * cx * u + CS * cy * v);
//                         inFlux +=
//                             (cDotNormal * WEIGHTS[xiIndex] *
//                              (1 + cu +
//                               0.5 * (cu * cu - (u * u + v * v))));  // i.e., the
//                         // equilibrium
//                     }
//                 }
//                 Real rho = outFlux / inFlux;
//                 for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
//                     const Real cx = XI[xiIndex * LATTDIM];
//                     const Real cy = (int)XI[xiIndex * LATTDIM + 1];
//                     Real cDotNormal = (CS * cx - u) * wallNormalVector[0] +
//                                       (CS * cy - v) * wallNormalVector[1];
//                     if (cDotNormal >= 0) {
//                         Real cu = (CS * cx * u + CS * cy * v);
//                         f(xiIndex, 0, 0) =
//                             (rho * WEIGHTS[xiIndex] *
//                              (1 + cu + 0.5 * (cu * cu - (u * u + v * v))));
//                     }
//                 }
//             }  // case Vertex_KineticDiffuseWall
//             break;
//             default:
// #ifdef debug
//                 ops_printf("%s\n",
//                            "Warning: KerCutCellImmersedBoundary: there seems a "
//                            "boundary condition that has note been implemented");
// #endif
//                 break;
//         }
//     }
// #endif OPS_2D
// }

void KerCutCellExtrapolPressure1ST(ACC<Real> &f, const ACC<int> &nodeType,
                                   const ACC<int> &geometryProperty,
                                   const Real *givenBoundaryVars,
                                   const int *surface, const int *lattIdx) {
#ifdef OPS_2D

    const BoundarySurface boundarySurface{(BoundarySurface)(*surface)};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    Real rhoGiven = givenBoundaryVars[0];
    Real rho = 0;
    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        switch (vg) {
            case VG_IP: {
                if (cx > 0) {
                    f(xiIdx, 0, 0) = f(xiIdx, 1, 0);
                }
            } break;
            case VG_IM: {
                if (cx < 0) {
                    f(xiIdx, 0, 0) = f(xiIdx, -1, 0);
                }
            } break;
            case VG_JP: {
                if (cy > 0) {
                    f(xiIdx, 0, 0) = f(xiIdx, 0, 1);
                }
            } break;
            case VG_JM: {
                if (cy < 0) {
                    f(xiIdx, 0, 0) = f(xiIdx, 0, -1);
                }
            } break;
            case VG_IPJP_I: {
                if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0) = f(xiIdx, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0) = f(xiIdx, 0, 1);
                    }
                }
            } break;
            case VG_IPJM_I: {
                if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0) = f(xiIdx, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0) = f(xiIdx, 0, -1);
                    }
                }
            } break;
            case VG_IMJP_I: {
                if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0) = f(xiIdx, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0) = f(xiIdx, 0, 1);
                    }
                }
            } break;
            case VG_IMJM_I: {
                if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0) = f(xiIdx, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0) = f(xiIdx, 0, -1);
                    }
                }
            } break;
            default:
                break;
        }
        rho += f(xiIdx, 0, 0);
    }
    Real ratio = rhoGiven / rho;
    for (int xiIdx = lattIdx[0]; xiIdx < lattIdx[1]; xiIdx++) {
        f(xiIdx, 0, 0) *= ratio;
    }
#endif  // OPS_2D
}

// void KerCutCellBounceBack(const ACC<int> &nodeType,
//                           const ACC<int> &geometryProperty, ACC<Real> &f) {
// #ifdef OPS_2D

//     VertexGeometryTypes vg = (VertexGeometryTypes)geometryProperty(0, 0);
//     switch (vg) {
//         case VG_IP: {
//             f(5, 0, 0) = f(7, 0, 0);
//             f(1, 0, 0) = f(3, 0, 0);
//             f(8, 0, 0) = f(6, 0, 0);
//         } break;
//         case VG_IM: {
//             f(7, 0, 0) = f(5, 0, 0);
//             f(3, 0, 0) = f(1, 0, 0);
//             f(6, 0, 0) = f(8, 0, 0);
//         } break;
//         case VG_JP: {
//             f(2, 0, 0) = f(4, 0, 0);
//             f(6, 0, 0) = f(8, 0, 0);
//             f(5, 0, 0) = f(7, 0, 0);
//         } break;
//         case VG_JM: {
//             f(4, 0, 0) = f(2, 0, 0);
//             f(8, 0, 0) = f(6, 0, 0);
//             f(7, 0, 0) = f(5, 0, 0);
//         } break;
//         case VG_IPJP_I:  // inner corner point
//             f(5, 0, 0) = f(7, 0, 0);
//             break;
//         case VG_IPJM_I:  // inner corner point
//             f(8, 0, 0) = f(6, 0, 0);
//             break;
//         case VG_IMJP_I:  // inner corner point
//             f(6, 0, 0) = f(8, 0, 0);
//             break;
//         case VG_IMJM_I:  // inner corner point
//             f(7, 0, 0) = f(5, 0, 0);
//             break;
//         default:
//             break;
//     }
// #endif  // OPS_2D
// }

void KerCutCellEQMDiffuseRefl(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars,
                                const int *lattIdx) {
#ifdef OPS_2D
    // This kernel is suitable for a single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    const int equilibriumOrder{2};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
#ifdef CPU
#if DebugLevel >= 2
    ops_printf(
        "KerCutCellEQMDiffuseRefl: We received the following "
        "conditions for the surface %i:\n",
        geometryProperty(0, 0));
    ops_printf("U=%f, V=%f\n", u, v);
#endif
#endif
    int numOutgoing{0};
    int numIncoming{0};
    int numParallel{0};
    int lattStart{lattIdx[0]};
    int lattEnd{lattIdx[1]};
    int *outgoing = new int[lattEnd - lattStart + 1];
    int *incoming = new int[lattEnd - lattStart + 1];
    int *parallel = new int[lattEnd - lattStart + 1];
    Real rhoIncoming{0};
    Real rhoParallel{0};
    Real deltaRho{0};
    for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        BndryDvType bdt = FindBdyDvType(vg, &XI[xiIdx * LATTDIM]);
        switch (bdt) {
            case BndryDv_Incoming: {
                incoming[numIncoming] = xiIdx;
                rhoIncoming += f(xiIdx, 0, 0);
                numIncoming++;
            } break;
            case BndryDv_Outgoing: {
                outgoing[numOutgoing] = xiIdx;
                deltaRho += (2 * WEIGHTS[xiIdx]) * (cx * u + cy * v);
                numOutgoing++;
            } break;
            case BndryDv_Parallel: {
                parallel[numParallel] = xiIdx;
                rhoParallel += CalcBGKFeq(xiIdx, 1, u, v, 1, equilibriumOrder);
                numParallel++;
            } break;
            default:
                break;
        }
    }
    Real rhoWall = 2 * rhoIncoming / (1 - deltaRho - rhoParallel);

#ifdef CPU
#if DebugLevel >= 2
    ops_printf("Calculated wall density =  %f\n", rhoWall);
#endif
#endif
    for (int idx = 0; idx < numParallel; idx++) {
        f(parallel[idx], 0, 0) =
            CalcBGKFeq(parallel[idx], rhoWall, u, v, 1, equilibriumOrder);
    }
    for (int idx = 0; idx < numOutgoing; idx++) {
        int xiIdx = outgoing[idx];
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        f(xiIdx, 0, 0) = f(OPP[xiIdx], 0, 0) +
                         2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
    }
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;

#endif  // OPS_2D
}

void KerCutCellPeriodic(ACC<Real> &f, const ACC<int> &nodeType,
                        const ACC<int> &geometryProperty, const int *lattIdx,
                        const int *surface) {
#ifdef OPS_2D
    const int xiStartPos{lattIdx[0]};
    const int xiEndPos{lattIdx[1]};
    const BoundarySurface boundarySurface{(BoundarySurface)(*surface)};

    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    switch (vg) {
        case VG_IP:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0) = f(xiIndex, -1, 0);
            }
            break;
        case VG_IM:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0) = f(xiIndex, 1, 0);
            }
            break;
        case VG_JP:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0) = f(xiIndex, 0, -1);
            }
            break;
        case VG_JM:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0) = f(xiIndex, 0, 1);
            }
            break;
            // There are only inner corners for block boundaries
        case VG_IPJP_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, -1, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 0, -1);
                }
            }
        } break;
        case VG_IPJM_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, -1, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 0, 1);
                }
            }
        } break;
        case VG_IMJP_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 1, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 0, -1);
                }
            }
        } break;
        case VG_IMJM_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 1, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0) = f(xiIndex, 0, 1);
                }
            }
        } break;
        default:
            break;
    }
#ifdef CPU
    for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
        const Real res{f(xiIndex, 0, 0)};
        if (isnan(res) || res <= 0 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n at the surface %i\n",
                res, xiIndex, geometryProperty(0, 0));
            assert(!(isnan(res) || res <= 0 || isinf(res)));
        }
    }
#endif

#endif  // OPS_2D
}

void KerCutCellZouHeVelocity(const Real *givenMacroVars,
                             const ACC<int> &nodeType,
                             const ACC<int> &geometryProperty,
                             const ACC<Real> &macroVars, ACC<Real> &f) {
#ifdef OPS_2D
    /*!
    Note: This boundary condition requires both stream and collision happenning
    at a boundary point.
    Note: This boundary condition is lattice specific.
    */

    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    Real rho{0};
    Real u{givenMacroVars[1]};
    Real v{givenMacroVars[2]};
    Real sqrt3 = sqrt(3);
    switch (vg) {
        case VG_IP: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f3 = f(3, 0, 0);
            Real f4 = f(4, 0, 0);
            Real f6 = f(6, 0, 0);
            Real f7 = f(7, 0, 0);
            rho =
                sqrt3 * (f0 + f2 + 2 * f3 + f4 + 2 * f6 + 2 * f7) / (sqrt3 - u);
            f(1, 0, 0) = (2 * sqrt3 * rho * u + 9 * f3) / 9.0;
            f(5, 0, 0) = (sqrt3 * rho * u + 3 * sqrt3 * rho * v - 9 * f2 +
                          9 * f4 + 18 * f7) /
                         18.0;
            f(8, 0, 0) = (sqrt3 * rho * u - 3 * sqrt3 * rho * v + 9 * f2 -
                          9 * f4 + 18 * f6) /
                         18.0;
        } break;
        case VG_IM: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f4 = f(4, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f5 = f(5, 0, 0);
            Real f8 = f(8, 0, 0);
            rho = (sqrt3 * f0 + 2 * sqrt3 * f1 + sqrt3 * f2 + sqrt3 * f4 +
                   2 * sqrt3 * f5 + 2 * sqrt3 * f8) /
                  (sqrt3 + u);
            f(3, 0, 0) = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
            f(6, 0, 0) = (-(sqrt3 * u * rho) + 3 * sqrt3 * v * rho - 9 * f2 +
                          9 * f4 + 18 * f8) /
                         18.0;
            f(7, 0, 0) = (-(sqrt3 * u * rho) - 3 * sqrt3 * v * rho + 9 * f2 -
                          9 * f4 + 18 * f5) /
                         18.0;
        } break;
        case VG_JP: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f3 = f(3, 0, 0);
            Real f4 = f(4, 0, 0);
            Real f7 = f(7, 0, 0);
            Real f8 = f(8, 0, 0);
            rho = (sqrt3 * f0 + sqrt3 * f1 + sqrt3 * f3 + 2 * sqrt3 * f4 +
                   2 * sqrt3 * f7 + 2 * sqrt3 * f8) /
                  (sqrt3 - v);
            f(2, 0, 0) = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
            f(5, 0, 0) = (3 * sqrt3 * u * rho + sqrt3 * v * rho - 9 * f1 +
                          9 * f3 + 18 * f7) /
                         18.0;
            f(6, 0, 0) = (-3 * sqrt3 * u * rho + sqrt3 * v * rho + 9 * f1 -
                          9 * f3 + 18 * f8) /
                         18.0;
        } break;
        case VG_JM: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f3 = f(3, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f5 = f(5, 0, 0);
            Real f6 = f(6, 0, 0);
            rho = (sqrt3 * f0 + sqrt3 * f1 + 2 * sqrt3 * f2 + sqrt3 * f3 +
                   2 * sqrt3 * f5 + 2 * sqrt3 * f6) /
                  (sqrt3 + v);
            f(4, 0, 0) = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
            f(7, 0, 0) = (-3 * sqrt3 * u * rho - sqrt3 * v * rho + 9 * f1 -
                          9 * f3 + 18 * f5) /
                         18.0;
            f(8, 0, 0) = (3 * sqrt3 * u * rho - sqrt3 * v * rho - 9 * f1 +
                          9 * f3 + 18 * f6) /
                         18.0;
        } break;
        case VG_IPJM_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f6 = f(6, 0, 0);
            Real f3 = f(3, 0, 0);
            rho = macroVars(0, 1, -1);
            f(1, 0, 0) = (2 * sqrt3 * u * rho + 9 * f3) / 9.0;
            f(5, 0, 0) = (9 * rho - 2 * sqrt3 * u * rho + 3 * sqrt3 * v * rho -
                          9 * f0 - 18 * f2 - 18 * f3 - 18 * f6) /
                         18.0;
            f(7, 0, 0) = (9 * rho - 3 * sqrt3 * u * rho + 2 * sqrt3 * v * rho -
                          9 * f0 - 18 * f2 - 18 * f3 - 18 * f6) /
                         18.0;
            f(4, 0, 0) = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
            f(8, 0, 0) = (sqrt3 * u * rho - sqrt3 * v * rho + 18 * f6) / 18.0;
        } break;
        case VG_IPJP_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f4 = f(4, 0, 0);
            Real f3 = f(3, 0, 0);
            Real f7 = f(7, 0, 0);
            rho = macroVars(0, 1, 1);
            f(1, 0, 0) = (2 * sqrt3 * u * rho + 9 * f3) / 9.0;
            f(5, 0, 0) = (sqrt3 * u * rho + sqrt3 * v * rho + 18 * f7) / 18.0;
            f(8, 0, 0) = (9 * rho - 2 * sqrt3 * u * rho - 3 * sqrt3 * v * rho -
                          9 * f0 - 18 * f3 - 18 * f4 - 18 * f7) /
                         18.0;
            f(2, 0, 0) = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
            f(6, 0, 0) = (9 * rho - 3 * sqrt3 * u * rho - 2 * sqrt3 * v * rho -
                          9 * f0 - 18 * f3 - 18 * f4 - 18 * f7) /
                         18.0;
        } break;
        case VG_IMJP_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f8 = f(8, 0, 0);
            Real f4 = f(4, 0, 0);
            rho = macroVars(0, -1, 1);
            f(5, 0, 0) = (9 * rho + 3 * sqrt3 * u * rho - 2 * sqrt3 * v * rho -
                          9 * f0 - 18 * f1 - 18 * f4 - 18 * f8) /
                         18.0;
            f(2, 0, 0) = (2 * sqrt3 * v * rho + 9 * f4) / 9.0;
            f(6, 0, 0) =
                (-(sqrt3 * u * rho) + sqrt3 * v * rho + 18 * f8) / 18.0;
            f(3, 0, 0) = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
            f(7, 0, 0) = (9 * rho + 2 * sqrt3 * u * rho - 3 * sqrt3 * v * rho -
                          9 * f0 - 18 * f1 - 18 * f4 - 18 * f8) /
                         18.0;
        } break;
        case VG_IMJM_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f5 = f(5, 0, 0);
            rho = macroVars(0, -1, -1);
            f(6, 0, 0) = (9 * rho + 2 * sqrt3 * u * rho + 3 * sqrt3 * v * rho -
                          9 * f0 - 18 * f1 - 18 * f2 - 18 * f5) /
                         18.0;
            f(3, 0, 0) = (-2 * sqrt3 * u * rho + 9 * f1) / 9.0;
            f(7, 0, 0) =
                (-(sqrt3 * u * rho) - sqrt3 * v * rho + 18 * f5) / 18.0;
            f(4, 0, 0) = (-2 * sqrt3 * v * rho + 9 * f2) / 9.0;
            f(8, 0, 0) = (9 * rho + 3 * sqrt3 * u * rho + 2 * sqrt3 * v * rho -
                          9 * f0 - 18 * f1 - 18 * f2 - 18 * f5) /
                         18.0;
        } break;
        default:
            break;
    }

#endif  // OPS_2D
}

#endif // OPS_2D outer

// Boundary conditions for three-dimensional problems
#ifdef OPS_3D
void KerCutCellExtrapolPressure1ST3D(ACC<Real> &f, const ACC<int> &nodeType,
                                     const ACC<int> &geometryProperty,
                                     const Real *givenBoundaryVars,
                                     const int *surface,
                                     const int *lattIdx) {
#ifdef OPS_3D

    const BoundarySurface boundarySurface{(BoundarySurface)(*surface)};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    Real rhoGiven = givenBoundaryVars[0];
    Real rho = 0;
    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        Real cz{CS * XI[xiIdx * LATTDIM + 2]};
        switch (vg) {
            case VG_IP: {
                if (cx > 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                }
            } break;
            case VG_IM: {
                if (cx < 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                }
            } break;
            case VG_JP: {
                if (cy > 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                }
            } break;
            case VG_JM: {
                if (cy < 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                }
            } break;
            case VG_KP: {
                if (cz > 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                }
            } break;
            case VG_KM: {
                if (cz < 0) {
                    f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                }
            } break;
            case VG_IPJP_I: {
                if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                }
            } break;
            case VG_IPJM_I: {
                if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                }
            } break;
            case VG_IMJP_I: {
                if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                }
            } break;
            case VG_IMJM_I: {
                if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                }
            } break;
            case VG_IPKP_I: {
                if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IPKM_I: {
                if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_IMKP_I: {
                if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IMKM_I: {
                if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_JPKP_I: {
                if ((cy >= 0 && cz > 0) || (cy > 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_JPKM_I: {
                if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_JMKP_I: {
                if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_JMKM_I: {
                if ((cy <= 0 && cz < 0) || (cy < 0 && cz == 0)) {
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_IPJPKP_I: {
                if ((cx >= 0 && cy >= 0 && cz >= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IPJPKM_I: {
                if ((cx >= 0 && cy >= 0 && cz <= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_IPJMKP_I: {
                if ((cx >= 0 && cy <= 0 && cz >= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IPJMKM_I: {
                if ((cx >= 0 && cy <= 0 && cz <= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Left) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_IMJPKP_I: {
                if ((cx <= 0 && cy >= 0 && cz >= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IMJPKM_I: {
                if ((cx <= 0 && cy >= 0 && cz <= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Bottom) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            case VG_IMJMKP_I: {
                if ((cx <= 0 && cy <= 0 && cz >= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Back) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, 1);
                    }
                }
            } break;
            case VG_IMJMKM_I: {
                if ((cx <= 0 && cy <= 0 && cz <= 0) &&
                    (cx != 0 || cy != 0 || cz != 0)) {
                    if (boundarySurface == BoundarySurface::Right) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, -1, 0, 0);
                    }
                    if (boundarySurface == BoundarySurface::Top) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, -1, 0);
                    }
                    if (boundarySurface == BoundarySurface::Front) {
                        f(xiIdx, 0, 0, 0) = f(xiIdx, 0, 0, -1);
                    }
                }
            } break;
            default:
                break;
        }
        rho += f(xiIdx, 0, 0, 0);
    }
    Real ratio = rhoGiven / rho;
    for (int xiIdx = lattIdx[0]; xiIdx < lattIdx[1]; xiIdx++) {
        f(xiIdx, 0, 0, 0) *= ratio;
    }
#endif  // OPS_3D
}

void KerCutCellEQMDiffuseRefl3D(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars,
                                const int *lattIdx) {
#ifdef OPS_3D
    // This kernel is suitable for any single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    const int equilibriumOrder{2};
    const VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
    Real w = givenMacroVars[2];
#ifdef CPU
#if DebugLevel >= 2
    ops_printf(
        "KerCutCellEQMDiffuseRefl3D: We received the following "
        "conditions for the surface %i:\n",
        geometryProperty(0, 0, 0));
    ops_printf("U=%f, V=%f, W=%f\n", u, v, w);
#endif
#endif
    int numOutgoing{0};
    int numIncoming{0};
    int numParallel{0};
    int lattStart{lattIdx[0]};
    int lattEnd{lattIdx[1]};
    int *outgoing = new int[lattEnd - lattStart + 1];
    int *incoming = new int[lattEnd - lattStart + 1];
    int *parallel = new int[lattEnd - lattStart + 1];
    Real rhoIncoming{0};
    Real rhoParallel{0};
    Real deltaRho{0};
    for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        Real cz{CS * XI[xiIdx * LATTDIM + 2]};
        BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
        switch (bdt) {
            case BndryDv_Incoming: {
                incoming[numIncoming] = xiIdx;
                rhoIncoming += f(xiIdx, 0, 0, 0);
                numIncoming++;
            } break;
            case BndryDv_Outgoing: {
                outgoing[numOutgoing] = xiIdx;
                deltaRho += (2 * WEIGHTS[xiIdx]) * (cx * u + cy * v + cz * w);
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
#ifdef CPU
#if DebugLevel >= 2
    ops_printf("Calculated wall density =  %f\n", rhoWall);
#endif
#endif
    for (int idx = 0; idx < numParallel; idx++) {
        f(parallel[idx], 0, 0, 0) =
            CalcBGKFeq(parallel[idx], rhoWall, u, v, w, 1, equilibriumOrder);
    }
    for (int idx = 0; idx < numOutgoing; idx++) {
        int xiIdx = outgoing[idx];
        Real cx{CS * XI[xiIdx * LATTDIM]};
        Real cy{CS * XI[xiIdx * LATTDIM + 1]};
        Real cz{CS * XI[xiIdx * LATTDIM + 2]};
        f(xiIdx, 0, 0, 0) =
            f(OPP[xiIdx], 0, 0, 0) +
            2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v + cz * w);
#ifdef CPU
        const Real res{f(xiIdx, 0, 0, 0)};
        if (isnan(res) || res <= 0 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n",
                res, xiIdx);
            assert(!(isnan(res) || res <= 0 || isinf(res)));
        }
#endif
    }
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
#endif //OPS_3D
}

void KerCutCellPeriodic3D(ACC<Real> &f, const ACC<int> &nodeType,
                          const ACC<int> &geometryProperty,
                          const int *lattIdx, const int* surface) {
#ifdef OPS_3D
    const int xiStartPos{lattIdx[0]};
    const int xiEndPos{lattIdx[1]};
    const BoundarySurface boundarySurface{(BoundarySurface)(*surface)};

    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    switch (vg) {
        case VG_IP:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
            }
            break;
        case VG_IM:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
            }
            break;
        case VG_JP:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
            }
            break;
        case VG_JM:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
            }
            break;
        case VG_KP:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
            }
            break;
        case VG_KM:
            for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
            }
            break;
            // There are only inner corners for block boundaries
        case VG_IPJP_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
        } break;
        case VG_IPJM_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
        } break;
        case VG_IMJP_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
        } break;
        case VG_IMJM_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
        } break;

        case VG_IPKP_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_IPKM_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_IMKP_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_IMKM_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_JPKP_I: {
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_JPKM_I: {
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_JMKP_I: {
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_JMKM_I: {
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_IPJPKP_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_IPJPKM_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_IPJMKP_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Left) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_IPJMKM_I: {
            // VG_IP
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, -1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_IMJPKP_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }
        } break;
        case VG_IMJPKM_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JP
            if (boundarySurface == BoundarySurface::Bottom) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, -1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        case VG_IMJMKP_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KP
            if (boundarySurface == BoundarySurface::Back) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, -1);
                }
            }

        } break;
        case VG_IMJMKM_I: {
            // VG_IM
            if (boundarySurface == BoundarySurface::Right) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 1, 0, 0);
                }
            }
            // VG_JM
            if (boundarySurface == BoundarySurface::Top) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 1, 0);
                }
            }
            // VG_KM
            if (boundarySurface == BoundarySurface::Front) {
                for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
                    f(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 1);
                }
            }
        } break;
        default:
            break;
    }
#ifdef CPU
    for (int xiIndex = xiStartPos; xiIndex <= xiEndPos; xiIndex++) {
        const Real res{f(xiIndex, 0, 0, 0)};
        if (isnan(res) || res <= 0 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n at the surface %i\n",
                res, xiIndex, geometryProperty(0, 0, 0));
            assert(!(isnan(res) || res <= 0 || isinf(res)));
        }
    }
#endif

#endif  // OPS_3D
}
#endif //OPS_3D
#endif // BOUNDARY_KERNEL_INC