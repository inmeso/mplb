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

void KerCutCellEQMDiffuseReflADF(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars, const ACC<Real> &acceleration,const Real *dt,
                                const int *lattIdx) {
#ifdef OPS_2D
    // This kernel is suitable for a single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Wall) {
    const int equilibriumOrder{2};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
    Real g[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
    Real CS=1/sqrt(3);
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
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
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
                rhoParallel += CalcBGKFeqFE(xiIdx, 1, u, v, g, dt, 1, equilibriumOrder);
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
            CalcBGKFeqFE(parallel[idx], rhoWall, u, v, g, dt, 1, equilibriumOrder);
    }
    for (int idx = 0; idx < numOutgoing; idx++) {
        int xiIdx = outgoing[idx];
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        f(xiIdx, 0, 0) = f(OPP[xiIdx], 0, 0) +
                         2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
    }
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
    }
#endif  // OPS_2D
}

void KerCutCellEQMDiffuseReflADG(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars, const ACC<Real> &Conc, const ACC<Real> &Mu,
                                const Real *dt, const int *lattIdx) {
#ifdef OPS_2D
    // This kernel is suitable for a single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Wall) {
    const int equilibriumOrder{2};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    Real mu{Mu(0,0)};
    Real C{Conc(0,0)};
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
    Real CS=1/sqrt(3);
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
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        BndryDvType bdt = FindBdyDvType(vg, &XI[xiIdx * LATTDIM]);
        switch (bdt) {
            case BndryDv_Incoming: {
                incoming[numIncoming] = xiIdx;
                //rhoIncoming += f(xiIdx, 0, 0);
                rhoIncoming += CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, equilibriumOrder);
                numIncoming++;
            } break;
            case BndryDv_Outgoing: {
                outgoing[numOutgoing] = xiIdx;
                //deltaRho += (2 * WEIGHTS[xiIdx]) * (cx * u + cy * v);
                deltaRho += CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, equilibriumOrder);
                //std::cout<<CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, equilibriumOrder)<<" ";
                numOutgoing++;
            } break;
            case BndryDv_Parallel: {
                parallel[numParallel] = xiIdx;
                //rhoParallel += CalcBGKGeqFE(xiIdx, 1, u, v, mu, dt, equilibriumOrder);
                //rhoParallel += f(xiIdx, 0, 0);
                if (xiIdx>lattIdx[0]){
                    rhoParallel += CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, equilibriumOrder);
                    //std::cout<<rhoParallel<<" ";
                }
                numParallel++;
            } break;
            default:
                break;
        }
    }
    //Real rhoWall = 2 * rhoIncoming / (1 - rhoParallel - deltaRho);
    Real rhoWall = rhoIncoming+deltaRho+rhoParallel;
    //std::cout<<rhoWall<<" ";
    //std::cout<<rhoWall<<"  ";
    //std::cout<<rhoWall<<" "<<rhoIncoming<<" "<<deltaRho<<" "<<rhoParallel<<"   ";
#ifdef CPU
#if DebugLevel >= 2
    ops_printf("Calculated wall density =  %f\n", rhoWall);
#endif
#endif
    for (int idx = 0; idx < numParallel; idx++) {
        //std::cout<<parallel[idx]<<"  ";
        if (idx>lattIdx[0]){
            
            f(parallel[idx], 0, 0) =
                CalcBGKGeqFE(parallel[idx], C, u, v, mu, dt, equilibriumOrder);
        }
        else if (idx==lattIdx[0]){
            f(parallel[idx], 0, 0) =C-rhoWall;
        }
    }
    for (int idx = 0; idx < numOutgoing; idx++) {
                int xiIdx = outgoing[idx];
                Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
                Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
                Real opps={OPP[xiIdx]};
                //std::cout<<f(opps,0,0)<<" ";
                f(xiIdx, 0, 0)=CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, equilibriumOrder);
                //f(xiIdx, 0, 0) = f(opps, 0, 0);
                //                    2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
    }
    /*switch (vg) {
            case VG_IP: {
                for (int idx = 0; idx < numOutgoing; idx++) {
                int xiIdx = outgoing[idx];
                //std::cout<<OPP[xiIdx]<<" ";
                Real opps={8+OPP[xiIdx]};
                Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
                Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
                if (cx>0) {
                    f(xiIdx, 0, 0) = f(opps, 0, 0)+
                                    2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                    }
                    //std::cout<<rhoWall<<"  ";
                    //f(xiIdx, 0, 0) +=
                    //                 2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                }
            } break;
            case VG_IM: {
                for (int idx = 0; idx < numOutgoing; idx++) {
                int xiIdx = outgoing[idx];
                Real opps={8+OPP[xiIdx]};
                Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
                Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
                if (cx<0) {
                    f(xiIdx, 0, 0) = f(opps, 0, 0)+
                                    2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                    }
                    //f(xiIdx, 0, 0) +=
                    //                 2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                }
            } break;
            case VG_JP: {
                for (int idx = 0; idx < numOutgoing; idx++) {
                int xiIdx = outgoing[idx];
                Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
                Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
                if (cy<0) {
                    f(xiIdx, 0, 0) = f(OPP[xiIdx], 0, 0)+
                                    2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                    }
                    //std::cout<<rhoWall<<"  ";
                    //f(xiIdx, 0, 0) +=
                    //                 2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                }
            } break;
            case VG_JM: {
                for (int idx = 0; idx < numOutgoing; idx++) {
                int xiIdx = outgoing[idx];
                Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
                Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
                if (cy>0) {
                    f(xiIdx, 0, 0) = f(OPP[xiIdx], 0, 0)+
                                    2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                    }
                    //f(xiIdx, 0, 0) +=
                    //                 2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v);
                }
            } break;
        }*/
            
    
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
    }
#endif  // OPS_2D
}

void KerCutCellZouHeVelocityADG(const Real *givenMacroVars,
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
    Real css = 1/sqrt(3);
    Real invcss=1/css;
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
                css * (f0 + f2 + 2 * f3 + f4 + 2 * f6 + 2 * f7) / (css - u);
            f(1, 0, 0) = (6 * invcss * rho * u + 9 * f3) / 9.0;
            f(5, 0, 0) = (3*invcss * rho * u + 9 * invcss * rho * v - 9 * f2 +
                          9 * f4 + 18 * f7) /
                         18.0;
            f(8, 0, 0) = (3*invcss * rho * u - 9 * invcss * rho * v + 9 * f2 -
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
            rho = (css * f0 + 2 * css * f1 + css * f2 + css * f4 +
                   2 * css * f5 + 2 * css * f8) /
                  (css + u);
            f(3, 0, 0) = (-6 * invcss * u * rho + 9 * f1) / 9.0;
            f(6, 0, 0) = (-3*(invcss * u * rho) + 9 * invcss * v * rho - 9 * f2 +
                          9 * f4 + 18 * f8) /
                         18.0;
            f(7, 0, 0) = (-3*(invcss * u * rho) - 9 * invcss * v * rho + 9 * f2 -
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
            rho = (css * f0 + css * f1 + css * f3 + 2 * css * f4 +
                   2 * css * f7 + 2 * css * f8) /
                  (css - v);
            f(2, 0, 0) = (6 * invcss * v * rho + 9 * f4) / 9.0;
            f(5, 0, 0) = (9 * invcss * u * rho + 3*invcss * v * rho - 9 * f1 +
                          9 * f3 + 18 * f7) /
                         18.0;
            f(6, 0, 0) = (-9 * invcss * u * rho + 3*invcss * v * rho + 9 * f1 -
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
            rho = (css * f0 + css * f1 + 2 * css * f2 + css * f3 +
                   2 * css * f5 + 2 * css * f6) /
                  (css + v);
            f(4, 0, 0) = (-6 * invcss * v * rho + 9 * f2) / 9.0;
            f(7, 0, 0) = (-9 * invcss * u * rho - 3*invcss * v * rho + 9 * f1 -
                          9 * f3 + 18 * f5) /
                         18.0;
            f(8, 0, 0) = (9 * invcss * u * rho - 3*invcss * v * rho - 9 * f1 +
                          9 * f3 + 18 * f6) /
                         18.0;
        } break;
        case VG_IPJM_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f6 = f(6, 0, 0);
            Real f3 = f(3, 0, 0);
            rho = macroVars(1, -1);
            f(1, 0, 0) = (6 * invcss * u * rho + 9 * f3) / 9.0;
            f(5, 0, 0) = (1 * invcss * u * rho + 1 * invcss * v * rho) /
                         12.0;
            f(7, 0, 0) = (- 1 * invcss * u * rho - 1 * invcss * v * rho) /
                         12.0;
            f(4, 0, 0) = (-6 * invcss * v * rho + 9 * f2) / 9.0;
            f(8, 0, 0) = (3*invcss * u * rho - 3*invcss * v * rho + 18 * f6) / 18.0;
        } break;
        case VG_IPJP_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f4 = f(4, 0, 0);
            Real f3 = f(3, 0, 0);
            Real f7 = f(7, 0, 0);
            rho = macroVars(1, 1);
            f(1, 0, 0) = (6 * invcss * u * rho + 9 * f3) / 9.0;
            f(5, 0, 0) = (3*invcss * u * rho + 3*invcss * v * rho + 18 * f7) / 18.0;
            f(8, 0, 0) = (-1 * invcss * u * rho + 1 * invcss * v * rho) /
                         12.0;
            f(2, 0, 0) = (6 * invcss * v * rho + 9 * f4) / 9.0;
            f(6, 0, 0) = (1 * invcss * u * rho - 1 * invcss * v * rho) /
                         12.0;
        } break;
        case VG_IMJP_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f8 = f(8, 0, 0);
            Real f4 = f(4, 0, 0);
            rho = macroVars(-1, 1);
            f(5, 0, 0) = (1 * invcss * u * rho + 1 * invcss * v * rho) /
                         12.0;
            f(2, 0, 0) = (6 * invcss * v * rho + 9 * f4) / 9.0;
            f(6, 0, 0) =
                (-3*(invcss * u * rho) + 3*invcss * v * rho + 18 * f8) / 18.0;
            f(3, 0, 0) = (-6 * invcss * u * rho + 9 * f1) / 9.0;
            f(7, 0, 0) = (- 1 * invcss * u * rho - 1 * invcss * v * rho) /
                         12.0;
        } break;
        case VG_IMJM_I: {
            // Knows
            Real f0 = f(0, 0, 0);
            Real f1 = f(1, 0, 0);
            Real f2 = f(2, 0, 0);
            Real f5 = f(5, 0, 0);
            rho = macroVars(-1, -1);
            f(6, 0, 0) = (1 * invcss * u * rho - 1 * invcss * v * rho) /
                         12.0;
            f(3, 0, 0) = (-6 * invcss * u * rho + 9 * f1) / 9.0;
            f(7, 0, 0) =
                (-3*(invcss * u * rho) - 3*invcss * v * rho + 18 * f5) / 18.0;
            f(4, 0, 0) = (-6 * invcss * v * rho + 9 * f2) / 9.0;
            f(8, 0, 0) = (- 1 * invcss * u * rho + 1 * invcss * v * rho) /
                         12.0;
        } break;
        case VG_IPJP_O: {
            Real f7 = f(7, 0, 0);
            f(5, 0, 0)=f7;
        }
        case VG_IPJM_O: {
            Real f8 = f(8, 0, 0);
            f(6, 0, 0)=f8;
        }
        case VG_IMJP_O: {
            Real f6 = f(6, 0, 0);
            f(8, 0, 0)=f6;
        }
        case VG_IMJM_O: {
            Real f5 = f(5, 0, 0);
            f(7, 0, 0)=f5;
        }
        default:
            break;
    }

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
        if (isnan(res) || res <= -2 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n",
                res, xiIdx);
            assert(!(isnan(res) || res <= -2 || isinf(res)));
        }
#endif
    }
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
    
#endif //OPS_3D
}

void KerCutCellEQMDiffuseReflADF3D(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars,const ACC<Real> &acceleration,const Real *dt,
                                const int *lattIdx) {
#ifdef OPS_3D
    // This kernel is suitable for any single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Wall) {
    const int equilibriumOrder{2};
    const VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
    Real w = givenMacroVars[2];
    Real g[]{acceleration(0, 0, 0), acceleration(1, 0, 0), acceleration(2,0,0)};
    Real CS=1/sqrt(3);
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
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        Real cz{XI[xiIdx * LATTDIM + 2]/(CS*CS)};
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
                    CalcBGKFeqFE3D(xiIdx, 1, u, v, w, g, dt, 1, equilibriumOrder);
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
            CalcBGKFeqFE3D(parallel[idx], rhoWall, u, v, w, g, dt, 1, equilibriumOrder);
    }
    for (int idx = 0; idx < numOutgoing; idx++) {
        int xiIdx = outgoing[idx];
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        Real cz{XI[xiIdx * LATTDIM + 2]/(CS*CS)};
        f(xiIdx, 0, 0, 0) =
            f(OPP[xiIdx], 0, 0, 0) +
            2 * rhoWall * WEIGHTS[xiIdx] * (cx * u + cy * v + cz * w);
#ifdef CPU
        const Real res{f(xiIdx, 0, 0, 0)};
        if (isnan(res) || res <= -2 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n",
                res, xiIdx);
            assert(!(isnan(res) || res <= -2 || isinf(res)));
        }
#endif
    }
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
    }
#endif //OPS_3D
}

void KerCutCellEQMDiffuseReflADG3D(ACC<Real> &f, const ACC<int> &nodeType,
                                const ACC<int> &geometryProperty,
                                const Real *givenMacroVars, const ACC<Real> &Conc, const ACC<Real> &Mu,
                                const Real *dt,
                                const int *lattIdx) {
#ifdef OPS_3D
    // This kernel is suitable for any single-speed lattice
    // but only for the second-order expansion at this moment
    // Therefore, the equilibrium function order is fixed at 2

    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Wall) {
    const int equilibriumOrder{2};
    const VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    Real mu{Mu(0,0,0)};
    Real C{Conc(0,0,0)};
    Real u = givenMacroVars[0];
    Real v = givenMacroVars[1];
    Real w = givenMacroVars[2];
    Real CS=1/sqrt(3);
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
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        Real cz{XI[xiIdx * LATTDIM + 2]/(CS*CS)};
        BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
        switch (bdt) {
            case BndryDv_Incoming: {
                incoming[numIncoming] = xiIdx;
                rhoIncoming += CalcBGKGeqFE3D(xiIdx, C, u, v, w, mu, dt, equilibriumOrder);
                numIncoming++;
            } break;
            case BndryDv_Outgoing: {
                outgoing[numOutgoing] = xiIdx;
                deltaRho += CalcBGKGeqFE3D(xiIdx, C, u, v, w, mu, dt, equilibriumOrder);
                numOutgoing++;
            } break;
            case BndryDv_Parallel: {
                parallel[numParallel] = xiIdx;
                if (xiIdx>lattIdx[0]){
                    rhoParallel +=CalcBGKGeqFE3D(xiIdx, C, u, v, w, mu, dt, equilibriumOrder);
                //    CalcBGKGeqFE3D(xiIdx, 1, u, v, w, 1, equilibriumOrder);
                }
                numParallel++;
            } break;
            default:
                break;
        }
    }
    Real rhoWall = rhoIncoming+deltaRho+rhoParallel;
#ifdef CPU
#if DebugLevel >= 2
    ops_printf("Calculated wall density =  %f\n", rhoWall);
#endif
#endif
    for (int idx = 0; idx < numParallel; idx++) {
        if (idx>lattIdx[0]){
        f(parallel[idx], 0, 0, 0) =CalcBGKGeqFE3D(parallel[idx], C, u, v, w, mu, dt, equilibriumOrder);
        }
        else if (idx==lattIdx[0]){
            f(parallel[idx], 0, 0, 0) =C-rhoWall;
        }
    }
    
    for (int idx = 0; idx < numOutgoing; idx++) {
        int xiIdx = outgoing[idx];
        Real cx{XI[xiIdx * LATTDIM]/(CS*CS)};
        Real cy{XI[xiIdx * LATTDIM + 1]/(CS*CS)};
        Real cz{XI[xiIdx * LATTDIM + 2]/(CS*CS)};
        Real opps={OPP[xiIdx]};
                //std::cout<<f(opps,0,0)<<" ";
                f(xiIdx, 0, 0, 0)=CalcBGKGeqFE3D(xiIdx, C, u, v, w, mu, dt, equilibriumOrder);
    }
    /*
#ifdef CPU
        const Real res{f(xiIdx, 0, 0, 0)};
        if (isnan(res) || res <= -10 || isinf(res)) {
            ops_printf(
                "Error! Distribution function %f becomes invalid  at the "
                "lattice %i\n",
                res, xiIdx);
            assert(!(isnan(res) || res <= -10 || isinf(res)));
        }
#endif
    */
    
    
    
    delete[] outgoing;
    delete[] incoming;
    delete[] parallel;
    }
#endif //OPS_3D
}

void KerCutCellZouHeVelocityADG3D(const Real *givenMacroVars,
                             const ACC<int> &nodeType,
                             const ACC<int> &geometryProperty,
                             ACC<Real> &f,
                                const int *lattIdx) {
#ifdef OPS_3D
    /*!
    Note: This boundary condition requires both stream and collision happenning
    at a boundary point.
    Note: This boundary condition is lattice specific.
    */
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Wall) {
    int numOutgoing{0};
    int numIncoming{0};
    int numParallel{0};
    int lattStart{lattIdx[0]};
    int lattEnd{lattIdx[1]};
    Real deltaRho{0};

    const VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0, 0);
    Real rho{0};
    Real u{0};
    Real v{0};
    Real w{0};
    Real css = 1/sqrt(3);
    Real invcss=1/css;
    
    //Real norm[3]={0,0,0};
    std::array<Real,3> norm = {findBdyNorm3D(vg)};
    
    switch (vg) {
        case VG_IP:
        case VG_IM:
        case VG_JP:
        case VG_JM:
        case VG_KP:
        case VG_KM:{
            Real Npref;
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
            BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                Real cjx{XI[xiIdx * LATTDIM]/(css*css)};
                Real cjy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                Real cjz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                Npref+=0.5*f(xiIdx,0,0,0)*(1-abs(cjx*norm[0]+cjy*norm[1]+cjz*norm[2]));
                
                //switch (bdt) {
                //case BndryDv_Incoming: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //case BndryDv_Parallel: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //deafult:
                //    break;
                //}
                
            }
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
                BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                if (bdt==BndryDv_Outgoing) {
                    Real cx{XI[xiIdx * LATTDIM]/(css*css)};
                    Real cy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                    Real cz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                    Real tangx=cx-(cx*norm[0])*norm[0];
                    Real tangy=cx-(cx*norm[2])*norm[1];
                    Real tangz=cx-(cx*norm[2])*norm[2];
                    
                    Real opps={OPP[xiIdx]};
                    Real N;
                    for (int xjIdx = lattStart; xjIdx <= lattEnd; xjIdx++) {
                        Real cjx{XI[xjIdx * LATTDIM]/(css*css)};
                        Real cjy{XI[xjIdx * LATTDIM + 1]/(css*css)};
                        Real cjz{XI[xjIdx * LATTDIM + 2]/(css*css)};
                        N+=tangx*cjx+tangy*cjy+tangz*cjz;
                    }
                    f(xiIdx,0,0,0)=f(opps,0,0,0)-rho/6*(cx*u+cy*v+cz*w)-rho/3*(tangx*u+tangy*v+tangz*w)+Npref*N;
                }
            }
        }break;
        
        case VG_IPJP_I:
        case VG_IPJM_I:
        case VG_IMJP_I:
        case VG_IMJM_I:
        case VG_IPKP_I:
        case VG_IPKM_I:
        case VG_IMKP_I:
        case VG_IMKM_I:
        case VG_JPKP_I:
        case VG_JPKM_I:
        case VG_JMKP_I:
        case VG_JMKM_I:{
            Real fP=0;
            Real Npref=0;
            Real norma[3]={0,0,0};
            Real normb[3]={0,0,0};
            if (norm[0]==0){
                norma[1]=sqrt(2)*norm[1];
                normb[2]=sqrt(2)*norm[2];
            }
            if (norm[1]==0){
                norma[0]=sqrt(2)*norm[0];
                normb[2]=sqrt(2)*norm[2];
            }
            if (norm[2]==0){
                norma[0]=sqrt(2)*norm[0];
                normb[0]=sqrt(2)*norm[2];
            }

            Real normcross[3]={norma[1]*normb[2]-norma[2]*normb[1],norma[2]*normb[0]-norma[0]*normb[2],norma[0]*normb[1]-norma[1]*normb[0]};
            Real dnorm[3]={(norma[0]-normb[0])/sqrt(2),(norma[1]-normb[1])/sqrt(2),(norma[2]-normb[2])/sqrt(2)};
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
            BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                Real cjx{XI[xiIdx * LATTDIM]/(css*css)};
                Real cjy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                Real cjz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                Npref+=0.25*f(xiIdx,0,0,0)*(1-abs(cjx*norma[0]+cjy*norma[1]+cjz*norma[2]))*(1-abs(cjx*normb[0]+cjy*normb[1]+cjz*normb[2]));
                
                if (xiIdx!=lattStart){
                    fP+=1/22*f(xiIdx,0,0,0)*(1-abs(cjx*normcross[0]+cjy*normcross[1]+cjz*normcross[2]))*(1-abs(cjx*dnorm[0]+cjy*dnorm[1]+cjz*dnorm[2]));
                    //std::cout<<f(xiIdx,0,0,0)<<" "<<abs(cjx*normcross[0]+cjy*normcross[1]+cjz*normcross[2])<<" "<<abs(cjx*dnorm[0]+cjy*dnorm[1]+cjz*dnorm[2])<<" ";
                }
                
                
                //switch (bdt) {
                //case BndryDv_Incoming: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //case BndryDv_Parallel: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //deafult:
                //    break;
                //}
                
            
            }
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
                BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                
                if (bdt==BndryDv_Outgoing) {
                    Real cx{XI[xiIdx * LATTDIM]/(css*css)};
                    Real cy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                    Real cz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                    Real tangx=cx-(cx*norm[0])*norm[0];
                    Real tangy=cx-(cx*norm[2])*norm[1];
                    Real tangz=cx-(cx*norm[2])*norm[2];
                    
                    Real opps={OPP[xiIdx]};
                    Real N;
                    for (int xjIdx = lattStart; xjIdx <= lattEnd; xjIdx++) {
                        Real cjx{XI[xjIdx * LATTDIM]/(css*css)};
                        Real cjy{XI[xjIdx * LATTDIM + 1]/(css*css)};
                        Real cjz{XI[xjIdx * LATTDIM + 2]/(css*css)};
                        N+=tangx*cjx+tangy*cjy+tangz*cjz;
                    }
                    f(xiIdx,0,0,0)=f(opps,0,0,0)-rho/6*(cx*u+cy*v+cz*w)-rho/3*(tangx*u+tangy*v+tangz*w)+Npref*N;
                }
                else if (bdt==BndryDv_Parallel) {
                    
                    if (xiIdx==lattStart) {
                        //f(xiIdx,0,0,0)=12*fP;
                    }
                    else {
                        f(xiIdx,0,0,0)=fP;
                    }
                }
            }
                
        }break;
        
        case VG_IPJP_O:
        case VG_IPJM_O:
        case VG_IMJP_O:
        case VG_IMJM_O:
        case VG_IPKP_O:
        case VG_IPKM_O:
        case VG_IMKP_O:
        case VG_IMKM_O:
        case VG_JPKP_O:
        case VG_JPKM_O:
        case VG_JMKP_O:
        case VG_JMKM_O:{
            Real Npref;
            Real norma[3]={0,0,0};
            Real normb[3]={0,0,0};
            if (norm[0]==0){
                norma[1]=sqrt(2)*norm[1];
                normb[2]=sqrt(2)*norm[2];
            }
            if (norm[1]==0){
                norma[0]=sqrt(2)*norm[0];
                normb[2]=sqrt(2)*norm[2];
            }
            if (norm[2]==0){
                norma[0]=sqrt(2)*norm[0];
                normb[0]=sqrt(2)*norm[2];
            }

            Real normcross[3]={norma[1]*normb[2]-norma[2]*normb[1],norma[2]*normb[0]-norma[0]*normb[2],norma[0]*normb[1]-norma[1]*normb[0]};
	    
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
            BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                Real cjx{XI[xiIdx * LATTDIM]/(css*css)};
                Real cjy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                Real cjz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                Npref+=0.25*f(xiIdx,0,0,0)*(cjx*normcross[0]+cjy*normcross[1]+cjz*normcross[2]);
                
                
                //switch (bdt) {
                //case BndryDv_Incoming: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //case BndryDv_Parallel: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //deafult:
                //    break;
                //}
                
            
            }
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
                BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                
                if (bdt==BndryDv_Outgoing) {
		    //std::cout<<xiIdx<<" ";
                    Real cx{XI[xiIdx * LATTDIM]/(css*css)};
                    Real cy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                    Real cz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                    Real pref=(cx*normcross[0]+cy*normcross[1]+cz*normcross[2]);
                    Real tangx=cx-(cx*norm[0])*norm[0];
                    Real tangy=cx-(cx*norm[2])*norm[1];
                    Real tangz=cx-(cx*norm[2])*norm[2];
                    
                    Real opps={OPP[xiIdx]};
                    f(xiIdx,0,0,0)=f(opps,0,0,0)-rho/6*(cx*u+cy*v+cz*w)-rho/3*(tangx*u+tangy*v+tangz*w)+pref*Npref;
		    //std::cout<<f(xiIdx,0,0,0)<<" ";
                }
            }
        }break;
        
        case VG_IPJPKP_I:
        case VG_IPJMKP_I:
        case VG_IMJPKP_I:
        case VG_IMJMKP_I:
        case VG_IPJPKM_I:
        case VG_IPJMKM_I:
        case VG_IMJPKM_I:
        case VG_IMJMKM_I:{
            Real fP=0;

            Real normTot[3]={norm[0],norm[1],norm[2]};
            for (int xiIdx = lattStart+1; xiIdx <= lattEnd; xiIdx++) {
            BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                Real cjx{XI[xiIdx * LATTDIM]/(css*css)};
                Real cjy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                Real cjz{XI[xiIdx * LATTDIM + 2]/(css*css)};
                fP+=1/18*f(xiIdx,0,0,0)*(1-abs(cjx*normTot[0]+cjy*normTot[1]+cjz*normTot[2]));
                //std::cout<<fP<<" ";
                
                //switch (bdt) {
                //case BndryDv_Incoming: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //case BndryDv_Parallel: {
                //    
                //    rho += f(xiIdx, 0, 0, 0);
                //} break;
                //deafult:
                //    break;
                //}
                
            
            }
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
                BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                if (bdt==BndryDv_Parallel) {
                    if (xiIdx!=lattStart) {
                        //f(xiIdx,0,0,0)=12*fP;
			f(xiIdx,0,0,0)=fP;
                    }
                }
            }
	    //std::cout<<" ";
        }

        case VG_IPJPKP_O:
        case VG_IPJMKP_O:
        case VG_IMJPKP_O:
        case VG_IMJMKP_O:
        case VG_IPJPKM_O:
        case VG_IPJMKM_O:
        case VG_IMJPKM_O:
        case VG_IMJMKM_O:{                
            
            for (int xiIdx = lattStart; xiIdx <= lattEnd; xiIdx++) {
                BndryDvType bdt = FindBdyDvType3D(vg, &XI[xiIdx * LATTDIM]);
                if (bdt==BndryDv_Outgoing) {
		    //std::cout<<xiIdx<<" ";
                    Real cx{XI[xiIdx * LATTDIM]/(css*css)};
                    Real cy{XI[xiIdx * LATTDIM + 1]/(css*css)};
                    Real cz{XI[xiIdx * LATTDIM + 2]/(css*css)};                    
                    Real opps={OPP[xiIdx]};
                    Real tangx=cx-(cx*norm[0])*norm[0];
                    Real tangy=cx-(cx*norm[2])*norm[1];
                    Real tangz=cx-(cx*norm[2])*norm[2];

                    f(xiIdx,0,0,0)=f(opps,0,0,0)-rho/6*(cx*u+cy*v+cz*w)-rho/3*(tangx*u+tangy*v+tangz*w);
		    //std::cout<<f(xiIdx,0,0,0)<<" ";
                }
            }
	    //std::cout<<" ";
        }

        
        default:
            break;
    }

    }
#endif  // OPS_3D
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