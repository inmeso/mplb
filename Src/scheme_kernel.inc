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

/*！ @brief Define kernels for numerical schemes.
 *  @author Jianping Meng
 *  @details Including various space and time scheme.
 *  For the stream-collision, we plan only to support the standard lattice,
 *  That is the speed is only 1, not multi-speed lattice
 *  For the finite-difference scheme, only the cut-cell mesh technique is
 *  supported currently
 **/

#ifndef SCHEME_KERNEL_INC
#define SCHEME_KERNEL_INC
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "boundary.h"
#include "model.h"
#ifdef OPS_2D  // two dimensional code

void KerStream(ACC<Real>&f, const ACC<Real>&fStage, const ACC<int>&nodeType,
    const ACC<int>&geometry, const int *lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        int cx = (int)XI[xiIndex * LATTDIM];
        int cy = (int)XI[xiIndex * LATTDIM + 1];

        if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
            vt == VertexType::VirtualBoundary) {
            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            continue;
        }
        // Block boundary and immersed boundary
        if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid &&
            vt != VertexType::VirtualBoundary && vt != VertexType::MDPeriodic) {
            // this is reserved in case some boundary conditions
            // might not allow stream
            bool streamRequired{true};

            if (streamRequired) {
                if ((cx == 0) && (cy == 0)) {
                    f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
                    continue;
                }
            }
            switch (vg) {
                case VG_IP:
                    // (cx=0 means stream is implemented at i=0,so here we
                    //  disable the step at boundary)
                    if (streamRequired) {
                        if (cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IM:
                    if (streamRequired) {
                        if (cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_JP:
                    if (streamRequired) {
                        if (cy <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_JM:
                    if (streamRequired) {
                        if (cy >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 && cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 && cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 && cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 && cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 || cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 || cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 || cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 || cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }
#endif  // OPS_2D
}

// void KerCutCellCVTUpwind1st(const ACC<Real>& coordinateXYZ,
//                             const ACC<int>& nodeType, const ACC<int>&
//                             geometry, const ACC<Real>& f, ACC<Real>&
//                             fGradient) {
// #ifdef OPS_2D
//     VertexType vt = (VertexType)nodeType(0, 0);
//     VertexGeometryTypes vg = (VertexGeometryTypes)geometry(0, 0);
//     for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
//         for (int xiIndex = COMPOINDEX[2 * compoIndex];
//              xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
//             Real cx{XI[xiIndex * LATTDIM]};
//             Real cy{XI[xiIndex * LATTDIM + 1]};
//             // setting a initial value
//             fGradient(xiIndex, 0, 0) = 0;
//             // make sure no calculation occurring at boundary when
//             unnecessary
//             // this can avoid access undefined memory if the halo point is
//             set
//             // incorrectly
//             bool needCalc{true};
//             if (vt == VertexType::ImmersedSolid) {
//                 needCalc = false;
//             }
//             // Block and immersed boundary nodes
//             if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid) {
//                 switch (vg) {
//                     case VG_IP:
//                         if (cx > 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IM:
//                         if (cx < 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_JP:
//                         if (cy > 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_JM:
//                         if (cy < 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IPJP_I:
//                         if (cy > 0 || cx > 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IPJM_I:
//                         if (cy < 0 || cx > 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IMJP_I:
//                         if (cy > 0 || cx < 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IMJM_I:
//                         if (cy < 0 || cx < 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     // It appears that the information is well defined for
//                     the
//                     // outer corners
//                     // which is quite different from the inner corners
//                     case VG_IPJP_O:
//                         if (cy >= 0 && cx >= 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IPJM_O:
//                         if (cy <= 0 && cx >= 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IMJP_O:
//                         if (cy >= 0 && cx <= 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     case VG_IMJM_O:
//                         if (cy <= 0 && cx <= 0) {
//                             needCalc = false;
//                         }
//                         break;
//                     default:
//                         break;
//                 }
//             }
//             if (needCalc) {
//                 if (cx > 0) {
//                     fGradient(xiIndex, 0, 0) +=
//                         CS * cx *
//                         ((f(xiIndex, 0, 0) - f(xiIndex, -1, 0)) /
//                          (coordinateXYZ(0, 0, 0) - coordinateXYZ(0, -1, 0)));
//                 }
//                 if (cx < 0) {
//                     fGradient(xiIndex, 0, 0) +=
//                         CS * cx *
//                         ((f(xiIndex, 1, 0) - f(xiIndex, 0, 0)) /
//                          (coordinateXYZ(0, 1, 0) - coordinateXYZ(0, 0, 0)));
//                 }
//                 if (cy > 0) {
//                     fGradient(xiIndex, 0, 0) +=
//                         CS * cy *
//                         ((f(xiIndex, 0, 0) - f(xiIndex, 0, -1)) /
//                          (coordinateXYZ(1, 0, 0) - coordinateXYZ(1, 0, -1)));
//                 }
//                 if (cy < 0) {
//                     fGradient(xiIndex, 0, 0) +=
//                         CS * cy *
//                         ((f(xiIndex, 0, 1) - f(xiIndex, 0, 0)) /
//                          (coordinateXYZ(1, 0, 1) - coordinateXYZ(1, 0, 0)));
//                 }
//             }
//         }
//     }
// #endif  // OPS_2D
// }

// void KerCutCellCVTUpwind2nd(const ACC<Real>& coordinateXYZ,
//                             const ACC<int> nodeType, const ACC<int>&
//                             geometry, const ACC<Real>& f, ACC<Real>&
//                             fGradient) {
// #ifdef OPS_2D
//     VertexType vt = (VertexType)nodeType(0, 0);
//     VertexGeometryTypes vg = (VertexGeometryTypes)geometry(0, 0);
//     for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
//         for (int xiIndex = COMPOINDEX[2 * compoIndex];
//              xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
//             Real cx{XI[xiIndex * LATTDIM]};
//             Real cy{XI[xiIndex * LATTDIM + 1]};
//             bool reduceOrderX{false};
//             bool reduceOrderY{false};
//             bool needCalc{true};
//             // setting a initial value
//             fGradient(xiIndex, 0, 0) = 0;
//             if (vt == VertexType::ImmersedSolid) {
//                 needCalc = false;
//             }
//             // make sure no calculation occurring at boundary when
//             unnecessary
//             // this can avoid access undefined memory if the halo point is
//             set
//             // incorrectly
//             if (vt == VertexType::Fluid) {
//                 // current node is a fluid point
//                 if (cx > 0) {
//                     VertexType vtUpwind = (VertexType)nodeType(-1, 0);
//                     VertexGeometryTypes vgUpwind =
//                         (VertexGeometryTypes)geometry(-1, 0);
//                     if (vtUpwind != VertexType::Fluid && vgUpwind == VG_IP) {
//                         // if the upwind node is boundary and VG_IP point
//                         reduceOrderX = true;
//                     }
//                 }  // cx >0
//                 if (cx < 0) {
//                     VertexType vtUpwind = (VertexType)nodeType(1, 0);
//                     VertexGeometryTypes vgUpwind =
//                         (VertexGeometryTypes)geometry(1, 0);
//                     if (vtUpwind != VertexType::Fluid && vgUpwind == VG_IM) {
//                         // if the upwind node is boundary and VG_IM point
//                         reduceOrderX = true;
//                     }
//                 }  // Cx<0
//                 if (cy > 0) {
//                     VertexType vtUpwind = (VertexType)nodeType(0, -1);
//                     VertexGeometryTypes vgUpwind =
//                         (VertexGeometryTypes)geometry(0, -1);
//                     if (vtUpwind != VertexType::Fluid && vgUpwind == VG_JP) {
//                         reduceOrderY = true;
//                     }
//                 }  // cy >0
//                 if (cy < 0) {
//                     VertexType vtUpwind = (VertexType)nodeType(0, 1);
//                     VertexGeometryTypes vgUpwind =
//                         (VertexGeometryTypes)geometry(0, 1);
//                     if (vtUpwind != VertexType::Fluid && vgUpwind == VG_JM) {
//                         reduceOrderY = true;
//                     }
//                 }  // cy < 0
//             }      // current node is a fluid point
//             // Current node is a boundary node or a immersed boundary node
//             if (vt != VertexType::Fluid && vt != VertexType::ImmersedSolid) {
//                 switch (vg) {
//                     case VG_IP:
//                         if (cx > 0) {
//                             needCalc = false;
//                         } else {
//                             if (cy > 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(0, -1);
//                                 if (((vgUpwind == VG_IPJM_O ||
//                                       vgUpwind == VG_IPJP_O) &&
//                                      (((VertexType)nodeType(0, -2)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IPJM_I) ||
//                                     (vgUpwind == VG_IPJP_I)) {
//                                     reduceOrderY = true;
//                                 }
//                             }  // cy > 0 VG_IP
//                             if (cy < 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(0, 1);
//                                 if (((vgUpwind == VG_IPJM_O ||
//                                       vgUpwind == VG_IPJP_O) &&
//                                      (((VertexType)nodeType(0, 2)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IPJM_I) ||
//                                     (vgUpwind == VG_IPJP_I)) {
//                                     reduceOrderY = true;
//                                 }
//                             }  // Cy < 0 VG_IP
//                         }
//                         break;
//                     case VG_IM:
//                         if (cx < 0) {
//                             needCalc = false;
//                         } else {
//                             if (cy > 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(0, -1);
//                                 if (((vgUpwind == VG_IMJM_O ||
//                                       vgUpwind == VG_IMJP_O) &&
//                                      (((VertexType)nodeType(0, -2)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJM_I) ||
//                                     (vgUpwind == VG_IMJP_I)) {
//                                     reduceOrderY = true;
//                                 }
//                             }  // cy > 0 VG_IM
//                             if (cy < 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(0, 1);
//                                 if (((vgUpwind == VG_IMJM_O ||
//                                       vgUpwind == VG_IMJP_O) &&
//                                      (((VertexType)nodeType(0, 2)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJM_I) ||
//                                     (vgUpwind == VG_IMJP_I)) {
//                                     reduceOrderY = true;
//                                 }
//                             }  // cy < 0 VG_IM
//                         }
//                         break;
//                     case VG_JP:
//                         if (cy > 0) {
//                             needCalc = false;
//                         } else {
//                             if (cx > 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(-1, 0);
//                                 if (((vgUpwind == VG_IMJP_O ||
//                                       vgUpwind == VG_IPJP_O) &&
//                                      (((VertexType)nodeType(-2, 0)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJP_I) ||
//                                     (vgUpwind == VG_IPJP_I)) {
//                                     reduceOrderX = true;
//                                 }
//                             }  // cx > 0  VG_JP
//                             if (cx < 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(1, 0);
//                                 if (((vgUpwind == VG_IMJP_O ||
//                                       vgUpwind == VG_IPJP_O) &&
//                                      (((VertexType)nodeType(2, 0)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJP_I) ||
//                                     (vgUpwind == VG_IPJP_I)) {
//                                     reduceOrderX = true;
//                                 }
//                             }  // cx <0 VG_JP
//                         }
//                         break;
//                     case VG_JM:
//                         if (cy < 0) {
//                             needCalc = false;
//                         } else {
//                             if (cx > 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(-1, 0);
//                                 if (((vgUpwind == VG_IMJM_O ||
//                                       vgUpwind == VG_IPJM_O) &&
//                                      (((VertexType)nodeType(-2, 0)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJM_I) ||
//                                     (vgUpwind == VG_IPJM_I)) {
//                                     reduceOrderX = true;
//                                 }
//                             }  // cx > 0 VG_JM
//                             if (cx < 0) {
//                                 VertexGeometryTypes vgUpwind =
//                                     (VertexGeometryTypes)geometry(1, 0);
//                                 if (((vgUpwind == VG_IPJM_O ||
//                                       vgUpwind == VG_IMJM_O) &&
//                                      (((VertexType)nodeType(2, 0)) ==
//                                       VertexType::ImmersedSolid)) ||
//                                     (vgUpwind == VG_IMJM_I) ||
//                                     (vgUpwind == VG_IPJM_I)) {
//                                     reduceOrderX = true;
//                                 }
//                             }  // cx > 0 VG_JM
//                         }
//                         break;
//                     case VG_IPJP_I:
//                         if (cy > 0 || cx > 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(1, 0);
//                             if (vgUpwind == VG_IMJP_I) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, 1);
//                             if (vgUpwind == VG_IPJM_I) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IPJM_I:
//                         if (cy < 0 || cx > 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(1, 0);
//                             if (vgUpwind == VG_IMJM_I) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, -1);
//                             if (vgUpwind == VG_IPJP_I) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IMJP_I:
//                         if (cy > 0 || cx < 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(-1, 0);
//                             if (vgUpwind == VG_IPJP_I) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, 1);
//                             if (vgUpwind == VG_IMJM_I) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IMJM_I:
//                         if (cy < 0 || cx < 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(-1, 0);
//                             if (vgUpwind == VG_IPJM_I) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, -1);
//                             if (vgUpwind == VG_IMJP_I) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IPJP_O:
//                         if (cy >= 0 && cx >= 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(1, 0);
//                             if (vgUpwind == VG_IM) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, 1);
//                             if (vgUpwind == VG_JM) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IPJM_O:
//                         if (cy <= 0 && cx >= 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(1, 0);
//                             if (vgUpwind == VG_IM) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, -1);
//                             if (vgUpwind == VG_JP) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IMJP_O:
//                         if (cy >= 0 && cx <= 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(-1, 0);
//                             if (vgUpwind == VG_IP) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, 1);
//                             if (vgUpwind == VG_JM) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     case VG_IMJM_O:
//                         if (cy <= 0 && cx <= 0) {
//                             needCalc = false;
//                         } else {
//                             VertexGeometryTypes vgUpwind =
//                                 (VertexGeometryTypes)geometry(-1, 0);
//                             if (vgUpwind == VG_IP) {
//                                 reduceOrderX = true;
//                             }
//                             vgUpwind = (VertexGeometryTypes)geometry(0, -1);
//                             if (vgUpwind == VG_JP) {
//                                 reduceOrderY = true;
//                             }
//                         }
//                         break;
//                     default:
//                         break;
//                 }
//             }  // if the current node is a boundary node
//             if (needCalc) {
//                 if (cx > 0) {
//                     if (reduceOrderX) {
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cx *
//                             ((f(xiIndex, 0, 0) - f(xiIndex, -1, 0)) /
//                              (coordinateXYZ(0, 0, 0) -
//                               coordinateXYZ(0, -1, 0)));
//                     } else {
//                         Real x_2{coordinateXYZ(0, -2, 0)};
//                         Real x_1{coordinateXYZ(0, -1, 0)};
//                         Real x_0{coordinateXYZ(0, 0, 0)};
//                         // ops_printf("x_2=%f x_1=%f x_0=%f\n",x_2,x_1,x_0);
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cx *
//                             (f(xiIndex, -2, 0) *
//                                  (1 / (x_2 - x_0) - 1 / (x_2 - x_1)) +
//                              f(xiIndex, -1, 0) *
//                                  (1 / (x_1 - x_0) - 1 / (x_1 - x_2)) +
//                              f(xiIndex, 0, 0) *
//                                  (1 / (x_0 - x_1) + 1 / (x_0 - x_2)));
//                     }
//                 }
//                 if (cx < 0) {
//                     if (reduceOrderX) {
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cx *
//                             ((f(xiIndex, 1, 0) - f(xiIndex, 0, 0)) /
//                              (coordinateXYZ(0, 1, 0) - coordinateXYZ(0, 0,
//                              0)));
//                     } else {
//                         Real xp2{coordinateXYZ(0, 2, 0)};
//                         Real xp1{coordinateXYZ(0, 1, 0)};
//                         Real x0{coordinateXYZ(0, 0, 0)};
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cx *
//                             (f(xiIndex, 0, 0) *
//                                  (1 / (x0 - xp2) + 1 / (x0 - xp1)) +
//                              f(xiIndex, 1, 0) *
//                                  (1 / (xp2 - xp1) + 1 / (xp1 - x0)) +
//                              f(xiIndex, 2, 0) *
//                                  (1 / (xp2 - x0) + 1 / (xp1 - xp2)));
//                     }
//                 }
//                 if (cy > 0) {
//                     if (reduceOrderY) {
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cy *
//                             ((f(xiIndex, 0, 0) - f(xiIndex, 0, -1)) /
//                              (coordinateXYZ(1, 0, 0) -
//                               coordinateXYZ(1, 0, -1)));
//                     } else {
//                         Real y_2{coordinateXYZ(1, 0, -2)};
//                         Real y_1{coordinateXYZ(1, 0, -1)};
//                         Real y_0{coordinateXYZ(1, 0, 0)};
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cy *
//                             (f(xiIndex, 0, -2) *
//                                  (1 / (y_2 - y_0) - 1 / (y_2 - y_1)) +
//                              f(xiIndex, 0, -1) *
//                                  (1 / (y_1 - y_0) - 1 / (y_1 - y_2)) +
//                              f(xiIndex, 0, 0) *
//                                  (1 / (y_0 - y_1) + 1 / (y_0 - y_2)));
//                     }
//                 }
//                 if (cy < 0) {
//                     if (reduceOrderY) {
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cy *
//                             ((f(xiIndex, 0, 1) - f(xiIndex, 0, 0)) /
//                              (coordinateXYZ(1, 0, 1) - coordinateXYZ(1, 0,
//                              0)));
//                     } else {
//                         Real yp2{coordinateXYZ(1, 0, 2)};
//                         Real yp1{coordinateXYZ(1, 0, 1)};
//                         Real y0{coordinateXYZ(1, 0, 0)};
//                         fGradient(xiIndex, 0, 0) +=
//                             CS * cy *
//                             (f(xiIndex, 0, 0) *
//                                  (1 / (y0 - yp2) + 1 / (y0 - yp1)) +
//                              f(xiIndex, 0, 1) *
//                                  (1 / (yp2 - yp1) + 1 / (yp1 - y0)) +
//                              f(xiIndex, 0, 2) *
//                                  (1 / (yp2 - y0) + 1 / (yp1 - yp2)));
//                     }
//                 }
//             }
//         }
//     }
// #endif  // OPS_2D
// }

// void KerCutCellSemiImplicitTimeMach(
//     const Real* dt, const Real* schemeCoeff, const ACC<int>& nodeType,
//     const ACC<int>& geometry, const ACC<Real>& fGradient, const ACC<Real>&
//     feq, const ACC<Real>& relaxationTime, const ACC<Real>& bodyForce, Real*
//     f) {
// #ifdef OPS_2D
//     /*
//     dt(0), schemeCoeff(1),nodeType(2),geometry(3), fGradient (4),
//     feq(5),
//     relaxationTime(6), bodyForce(7), f(8)
//     */
//     VertexType vt = (VertexType)nodeType(0, 0);
//     VertexGeometryTypes vg = (VertexGeometryTypes)geometry(0, 0);
//     for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
//         Real tau = relaxationTime(compoIndex, 0, 0);
//         for (int xiIndex = COMPOINDEX[2 * compoIndex];
//              xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
//             Real cx{XI[xiIndex * LATTDIM]};
//             Real cy{XI[xiIndex * LATTDIM + 1]};
//             bool needMarch{true};
//             if (vt == VertexType::ImmersedSolid) {
//                 needMarch = false;
//             }
//             // Block or immersed boundary
//             if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid) {
//                 switch (vg) {
//                     case VG_IP:
//                         if (cx > 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IM:
//                         if (cx < 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_JP:
//                         if (cy > 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_JM:
//                         if (cy < 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IPJP_I:
//                         if (cy > 0 || cx > 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IPJM_I:
//                         if (cy < 0 || cx > 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IMJP_I:
//                         if (cy > 0 || cx < 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IMJM_I:
//                         if (cy < 0 || cx < 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IPJP_O:
//                         if (cy >= 0 && cx >= 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IPJM_O:
//                         if (cy <= 0 && cx >= 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IMJP_O:
//                         if (cy >= 0 && cx <= 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     case VG_IMJM_O:
//                         if (cy <= 0 && cx <= 0) {
//                             needMarch = false;
//                         }
//                         break;
//                     default:
//                         break;
//                 }
//             }
//             if (needMarch) {
//                 f(xiIndex, 0, 0) = (feq(xiIndex, 0, 0) * (*dt) +
//                                     bodyForce(xiIndex, 0, 0) * (*dt) * tau -
//                                     fGradient(xiIndex, 0, 0) * (*dt) * tau +
//                                     f(xiIndex, 0, 0) * tau) /
//                                    ((*dt) + tau);
//             }
//         }
//     }
// #endif  // OPS_2D
// }

// void KerCutCellExplicitTimeMach(
//     const Real* dt, const Real* schemeCoeff, const ACC<int>& nodeType,
//     const ACC<int>& geometry, const ACC<Real>& fGradient, const ACC<Real>&
//     feq, const ACC<Real>& relaxationTime, const ACC<Real>& bodyForce,
//     ACC<Real>& f) {
// #ifdef OPS_2D
//     /*
//     dt(0), schemeCoeff(1),nodeType(2),geometry(3), fGradient (4),
//     feq(5),
//     relaxationTime(6), bodyForce(7), f(8)
//     */
//     VertexType vt = (VertexType)nodeType(0, 0);
//     VertexGeometryTypes vg = (VertexGeometryTypes)geometry(0, 0);
//     for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
//         Real tau = relaxationTime(compoIndex, 0, 0);
//         for (int xiIndex = COMPOINDEX[2 * compoIndex];
//              xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
//             Real cx{XI[xiIndex * LATTDIM]};
//             Real cy{XI[xiIndex * LATTDIM + 1]};
//             bool needMarch{true};
//             if (vt == VertexType::ImmersedSolid) {
//                 needMarch = false;
//             }
//             // Block or immersed boundary
//             if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid) {
//                 {
//                     switch (vg) {
//                         case VG_IP:
//                             if (cx > 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IM:
//                             if (cx < 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_JP:
//                             if (cy > 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_JM:
//                             if (cy < 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IPJP_I:
//                             if (cy > 0 || cx > 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IPJM_I:
//                             if (cy < 0 || cx > 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IMJP_I:
//                             if (cy > 0 || cx < 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IMJM_I:
//                             if (cy < 0 || cx < 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IPJP_O:
//                             if (cy >= 0 && cx >= 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IPJM_O:
//                             if (cy <= 0 && cx >= 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IMJP_O:
//                             if (cy >= 0 && cx <= 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         case VG_IMJM_O:
//                             if (cy <= 0 && cx <= 0) {
//                                 needMarch = false;
//                             }
//                             break;
//                         default:
//                             break;
//                     }
//                 }
//                 if (needMarch) {
//                     f(xiIndex, 0, 0) =
//                         f(xiIndex, 0, 0) +
//                         (*dt) * (*schemeCoeff) *
//                             (bodyForce(xiIndex, 0, 0) -
//                              fGradient(xiIndex, 0, 0) +
//                              (feq(xiIndex, 0, 0) - f(xiIndex, 0, 0)) / tau);
//                 }
//             }
//         }
// #endif  // OPS_2D
//     }

#endif  // OPS_2D outter

#ifdef OPS_3D  // three dimensional code

void KerStream3D(ACC<Real>& f, const ACC<Real>& fStage,
                 const ACC<int>& nodeType, const ACC<int>& geometry,
                 const int* lattIdx) {
#ifdef OPS_3D
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    VertexType vt = (VertexType)nodeType(0, 0, 0);

    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        int cx = (int)XI[xiIndex * LATTDIM];
        int cy = (int)XI[xiIndex * LATTDIM + 1];
        int cz = (int)XI[xiIndex * LATTDIM + 2];

        if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
            vt == VertexType::VirtualBoundary) {
            f(xiIndex, 0, 0, 0) = fStage(xiIndex, -cx, -cy, -cz);
            continue;
        }

        if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid &&
            vt != VertexType::VirtualBoundary && vt != VertexType::MDPeriodic) {
            bool streamRequired{true};
            if (streamRequired) {
                if ((cx == 0) && (cy == 0) && (cz == 0)) {
                    f(xiIndex, 0, 0, 0) = fStage(xiIndex, 0, 0, 0);
                    continue;
                }
            }

            switch (vg) {
                    // faces six types
                case VG_IP:
                    // (cx=0 means stream is implemented at i=0,so here we
                    //  disable the step at boundary)
                    if (streamRequired) {
                        if (cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IM:
                    if (streamRequired) {
                        if (cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JP:
                    if (streamRequired) {
                        if (cy <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JM:
                    if (streamRequired) {
                        if (cy >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_KP:
                    // (cx=0 means stream is implemented at i=0,so here we
                    //  disable the step at boundary)
                    if (streamRequired) {
                        if (cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_KM:
                    if (streamRequired) {
                        if (cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }

                    } else {
                        if (cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                // faces six types end
                // 12 edges
                case VG_IPJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy < 0 && cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy > 0 && cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy < 0 && cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy > 0 && cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                // k
                case VG_IPKP_I:
                    if (streamRequired) {
                        if (cz <= 0 && cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 && cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPKM_I:
                    if (streamRequired) {
                        if (cz >= 0 && cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 && cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMKP_I:
                    if (streamRequired) {
                        if (cz <= 0 && cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 && cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMKM_I:
                    if (streamRequired) {
                        if (cz >= 0 && cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 && cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JPKP_I:
                    if (streamRequired) {
                        if (cz <= 0 && cy <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 && cy < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JPKM_I:
                    if (streamRequired) {
                        if (cz >= 0 && cy <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 && cy < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JMKP_I:
                    if (streamRequired) {
                        if (cz <= 0 && cy >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 && cy > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JMKM_I:
                    if (streamRequired) {
                        if (cz >= 0 && cy >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 && cy > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                // K
                // k_out
                case VG_IPJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy < 0 || cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy > 0 || cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy < 0 || cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cy > 0 || cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                    // IK
                case VG_IPKP_O:
                    if (streamRequired) {
                        if (cz <= 0 || cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 || cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPKM_O:
                    if (streamRequired) {
                        if (cz >= 0 || cx <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 || cx < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMKP_O:
                    if (streamRequired) {
                        if (cz <= 0 || cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 || cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMKM_O:
                    if (streamRequired) {
                        if (cz >= 0 || cx >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 || cx > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                    // JK
                case VG_JPKP_O:
                    if (streamRequired) {
                        if (cz <= 0 || cy <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 || cy < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JPKM_O:
                    if (streamRequired) {
                        if (cz >= 0 || cy <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 || cy < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JMKP_O:
                    if (streamRequired) {
                        if (cz <= 0 || cy >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz < 0 || cy > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_JMKM_O:
                    if (streamRequired) {
                        if (cz >= 0 || cy >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cz > 0 || cy > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                // k_out end
                // 12 edges end
                // 8 corners
                // inner corners
                case VG_IPJPKP_I:
                    if (streamRequired) {
                        if (cx <= 0 && cy <= 0 && cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 && cy < 0 && cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJPKM_I:
                    if (streamRequired) {
                        if (cx <= 0 && cy <= 0 && cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 && cy < 0 && cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJMKP_I:
                    if (streamRequired) {
                        if (cx <= 0 && cy >= 0 && cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 && cy > 0 && cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJMKM_I:
                    if (streamRequired) {
                        if (cx <= 0 && cy >= 0 && cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 && cy > 0 && cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJPKP_I:
                    if (streamRequired) {
                        if (cx >= 0 && cy <= 0 && cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 && cy < 0 && cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJPKM_I:
                    if (streamRequired) {
                        if (cx >= 0 && cy <= 0 && cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 && cy < 0 && cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJMKP_I:
                    if (streamRequired) {
                        if (cx >= 0 && cy >= 0 && cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 && cy > 0 && cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJMKM_I:
                    if (streamRequired) {
                        if (cx >= 0 && cy >= 0 && cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 && cy > 0 && cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                // out corner
                case VG_IPJPKP_O:
                    if (streamRequired) {
                        if (cx <= 0 || cy <= 0 || cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 || cy < 0 || cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJPKM_O:
                    if (streamRequired) {
                        if (cx <= 0 || cy <= 0 || cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 || cy < 0 || cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJMKP_O:
                    if (streamRequired) {
                        if (cx <= 0 || cy >= 0 || cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 || cy > 0 || cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IPJMKM_O:
                    if (streamRequired) {
                        if (cx <= 0 || cy >= 0 || cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx < 0 || cy > 0 || cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJPKP_O:
                    if (streamRequired) {
                        if (cx >= 0 || cy <= 0 || cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 || cy < 0 || cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJPKM_O:
                    if (streamRequired) {
                        if (cx >= 0 || cy <= 0 || cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 || cy < 0 || cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJMKP_O:
                    if (streamRequired) {
                        if (cx >= 0 || cy >= 0 || cz <= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 || cy > 0 || cz < 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                case VG_IMJMKM_O:
                    if (streamRequired) {
                        if (cx >= 0 || cy >= 0 || cz >= 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    } else {
                        if (cx > 0 || cy > 0 || cz > 0) {
                            f(xiIndex, 0, 0, 0) =
                                fStage(xiIndex, -cx, -cy, -cz);
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }

#endif  // OPS_3D
}

#endif  // OPS_3D outter

#endif  // SCHEME_KERNEL.inc
