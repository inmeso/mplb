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

#ifndef CAVITY3D_SWAP_KERNEL_INC
#define CAVITY3D_SWAP_KERNEL_INC

void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            ACC<Real>& w, const ACC<Real>& coordinates,
                            const int* idx) {
    rho(0, 0, 0) = 1;
    u(0, 0, 0) = 0;
    v(0, 0, 0) = 0;
    w(0, 0, 0) = 0;
}

// This kernel function needs lattices sorted in a special order
// see Jonas Latt: Technical report: How to implement your DdQq dynamics with
// only q variables per node (instead of 2q)
void KerSwapCollide3D(ACC<Real>& f, const ACC<int>& nodeType,
                      const ACC<int>& geometry, const ACC<Real>& Rho,
                      const ACC<Real>& U, const ACC<Real>& V,
                      const ACC<Real>& W, const Real* tauRef, const Real* dt,
                      const int* lattIdx) {
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    VertexType vt = (VertexType)nodeType(0, 0, 0);

    Real rho{Rho(0, 0, 0)};
    Real u{U(0, 0, 0)};
    Real v{V(0, 0, 0)};
    Real w{W(0, 0, 0)};
    const Real T{1};
    const int polyOrder{2};
    Real tau = (*tauRef);
    Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        const Real feq{CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder)};
        f(xiIndex, 0, 0, 0) =
            feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0, 0) - feq);
    }

    const int half{(lattIdx[1] - lattIdx[0]) / 2};

    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
        vt == VertexType::VirtualBoundary) {
        for (int xiIndex = lattIdx[0] + 1; xiIndex <= (lattIdx[0] + half);
             xiIndex++) {
            Swap(f(xiIndex, 0, 0, 0), f(xiIndex + half, 0, 0, 0));
        }
    }
    if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid &&
        vt != VertexType::VirtualBoundary && vt != VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0] + 1; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            int cz = (int)XI[xiIndex * LATTDIM + 2];
            switch (vg) {
                    // faces six types
                case VG_IP: {
                    if (cx == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cx < 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                } break;
                case VG_IM: {
                    if (cx == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cx > 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                }

                break;
                case VG_JP: {
                    if (cy == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cy < 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                } break;
                case VG_JM: {
                    if (cy == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cy > 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                } break;
                case VG_KP: {
                    if (cz == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cz < 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }

                } break;
                case VG_KM: {
                    if (cz == 0) {
                        Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                    }
                    if (cz > 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                } break;
                    faces six types end 12 edges case VG_IPJP_I : {
                        if (cy <= 0 && cx <= 0) {
                            if (cy == 0 && cx == 0) {
                                Swap(f(xiIndex, 0, 0, 0),
                                     f(OPP[xiIndex], 0, 0, 0));
                            } else {
                                f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                            }
                        }
                    }
                    break;
                case VG_IPJM_I:
                    if (cy >= 0 && cx <= 0) {
                        if (cy == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_IMJP_I:
                    if (cy <= 0 && cx >= 0) {
                        if (cy == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_IMJM_I:
                    if (cy >= 0 && cx >= 0) {
                        if (cy == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                // k
                case VG_IPKP_I:
                    if (cz <= 0 && cx <= 0) {
                        if (cz == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_IPKM_I:
                    if (cz >= 0 && cx <= 0) {
                        if (cz == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_IMKP_I:
                    if (cz <= 0 && cx >= 0) {
                        if (cz == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_IMKM_I:
                    if (cz >= 0 && cx >= 0) {
                        if (cz == 0 && cx == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_JPKP_I:
                    if (cz <= 0 && cy <= 0) {
                        if (cz == 0 && cy == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_JPKM_I:
                    if (cz >= 0 && cy <= 0) {
                        if (cz == 0 && cy == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_JMKP_I:
                    if (cz <= 0 && cy >= 0) {
                        if (cz == 0 && cy == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                case VG_JMKM_I:
                    if (cz >= 0 && cy >= 0) {
                        if (cz == 0 && cy == 0) {
                            Swap(f(xiIndex, 0, 0, 0), f(OPP[xiIndex], 0, 0, 0));
                        } else {
                            f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                        }
                    }
                    break;
                // k_out end
                // 12 edges end
                // 8 corners
                // inner corners
                case VG_IPJPKP_I:
                    if (cx <= 0 && cy <= 0 && cz <= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IPJPKM_I:
                    if (cx <= 0 && cy <= 0 && cz >= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IPJMKP_I:
                    if (cx <= 0 && cy >= 0 && cz <= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IPJMKM_I:
                    if (cx <= 0 && cy >= 0 && cz >= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IMJPKP_I:
                    if (cx >= 0 && cy <= 0 && cz <= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IMJPKM_I:
                    if (cx >= 0 && cy <= 0 && cz >= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IMJMKP_I:
                    if (cx >= 0 && cy >= 0 && cz <= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                case VG_IMJMKM_I:
                    if (cx >= 0 && cy >= 0 && cz >= 0) {
                        f(xiIndex, 0, 0, 0) = f(OPP[xiIndex], 0, 0, 0);
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

void KerSwapStream3D(ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<int>& geometry, const int* lattIdx) {
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    VertexType vt = (VertexType)nodeType(0, 0, 0);

    const int half{(lattIdx[1] - lattIdx[0]) / 2};

    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
        vt == VertexType::VirtualBoundary) {
        for (int xiIndex = lattIdx[0] + 1; xiIndex <= (lattIdx[0] + half);
             xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            int cz = (int)XI[xiIndex * LATTDIM + 2];
            Swap(f(xiIndex + half, 0, 0, 0), f(xiIndex, cx, cy, cz));
        }
    }
    if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid &&
        vt != VertexType::VirtualBoundary && vt != VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0] + 1; xiIndex <= (lattIdx[0] + half);
             xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            int cz = (int)XI[xiIndex * LATTDIM + 2];
            // opposite lattice velocity
            int rcx = (int)XI[OPP[xiIndex] * LATTDIM];
            int rcy = (int)XI[OPP[xiIndex] * LATTDIM + 1];
            int rcz = (int)XI[OPP[xiIndex] * LATTDIM + 2];
            switch (vg) {
                    // faces six types
                case VG_IP:
                    if (rcx <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IM:
                    if (rcx >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JP:
                    if (rcy <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JM:
                    if (rcy >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_KP:
                    // (cx=0 means stream is implemented at i=0,so here we
                    //  disable the step at boundary)
                    if (rcz <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_KM:
                    if (rcz >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                // faces six types end
                // 12 edges
                case VG_IPJP_I:
                    if (rcy <= 0 && rcx <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IPJM_I:
                    if (rcy >= 0 && rcx <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJP_I:
                    if (rcy <= 0 && rcx >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJM_I:
                    if (rcy >= 0 && rcx >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                // k
                case VG_IPKP_I:
                    if (rcz <= 0 && rcx <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IPKM_I:
                    if (rcz >= 0 && rcx <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMKP_I:
                    if (rcz <= 0 && rcx >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMKM_I:
                    if (rcz >= 0 && rcx >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JPKP_I:
                    if (rcz <= 0 && rcy <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JPKM_I:
                    if (rcz >= 0 && rcy <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JMKP_I:
                    if (rcz <= 0 && rcy >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_JMKM_I:
                    if (rcz >= 0 && rcy >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                // k_out end
                // 12 edges end
                // 8 corners
                // inner corners
                case VG_IPJPKP_I:
                    if (rcx <= 0 && rcy <= 0 && rcz <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IPJPKM_I:
                    if (rcx <= 0 && rcy <= 0 && rcz >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IPJMKP_I:
                    if (rcx <= 0 && rcy >= 0 && rcz <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IPJMKM_I:
                    if (rcx <= 0 && rcy >= 0 && rcz >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJPKP_I:
                    if (rcx >= 0 && rcy <= 0 && rcz <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJPKM_I:
                    if (rcx >= 0 && rcy <= 0 && rcz >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJMKP_I:
                    if (rcx >= 0 && rcy >= 0 && rcz <= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                case VG_IMJMKM_I:
                    if (rcx >= 0 && rcy >= 0 && rcz >= 0) {
                        Swap(f(xiIndex + half, 0, 0, 0),
                             f(xiIndex, cx, cy, cz));
                    }
                    break;
                default:
                    break;
            }
        }
    }
}
#endif  // CAVITY3D_SWAP_KERNEL_INC
