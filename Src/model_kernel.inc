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

/*! @brief   Define kernel functions related to discrete velocity model
 * @author  Jianping Meng
 * @details Define kernel functions for calculating the equilibrium functions,
 * the relaxation time, body force term, and macroscopic variables.
 */

#ifndef MODEL_KERNEL_INC
#define MODEL_KERNEL_INC
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "model.h"
#include "flowfield_host_device.h"
#include "model_host_device.h"

/*!
 * We assume that the layout of MacroVars is rho, u, v, w, T, ...
 * In the macroVars, all variables are conserved, i.e., with density.
 * @todo how to deal with overflow in a kernel function? in particular, GPU
 */
#ifdef OPS_2D
void KerCalcDensity(ACC<Real>& Rho, const ACC<Real>& f,
                    const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            rho += f(xiIdx, 0, 0);
        }
#ifdef CPU
        if (isnan(rho) || rho <= 0 || isinf(rho)) {
            ops_printf(
                "Error! Density %f becomes invalid！Something "
                "wrong...",
                rho);
            assert(!(isnan(rho) || rho <= 0 || isinf(rho)));
        }
#endif
        Rho(0, 0) = rho;
    }
#endif  // OPS_2D
}

void KerCalcU(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
              const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! "
                "Maybe something wrong...\n",
                u);
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0) = u;
    }
#endif  // OPS_2D
}

void KerCalcV(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
              const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! "
                "Maybe something wrong...\n",
                v);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0) = v;
    }
#endif  // OPS_2D
}

void KerCalcUForce(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            u += ((*dt) * acceleration(0, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f\n",
                u, x, y);
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0) = u;
    }
#endif  // OPS_2D
}

void KerCalcVForce(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            v += ((*dt) * acceleration(1, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f\n",
                v, x, y);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0) = v;
    }
#endif  // OPS_2D
}

/*!
 * If a Newton-Cotes quadrature is used, it can be converted to the way
 * similar to the Gauss-Hermite quadrature *
 */

void KerInitialiseBGK2nd(ACC<Real>& f, const ACC<int>& nodeType,
                         const ACC<Real>& Rho, const ACC<Real>& U,
                         const ACC<Real>& V, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0) = CalcBGKFeq(xiIdx, rho, u, v, T, polyOrder);
#ifdef CPU
            const Real res{f(xiIdx, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCollideBGKIsothermal(ACC<Real>& fStage, const ACC<Real>& f,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Rho,
                             const ACC<Real>& U, const ACC<Real>& V,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
            } else {
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, feq, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCollideBGKThermal(ACC<Real>& fStage, const ACC<Real>& f,
                          const ACC<int>& nodeType, const ACC<Real>& Rho,
                          const ACC<Real>& U, const ACC<Real>& V,
                          const ACC<Real>& Temperature, const Real* tauRef,
                          const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real T{Temperature(0, 0)};
        const int polyOrder{4};
        Real tau = (*tauRef) / (rho * sqrt(T));
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0) =
                    f(xiIndex, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
            } else {
                fStage(xiIndex, 0, 0) =
                    f(xiIndex, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIndex);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCalcBodyForce1ST(ACC<Real>& fStage, const ACC<Real>& acceration,
                         const ACC<Real>& Rho, const ACC<int>& nodeType,
                         const int* lattIdx) {
#ifdef OPS_2D

    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real rho{Rho(0, 0)};
        Real g[]{acceration(0, 0, 0), acceration(1, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcBodyForce(xiIndex, rho, g)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            fStage(xiIndex, 0, 0) = bodyForce;
        }
    }
#endif  // OPS_2D
}

void KerCalcBodyForceNone(ACC<Real>& fStage, const ACC<Real>& acceration,
                          const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fStage(xiIndex, 0, 0) = 0;
        }
    }
#endif  // OPS_2D
}
#endif // OPS_2D outter

#ifdef OPS_3D
void KerInitialiseBGK2nd3D(ACC<Real>& f, const ACC<int>& nodeType,
                           const ACC<Real>& Rho, const ACC<Real>& U,
                           const ACC<Real>& V, const ACC<Real>& W,
                           const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0, 0) = CalcBGKFeq(xiIdx, rho, u, v, w, T, polyOrder);
#ifdef CPU
            const Real res{f(xiIdx, 0, 0, 0)};
            if (std::isnan(res) || res <= 0 || std::isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(std::isnan(res) || res <= 0 || std::isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerCollideBGKIsothermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                               const ACC<Real>& coordinates,
                               const ACC<int>& nodeType, const ACC<Real>& Rho,
                               const ACC<Real>& U, const ACC<Real>& V,
                               const ACC<Real>& W, const Real* tauRef,
                               const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
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
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0, 0);
            } else {
                fStage(xiIndex, 0, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0, 0)};
            if (std::isnan(res) || res <= 0 || std::isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e w=%e at "
                    "x=%e y=%e z=%e\n",
                    res, xiIndex, feq, rho, u, v, w, coordinates(0, 0, 0, 0),
                    coordinates(1, 0, 0, 0), coordinates(2, 0, 0, 0));
                assert(!(std::isnan(res) || res <= 0 || std::isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerCollideBGKThermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                            const ACC<int>& nodeType, const ACC<Real>& Rho,
                            const ACC<Real>& U, const ACC<Real>& V,
                            const ACC<Real>& W, const ACC<Real>& Temperature,
                            const Real* tauRef, const Real* dt,
                            const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real T{Temperature(0, 0, 0)};
        const int polyOrder{4};
        Real tau = (*tauRef) / (rho * sqrt(T));
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0, 0) =
                    f(xiIndex, 0, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0, 0);
            } else {
                fStage(xiIndex, 0, 0, 0) =
                    f(xiIndex, 0, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0, 0)};
            if (std::isnan(res) || res <= 0 || std::isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIndex);
                assert(!(std::isnan(res) || res <= 0 || std::isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerCalcBodyForce1ST3D(ACC<Real>& fStage, const ACC<Real>& acceration,
                           const ACC<Real>& Rho, const ACC<int>& nodeType,
                           const int* lattIdx) {
#ifdef OPS_3D

    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real rho{Rho(0, 0, 0)};
        Real g[]{acceration(0, 0, 0, 0), acceration(1, 0, 0, 0),
                 acceration(2, 0, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcBodyForce(xiIndex, rho, g)};
#ifdef CPU
            if (std::isnan(bodyForce) || std::isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(std::isnan(bodyForce) || std::isinf(bodyForce)));
            }
#endif
            fStage(xiIndex, 0, 0, 0) = bodyForce;
        }
    }
#endif  // OPS_3D
}

void KerCalcBodyForceNone3D(ACC<Real>& fStage, const ACC<Real>& acceration,
                            const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fStage(xiIndex, 0, 0, 0) = 0;
        }
    }
#endif  // OPS_3D
}

void KerCalcDensity3D(ACC<Real>& Rho, const ACC<Real>& f,
                      const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            rho += f(xiIdx, 0, 0, 0);
        }
#ifdef CPU
        if (std::isnan(rho) || rho <= 0 || std::isinf(rho)) {
            ops_printf(
                "Error! Density %f becomes invalid！Something "
                "wrong...",
                rho);
            assert(!(std::isnan(rho) || rho <= 0 || std::isinf(rho)));
        }
#endif
        Rho(0, 0, 0) = rho;
    }
#endif // OPS_3D
}

void KerCalcU3D(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0, 0);
        }
        u /= Rho(0, 0, 0);
#ifdef CPU
        if (std::isnan(u) || std::isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! "
                "Maybe something wrong...\n",
                u);
            assert(!(std::isnan(u) || std::isinf(u)));
        }
#endif
        U(0, 0, 0) = u;
    }
#endif // OPS_3D
}

void KerCalcV3D(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0, 0);
        }
        v /= Rho(0, 0, 0);
#ifdef CPU
        if (std::isnan(v) || std::isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! "
                "Maybe something wrong...\n",
                v);
            assert(!(std::isnan(v) || std::isinf(v)));
        }
#endif
        V(0, 0, 0) = v;
    }
#endif // OPS_3D
}

void KerCalcW3D(ACC<Real>& W, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real w{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            w += CS * XI[xiIdx * LATTDIM + 2] * f(xiIdx, 0, 0, 0);
        }
        w /= Rho(0, 0, 0);
#ifdef CPU
        if (std::isnan(w) || std::isinf(w)) {
            ops_printf(
                "Error! Velocity W=%f becomes invalid! "
                "Maybe something wrong...\n",
                w);
            assert(!(std::isnan(w) || std::isinf(w)));
        }
#endif
        W(0, 0, 0) = w;
    }
#endif // OPS_3D
}

void KerCalcUForce3D(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0, 0);
        }
        u /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            u += ((*dt) * acceleration(0, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (std::isnan(u) || std::isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f\n",
                u, x, y, z);
            assert(!(std::isnan(u) ||std::isinf(u)));
        }
#endif
        U(0, 0, 0) = u;
    }
#endif // OPS_3D
}

void KerCalcVForce3D(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0, 0);
        }
        v /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            v += ((*dt) * acceleration(1, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (std::isnan(v) || std::isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f\n",
                v, x, y, z);
            assert(!(std::isnan(v) || std::isinf(v)));
        }
#endif
        V(0, 0, 0) = v;
    }
#endif // OPS_3D
}

void KerCalcWForce3D(ACC<Real>& W, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real w{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            w += CS * XI[xiIdx * LATTDIM + 2] * f(xiIdx, 0, 0, 0);
        }
        w /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            w += ((*dt) * acceleration(2, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (std::isnan(w) || std::isinf(w)) {
            ops_printf(
                "Error! Velocity W=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f\n",
                w, x, y, z);
            assert(!(std::isnan(w) || std::isinf(w)));
        }
#endif
        W(0, 0, 0) = w;
    }
#endif //OPS_3D
}
#endif //OPS_3D outter

#endif //MODEL_KERNEL_INC