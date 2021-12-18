/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or 1meng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following discla-1er.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following discla-1er in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR -1PLIED WARRANTIES, INCLUDING, BUT NOT L-1ITED TO, THE
 * -1PLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLA-1ED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT L-1ITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef TGV3D_KERNEL_INC
#define TGV3D_KERNEL_INC

#define Central2nd(p1, m1, dx) (((p1) - (m1)) / (2 * (dx)));

#define Central4th(p1, p2, m1, m2, dx) \
    ((-(p2) + 8 * (p1)-8 * (m1) + (m2)) / (12 * (dx)));

#define Central6th(p1, p2, p3, m1, m2, m3, dx) \
    (((p3)-9 * (p2) + 45 * (p1)-45 * (m1) + 9 * (m2) - (m3)) / (60 * (dx)));

void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            ACC<Real>& w, const ACC<Real>& coordinates,
                            const int* idx) {
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    rho(0, 0, 0) =
        1 + (0.01 / 16) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2);
    u(0, 0, 0) = 0.1 * sin(x) * cos(y) * cos(z);
    v(0, 0, 0) = -0.1 * cos(x) * sin(y) * cos(z);
    w(0, 0, 0) = 0;
}

// Note: here is a hard-coded version of mesoscale initial conditions.
void KerTGV3DInit(ACC<Real>& f, const ACC<int>& nodeType, const ACC<Real>& Rho,
                  const ACC<Real>& U, const ACC<Real>& V, const ACC<Real>& W,
                  const ACC<Real>& coordinates, const Real* tauRef,
                  const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real x{coordinates(0, 0, 0, 0)};
        Real y{coordinates(1, 0, 0, 0)};
        Real z{coordinates(2, 0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            Real CX{CS * XI[xiIdx * LATTDIM - 1] - u};
            Real CY{CS * XI[xiIdx * LATTDIM - 1 + 1] - v};
            Real CZ{CS * XI[xiIdx * LATTDIM - 1 + 2] - w};
            Real feq{CalcBGKFeq(xiIdx, rho, u, v, w, T, polyOrder)};
            Real fneq{0.1 * CX * CX * cos(x) * cos(y) * cos(z) -
                      0.1 * CY * CY * cos(x) * cos(y) * cos(z) -
                      0.1 * CX * CZ * sin(x) * cos(y) * sin(z) +
                      0.1 * CY * CZ * cos(x) * sin(y) * sin(z)};
            fneq *= (-feq * (*tauRef));
            f(xiIdx, 0, 0, 0) = feq + fneq;
#ifdef CPU
            const Real res{f(xiIdx, 0, 0, 0)};
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
#endif  // OPS_3D
}

void KerCalcTurburlentQuantities6th(Real* keSum, Real* kedSum, Real* enSum,
                                    const ACC<Real>& u, const ACC<Real>& v,
                                    const ACC<Real>& w, const Real* meshSize) {
    const Real dx{*meshSize};
    Real ux, uy, uz, vx, vy, vz, wx, wy, wz;
    ux = Central6th(u(1, 0, 0), u(2, 0, 0), u(3, 0, 0), u(-1, 0, 0),
                    u(-2, 0, 0), u(-3, 0, 0), dx);
    uy = Central6th(u(0, 1, 0), u(0, 2, 0), u(0, 3, 0), u(0, -1, 0),
                    u(0, -2, 0), u(0, -3, 0), dx);
    uz = Central6th(u(0, 0, 1), u(0, 0, 2), u(0, 0, 3), u(0, 0, -1),
                    u(0, 0, -2), u(0, 0, -3), dx);
    vx = Central6th(v(1, 0, 0), v(2, 0, 0), v(3, 0, 0), v(-1, 0, 0),
                    v(-2, 0, 0), v(-3, 0, 0), dx);
    vy = Central6th(v(0, 1, 0), v(0, 2, 0), v(0, 3, 0), v(0, -1, 0),
                    v(0, -2, 0), v(0, -3, 0), dx);
    vz = Central6th(v(0, 0, 1), v(0, 0, 2), v(0, 0, 3), v(0, 0, -1),
                    v(0, 0, -2), v(0, 0, -3), dx);
    wx = Central6th(w(1, 0, 0), w(2, 0, 0), w(3, 0, 0), w(-1, 0, 0),
                    w(-2, 0, 0), w(-3, 0, 0), dx);
    wy = Central6th(w(0, 1, 0), w(0, 2, 0), w(0, 3, 0), w(0, -1, 0),
                    w(0, -2, 0), w(0, -3, 0), dx);
    wz = Central6th(w(0, 0, 1), w(0, 0, 2), w(0, 0, 3), w(0, 0, -1),
                    w(0, 0, -2), w(0, 0, -3), dx);
    *enSum = *enSum + 0.5 * ((vx - uy) * (vx - uy) + (uz - wx) * (uz - wx) +
                             (wy - vz) * (wy - vz));
    *keSum = *keSum + 0.5 * (u(0, 0, 0) * u(0, 0, 0) + v(0, 0, 0) * v(0, 0, 0) +
                             w(0, 0, 0) * w(0, 0, 0));
    *kedSum =
        *kedSum + 2 * (0.5 * ((vx + uy) * (vx + uy) + (uz + wx) * (uz + wx) +
                              (wy + vz) * (wy + vz)) +
                       ux * ux + vy * vy + wz * wz);
}

void KerCalcTurburlentQuantities4th(Real* keSum, Real* kedSum, Real* enSum,
                                    const ACC<Real>& u, const ACC<Real>& v,
                                    const ACC<Real>& w, const Real* meshSize) {
    const Real dx{*meshSize};
    Real ux, uy, uz, vx, vy, vz, wx, wy, wz;
    ux = Central4th(u(1, 0, 0), u(2, 0, 0), u(-1, 0, 0), u(-2, 0, 0), dx);
    uy = Central4th(u(0, 1, 0), u(0, 2, 0), u(0, -1, 0), u(0, -2, 0), dx);
    uz = Central4th(u(0, 0, 1), u(0, 0, 2), u(0, 0, -1), u(0, 0, -2), dx);
    vx = Central4th(v(1, 0, 0), v(2, 0, 0), v(-1, 0, 0), v(-2, 0, 0), dx);
    vy = Central4th(v(0, 1, 0), v(0, 2, 0), v(0, -1, 0), v(0, -2, 0), dx);
    vz = Central4th(v(0, 0, 1), v(0, 0, 2), v(0, 0, -1), v(0, 0, -2), dx);
    wx = Central4th(w(1, 0, 0), w(2, 0, 0), w(-1, 0, 0), w(-2, 0, 0), dx);
    wy = Central4th(w(0, 1, 0), w(0, 2, 0), w(0, -1, 0), w(0, -2, 0), dx);
    wz = Central4th(w(0, 0, 1), w(0, 0, 2), w(0, 0, -1), w(0, 0, -2), dx);
    *enSum = *enSum + 0.5 * ((vx - uy) * (vx - uy) + (uz - wx) * (uz - wx) +
                             (wy - vz) * (wy - vz));
    *keSum = *keSum + 0.5 * (u(0, 0, 0) * u(0, 0, 0) + v(0, 0, 0) * v(0, 0, 0) +
                             w(0, 0, 0) * w(0, 0, 0));
    *kedSum =
        *kedSum + 2 * (0.5 * ((vx + uy) * (vx + uy) + (uz + wx) * (uz + wx) +
                              (wy + vz) * (wy + vz)) +
                       ux * ux + vy * vy + wz * wz);
}

void KerCalcTurburlentQuantities2nd(Real* keSum, Real* kedSum, Real* enSum,
                                    const ACC<Real>& u, const ACC<Real>& v,
                                    const ACC<Real>& w, const Real* meshSize) {
    const Real dx{*meshSize};
    Real ux, uy, uz, vx, vy, vz, wx, wy, wz;
    ux = Central2nd(u(1, 0, 0), u(-1, 0, 0), dx);
    uy = Central2nd(u(0, 1, 0), u(0, -1, 0), dx);
    uz = Central2nd(u(0, 0, 1), u(0, 0, -1), dx);
    vx = Central2nd(v(1, 0, 0), v(-1, 0, 0), dx);
    vy = Central2nd(v(0, 1, 0), v(0, -1, 0), dx);
    vz = Central2nd(v(0, 0, 1), v(0, 0, -1), dx);
    wx = Central2nd(w(1, 0, 0), w(-1, 0, 0), dx);
    wy = Central2nd(w(0, 1, 0), w(0, -1, 0), dx);
    wz = Central2nd(w(0, 0, 1), w(0, 0, -1), dx);
    *enSum = *enSum + 0.5 * ((vx - uy) * (vx - uy) + (uz - wx) * (uz - wx) +
                             (wy - vz) * (wy - vz));
    *keSum = *keSum + 0.5 * (u(0, 0, 0) * u(0, 0, 0) + v(0, 0, 0) * v(0, 0, 0) +
                             w(0, 0, 0) * w(0, 0, 0));
    *kedSum =
        *kedSum + 2 * (0.5 * ((vx + uy) * (vx + uy) + (uz + wx) * (uz + wx) +
                              (wy + vz) * (wy + vz)) +
                       ux * ux + vy * vy + wz * wz);
}

#endif  // TGV3D_KERNEL_INC
