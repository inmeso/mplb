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
#ifndef FLOWFIELD_KERNEL_INC
#define FLOWFIELD_KERNEL_INC
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"

void KerSetIntField(const int* value, ACC<int>& var) {
#ifdef OPS_2D
    var(0, 0) = (*value);
#endif
#ifdef OPS_3D
    var(0, 0, 0) = (*value);
#endif
}

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

void KerCopyf(ACC<Real>& dest, const ACC<Real>& src) {
    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
#ifdef OPS_2D
        dest(xiIndex, 0, 0) = src(xiIndex, 0, 0);
#endif
#ifdef OPS_3D
        dest(xiIndex, 0, 0, 0) = src(xiIndex, 0, 0, 0);
#endif
    }
}

void KerCopyIntField(const ACC<int>& src, ACC<int>& dest) {
#ifdef OPS_2D
    dest(0, 0) = src(0, 0);
#endif
#ifdef OPS_3D
    dest(0, 0, 0) = src(0, 0, 0);
#endif
}

void KerCopyMacroVars(const ACC<Real>& src, ACC<Real>& dest) {
#ifdef OPS_2D
    dest(0, 0) = src(0, 0);
#endif
#ifdef OPS_3D
    dest(0, 0, 0) = src(0, 0, 0);
#endif
}

void KerCopyCoordinateXYZ(const ACC<Real>& src, ACC<Real>& dest) {
#ifdef OPS_2D
    for (int idx = 0; idx < 2; idx++) {
        dest(idx, 0, 0) = src(idx, 0, 0);
    }
#endif
#ifdef OPS_3D
    for (int idx = 0; idx < 3; idx++) {
        dest(idx, 0, 0, 0) = src(idx, 0, 0, 0);
    }
#endif
}

void KerCopyDispf(const ACC<Real>& src, ACC<Real>& dest, const int* disp) {
    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
#ifdef OPS_2D
        dest(xiIndex, disp[0], disp[1]) = src(xiIndex, 0, 0);
#endif
#ifdef OPS_3D
        dest(xiIndex, disp[0], disp[1], disp[2]) = src(xiIndex, 0, 0, 0);
#endif
    }
}

void KerNormaliseF(const Real* ratio, ACC<Real>& f) {
    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
#ifdef OPS_2D
        f(xiIndex, 0, 0) /= (*ratio);
#endif
#ifdef OPS_3D
        f(xiIndex, 0, 0, 0) /= (*ratio);
#endif
    }
}

void KerCalcMacroVarSquareofDifference(const ACC<Real>& macroVars,
                                       const ACC<Real>& macroVarsCopy,
                                       double* sumSquareDiff) {
#ifdef OPS_2D
    *sumSquareDiff =
        *sumSquareDiff + (macroVars(0, 0) - macroVarsCopy(0, 0)) *
                             (macroVars(0, 0) - macroVarsCopy(0, 0));
#endif
#ifdef OPS_3D
    *sumSquareDiff =
        *sumSquareDiff + (macroVars(0, 0, 0) - macroVarsCopy(0, 0, 0)) *
                             (macroVars(0, 0, 0) - macroVarsCopy(0, 0, 0));
#endif
}

void KerCalcMacroVarSquare(const ACC<Real>& macroVars, double* sumSquare) {
#ifdef OPS_2D
    *sumSquare = *sumSquare + (macroVars(0, 0)) * (macroVars(0, 0));
#endif
#ifdef OPS_3D
    *sumSquare = *sumSquare + (macroVars(0, 0, 0)) * (macroVars(0, 0, 0));
#endif
}

void KerSetfFixValue(const Real* value, ACC<Real>& f) {
    for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
#ifdef OPS_2D
        f(xiIndex, 0, 0) = (*value);
#endif
#ifdef OPS_3D
        f(xiIndex, 0, 0, 0) = (*value);
#endif
    }
}

#endif //FLOWFIELD_KERNEL_INC