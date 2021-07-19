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

#ifndef CAVITY2D_KERNEL_INC
#define CAVITY2D_KERNEL_INC

void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            const ACC<Real>& coordinates, const int* idx) {
    rho(0, 0) = 1;
    u(0, 0) = 0;
    v(0, 0) = 0;
}

//ADCase

void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
    Real Sig02{200.0};
    Real DiffCFact{3.0};
    Real dt = {*Dt};
    Real time = {*TIME};
    Real DiffC{DiffCFact/dt};
    Real x0[2] = {200.0,200.0};
    Real Vel[2] = {U(0,0),V(0,0)};
    Real SigD2{2*DiffC*time*0};
    Real Pref{Sig02/(Sig02+SigD2)};
    Real xiIndexVar1 = {xiIndex(0,0,0)};
    Real xiIndexVar2 = {xiIndex(1,0,0)};
    Real numeru{(24.0*xiIndexVar1-x0[0]-Vel[0]*time*0)};
    Real numerv{(24.0*xiIndexVar2-x0[1]-Vel[1]*time*0)};
    Real numer2{numeru*numeru+numerv*numerv};
    Real Gauss{exp(-numer2/(2*(Sig02+SigD2)))};

    concentration(0,0) = Pref * Gauss;
}

#endif // CAVITY2D_KERNEL_INC
