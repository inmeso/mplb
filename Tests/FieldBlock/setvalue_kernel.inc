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

#ifndef SETVALUE_KERNEL_INC
#define SETVALUE_KERNEL_INC
#ifdef OPS_3D
void KerSetValues(ACC<double>& geo) {
    geo(0, 0, 0, 0) = 1;
    geo(1, 0, 0, 0) = 0;
    geo(2, 0, 0, 0) = 0;
    geo(3, 0, 0, 0) = 0;
}
void KerSetValuesSingle(ACC<int>& macroVars) { macroVars(0, 0, 0) = 1; }

#endif  // OPS_3D
#ifdef OPS_2D
void KerSetValues(ACC<double>& geo) {
    geo(0, 0, 0) = 1;
    geo(1, 0, 0) = 0;
    geo(2, 0, 0) = 0;
    geo(3, 0, 0) = 0;
}

void KerSetValuesSingle(ACC<int>& macroVars) { macroVars(0, 0) = 1; }
#endif  // OPS_2D
#endif  // SETVALUE_KERNEL_INC
