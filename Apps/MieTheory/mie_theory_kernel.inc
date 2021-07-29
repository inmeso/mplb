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

#ifndef MIE_THEORY_INC
#define MIE_THEORY_INC

static besselFunctions* bessels;
double blkExtent[3];

void KerSetSphericalCoord(ACC<double>& coordinates, const double* blkExtent,
                            const int* idx) {
    // get the index i, j, k for the current grid point
    const int i{idx[0]};
    const int j{idx[1]};
    const int k{idx[2]};
    // calculate the coorindate based on index i, j,k
    std::cout << blkExtent[0] << blkExtent[2] << blkExtent[2] << std::endl;
    coordinates(0, 0, 0, 0) = i;
    coordinates(1, 0, 0, 0) = j;
    coordinates(2, 0, 0, 0) = k;
}

void KerCalculateMieSolution(ACC<double>& E, ACC<double>& H,
                             const ACC<double>& coordinates) {
    // get the index i, j, k for the current grid point
    const double x{coordinates(0, 0, 0, 0)};
    const double y{coordinates(1, 0, 0, 0)};
    const double z{coordinates(2, 0, 0, 0)};
    // calculate the coorindate based on index i, j,k

    E(0, 0, 0, 0) = x;
    E(1, 0, 0, 0) = y;
    E(2, 0, 0, 0) = z;
    H(0, 0, 0, 0) = x;
    H(1, 0, 0, 0) = y;
    H(2, 0, 0, 0) = z;
}

#endif  // MIE_THEORY_INC
