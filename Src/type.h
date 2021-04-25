
/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]s
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

/*! @brief   Define constants and enumeration types.
 * @author  Jianping Meng
 * @details Define important constants and enumeration types including
 * boundary type, geometry types etc.
 * cycle
 */

#ifndef TYPE_H
#define TYPE_H
#include <limits>
#define DP
#ifdef DP
typedef double Real;
//const char* RealC = "double";
#else
typedef float Real;
//const char* RealC = "float";
#endif
const Real PI{3.1415926535897932384626433832795};
const Real EPS{std::numeric_limits<Real>::epsilon()};
const Real BOLTZ{1.3806488e-23};
const int xaxis = 0;
const int yaxis = 1;
#ifdef OPS_3D
const int zaxis = 2;
#endif
// ZERO is the zero constant with the desired precision, i.e., float or double
const Real ZERO{(Real)((int)0)};
typedef std::size_t SizeType;
#include "assert.h"
#include <vector>
#include <cmath>

// It looks that OPS always fills the uninitialised storage with 0 so
// we try to avoid 0 value for these types
// Use this type to describe different node type in terms of boundary types
// or a fluid node
/*
 * Defining the position of a grid point relative to a polygon or polyhedron
 */
enum PointPosition {
    IsVertex = 1,
    StrictlyInterior = 2,
    StrictlyExterior = 3,
    RelativelyInteriorToEdge = 4,
    RelativelyInteriorToFace = 5
};



enum SolidBodyType { SolidBody_circle = 0, SolidBody_ellipse = 1 };

enum SpaceSchemeType {
    sstupwind2nd = 10,
    sstcentral2nd = 11,
    sstcentral4th = 12,
    sstcentral6th = 13,
    sstupwind3rd = 14,
    sstupwind4th = 15,
    sstupwind5h = 16,
    sstupwind6th = 17
};

typedef enum { tstrk4th = 10, tsteuler1st = 11 } TimeSchemeType;

typedef enum {
    qtdensity = 0,
    qtvelocity = 1,
    qttemperature = 2,
    qtpressure = 3,
    qtheatflux = 4,
    qtaccerlation = 5,
    qttime = 6,
    qtlength = 7,

} QuantityType;

typedef enum {
    fthardspheregas = 0,
    ftmaxwellgas = 1,
} FluidType;


typedef enum {
    bpanormal = 0,
    bpacorner = 1,
} BoundaryPointAttribute;

typedef enum {
    sbpif = 0,  // f first index of i
    sbpie = 1,  // e the end of i
    sbpjf = 2,
    sbpje = 3,
    sbpkf = 4,
    sbpke = 5,
} StructuredBoundaryPosition;

enum SchemeType {
    Scheme_E1st2nd = 1,
    Scheme_I1st2nd = -1,
    Scheme_StreamCollision = 10,
} ;

inline bool EssentiallyEqual(const Real* a, const Real* b, const Real epsilon) {
    return fabs(*a - *b) <=
           ((fabs(*a) > fabs(*b) ? fabs(*b) : fabs(*a)) * epsilon);
}

inline bool DefinitelyGreaterThan(const Real* a, const Real* b,
                                  const Real epsilon) {
    return (*a - *b) > ((fabs(*a) < fabs(*b) ? fabs(*b) : fabs(*a)) * epsilon);
}

inline bool DefinitelyLessThan(const Real* a, const Real* b,
                               const Real epsilon) {
    return (*b - *a) > ((fabs(*a) < fabs(*b) ? fabs(*b) : fabs(*a)) * epsilon);
}

template <typename T>
inline void FreeArrayMemory(T* pointerToArray) {
    if (pointerToArray != nullptr) {
        delete[] pointerToArray;
    }
}
#endif  // TYPE_H
