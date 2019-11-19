
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

/*! @brief   Head files for importing geometry from HDF5 file
  * @author  Jianping Meng
  * @details Declaring kernel functions related to create computing domain
  */
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <string>
  #include <vector>
  #include "flowfield.h"
  #include "scheme.h"
#ifdef OPS_2D
void KerSetCoordinates(const Real* coordX, const Real* coordY, const int* idx,
                       Real* coordinates);

/*!
 * set the boundary condition type to boundary surface points
 */
void KerSetEmbeddedBodyBoundary(int* surfaceBoundary,
                               const int* geometryProperty, int* nodeType);
/*!
 * set geometry property for boundary surface points
 */
void KerSetEmbeddedBodyGeometry(const int* nodeType, int* geometryProperty);
/*!
 * set the initial conditions
 */
void KerSetInitialMacroVars(const Real* coordinates, const int* idx,
                            Real* macroVars);

/*!
 * mark solid points enclosed in a circle
 */
void KerSetEmbeddedCircle(Real* diameter, Real* centerPos,
                         const Real* coordinates, int* nodeType,
                         int* geometryProperty);
/*!
 * wipe off all dangling solid points
 */
void KerSweep(const int* geometryProperty, int* nodeType);
/*!
 * Sync the geometryProperty with the nodeType
 * Use this function after KerSweep
 */
void KerSyncGeometryProperty(const int* nodeType, int* geometryProperty);
#endif //OPS_2D
#ifdef OPS_3D
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx, Real* coordinates);
#endif //OPS_3D
