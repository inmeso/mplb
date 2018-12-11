// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

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
void KerSetEmbededBodyBoundary(int* surfaceBoundary,
                               const int* geometryProperty, int* nodeType);
/*!
 * set geometry property for boundary surface points 
 */
void KerSetEmbededBodyGeometry(const int* nodeType, int* geometryProperty);
/*!
 * set the initial conditions
 */
void KerSetInitialMacroVars(const Real* coordinates, const int* idx,
                            Real* macroVars);

/*!
 * mark solid points enclosed in a circle
 */
void KerSetEmbededCircle(Real* diameter, Real* centerPos,
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
