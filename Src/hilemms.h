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

#ifndef Hilemms_H
#define Hilemms_H

#include <string.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include "boundary.h"
#include "flowfield.h"
#include "model.h"
#include "scheme.h"
#include "type.h"

// Add 2D polygon.
// vertexNum: total number of vertexes.
// vertexCoords: Coordinates of each vertex.
void AddEmbeddedBody(int vertexNum, Real* vertexCoords);

// blockIndex: block Index
// compoId: component ID whose BC we want to set.
// surface: which surface to set.
// type: boundary condition type.
// MacroVarsComp: which all macrovars are to be used in specifying the BC.
// valueMacroVarsComp: specified value for the boundary condition for the macro
// vars which are defined using MacroVarsComp.

void SetupGeomPropAndNodeType(int blockIndex, BoundaryScheme* boundType);
void SetupGeomPropAndNodeType(int blockIndex, BoundaryScheme* boundType,
                              BoundaryScheme* edgeType,
                              BoundaryScheme* cornerType);

// type: Circle/Sphere, Ellipse/Ellipsoid, superquadrics, ...
// centerPos: the position vector of the center point.
// controlParas: control parameters, e.g. radius for Circle/Sphere, ...
void AddEmbeddedBody(SolidBodyType type, int blockIndex,
                     std::vector<Real> centerPos,
                     std::vector<Real> controlParas);

/**********************************************************/
/* Functions for embedded body.                           */
/**********************************************************/

void KerSetEmbeddedBodyBoundary(ACC<int>& nodeType,
                                const ACC<int>& geometryProperty,
                                int* surfaceBoundary);

void KerSetEmbeddedCircle(ACC<int>& nodeType, ACC<int>& geometryProperty,
                          const Real* coordinates, Real* diameter,
                          Real* centerPos);

void KerSetEmbeddedEllipse(ACC<int>& nodeType, ACC<int>& geometryProperty,
                           const ACC<Real>& coordinates, Real* semiMajorAxes,
                           Real* semiMinorAxis, Real* centerPos);

void KerSweep(ACC<int> nodeType, const ACC<int>& geometryProperty);

void KerSetEmbeddedBodyGeometry(ACC<int>& geometryProperty,
                                const ACC<int>& nodeType);

void KerSyncGeometryProperty(ACC<int>& geometryProperty,
                             const ACC<int>& nodeType);

void WipeSolidPtsBasedNeigbours();
void UpdateGeometryAfterWiping();
void MarkSurfacePoints();
void SetBoundaryTypeofImmersedBody();

#endif  // Hilemms_H