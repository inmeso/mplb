#ifndef Hilemms_H
#define Hilemms_H

#include <string.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include "boundary.h"
#include "evolution.h"
#include "evolution3d.h"
#include "flowfield.h"
#include "model.h"
#include "ops_seq.h"
#include "scheme.h"
#include "type.h"
// blockNum: total number if blocks.
// blockSize: array of integers specifying the block blocksize.
// meshSize: The size of mesh i.e. dx (At present dx = dy = dz).
// startPos: Starting position of each block.
void DefineProblemDomain(const int blockNum, const std::vector<int> blockSize,
                         const Real meshSize, const std::vector<Real> startPos);

// Iterator for transient simulations.
void Iterate(const int steps, const int checkPointPeriod);

// Iterator for steady simulations.
void Iterate(const Real convergenceCriteria, const int checkPointPeriod);

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
void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         BoundaryType boundaryType,
                         const std::vector<VariableTypes>& macroVarTypes,
                         const std::vector<Real>& macroVarValues);

// blockIndex: Block Index
// blockStartPos: Starting position of the block in x, y, z directions to find
// coordinates. meshSize: Size of mesh (Assuming a constant in all 3 directions
// for stream collision scheme).
void CalBlockCoordinates(const int blockIndex, Real* blockStartPos,
                         Real meshsize);


void KerSetCoordinates(const Real* coordX, const Real* coordY, const int* idx,
                       Real* coordinates);
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx, Real* coordinates);

// Interface for Setting the initial values.
// coordinates: Coordinate array of nodes.
// idx: An array used by OPS which gives the index of the current grid point.
void KerSetInitialMacroVarsHilemms(const Real* coordinates, const int* idx,
                                   Real* macroVars, Real* macroVarsInitVal,
                                   const int* componentId);
// Kernel which will call a user-defined function for inital conditions
void KerSetInitialMacroVars(Real* macroVars, const Real* coordinates,
                            const int* idx);
//void AssignCoordinates(int blockIndex, Real* coordinates[SPACEDIM]);
void AssignCoordinates(int blockIndex, Real** coordinates);
// blockIndex: Block id.
// Coordinates: Array to store coordinates of various points inside the domain.
// This function was defined in setup_comput_domain and has been declared here
// as we are not using the preprocessor code separately.

// blockIndex: Block Id.
// componentId: Id of the component.
// initialMacroValues: Initial values of the macroscopic variables for a given
// component in a particular block.
void DefineInitialCondition(int blockIndex, int componentId,
                      std::vector<Real> initialMacroValues);

//User-defined function for initialising macroscopic variables
void InitialiseNodeMacroVars(Real* nodeMacroVars, const Real* nodeCoordinates);
//Defining the initial conditions by using user-defined functions
void DefineInitialCondition();

void SetupGeomPropAndNodeType(int blockIndex, BoundaryType* boundType);
void SetupGeomPropAndNodeType(int blockIndex, BoundaryType* boundType,
                              BoundaryType* edgeType, BoundaryType* cornerType);
void SetupDomainGeometryProperty(int blockIndex);

void DefineHaloNumber(int Halo_Number, int Halo_Depth, int Scheme_Halo_points,
                      int Num_Bound_Halo_Points);
// Functions to check for inclusion.

// A wrapper Function which implements all the boundary conditions.
void ImplementBoundaryConditions();


// type: Circle/Sphere, Ellipse/Ellipsoid, superquadrics, ...
// centerPos: the position vector of the center point.
// controlParas: control parameters, e.g. radius for Circle/Sphere, ...
void EmbeddedBody(SolidBodyType type, int blockIndex,
                  std::vector<Real> centerPos, std::vector<Real> controlParas);


/**********************************************************/
/* Functions for embedded body.                           */
/**********************************************************/

void KerSetEmbeddedBodyBoundary(int* surfaceBoundary,
                               const int* geometryProperty, int* nodeType);

void KerSetEmbeddedCircle(Real* diameter, Real* centerPos,
                         const Real* coordinates, int* nodeType,
                         int* geometryProperty);

void KerSetEmbeddedEllipse(Real* semiMajorAxes, Real* semiMinorAxis,
                           Real* centerPos, const Real* coordinates,
                           int* nodeType, int* geometryProperty);

void KerSweep(const int* geometryProperty, int* nodeType);

void KerSyncGeometryProperty(const int* nodeType, int* geometryProperty);

void KerSetEmbeddedBodyGeometry(const int* nodeType, int* geometryProperty);

void HandleImmersedSolid();

#endif  // Hilemms_H