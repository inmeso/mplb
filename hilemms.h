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

extern int SPACEDIM;

void DefineCase(std::string caseName, const int spaceDim);
// caseName: case name
// spaceDim: 2D or 3D application

void DefineProblemDomain(const int blockNum, const std::vector<int> blockSize,
                         const Real meshSize, const std::vector<Real> startPos);
// blockNum: total number if blocks.
// blockSize: array of integers specifying the block blocksize.
// meshSize: The size of mesh i.e. dx (At present dx = dy = dz).
// startPos: Starting position of each block.

/*
void DefineForceTerm(std::vector<ForceType> types, std::vector<int> compoId);
// types: which kind of force function to use.
// compoID: which component to act on.
*/

void Iterate(SchemeType scheme, const int steps, const int checkPointPeriod);
// scheme: which schemes to use for implementing such as finite difference
// scheme. for steady state simulations.

void Iterate(SchemeType scheme, const Real convergenceCriteria,
             const int checkPointPeriod);
// Same as iterate but for transient simulations.

void AddEmbededBody(int vertexNum, Real* vertexCoords);
// Add 2D polygon.
// vertexNum: total number of vertexes.
// vertexCoords: Coordinates of each vertex.

void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface surface, BoundaryType type,
                         std::vector<VariableTypes> MacorVarsComp,
                         std::vector<Real> valuesMacroVarsComp);
// blockIndex: block Index
// compoId: component ID whose BC we want to set.
// surface: which surface to set.
// type: boundary condition type.
// MacroVarsComp: which all macrovars are to be used in specifying the BC.
// valueMacroVarsComp: specified value for the boundary condition for the macro
// vars which are defined using MacroVarsComp.

void CalBlockCoordinates(const int blockIndex, Real* blockStartPos,
                         Real meshsize);
// blockIndex: Block Index
// blockStartPos: Starting position of the block in x, y, z directions to find
// coordinates. meshSize: Size of mesh (Assuming a constant in all 3 directions
// for stream collision scheme).

void KerSetCoordinates(const Real* coordX, const Real* coordY, const int* idx,
                       Real* coordinates);
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx, Real* coordinates);

void KerSetInitialMacroVarsHilemms(const Real* coordinates, const int* idx,
                                   Real* macroVars, Real* macroVarsInitVal,
                                   const int* componentId);
// Interface for Setting the initial values.
// coordinates: Coordinate array of nodes.
// idx: An array used by OPS which gives the index of the current grid point.

void AssignCoordinates(int blockIndex, Real* coordinates[SPACEDIM]);
// blockIndex: Block id.
// Coordinates: Array to store coordinates of various points inside the domain.
// This function was defined in setup_comput_domain and has been declared here
// as we are not using the preprocessor code separately.

void DefineIntialCond(int blockIndex, int componentId,
                      std::vector<Real> initialMacroValues);
// blockIndex: Block Id.
// componentId: Id of the component.
// initialMacroValues: Initial values of the macroscopic variables for a given
// component in a particluar block.

void SetupGeomPropAndNodeType(int blockIndex, BoundaryType* boundType);
void SetupDomainGeometryProperty(int blockIndex);
void SetupDomainNodeType(int blockIndex, VertexTypes* faceType);

// Functions to check for inclusion.
void DefineHaloNumber(int Halo_Number, int Halo_Depth, int Scheme_Halo_points,
                      int Num_Bound_Halo_Points);

// A wrapper Function which implements all the boundary conditions.
void ImplementBoundaryConditions();
#endif  // Hilemms_H