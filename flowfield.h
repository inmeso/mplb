// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
/*! @brief   Implementing functions related to the flow field
  * @author  Jianping Meng
  * @details Implementing functions related to create the flow
  * field (allocate memory), set up the geometry and the boundary
  * property, and deallocate the memory.
  */


#ifndef FLOWFIELD_H
#define FLOWFIELD_H
#include <cmath>
#include <algorithm>
#include <string>
#include "boundary.h"
#include "model.h"
#include "scheme.h"
#include "type.h"
/*!
 * This module is set for defining blocks and variables defined on a block
 * including distribution functions, macroscopic variables, node properties,
 * and relevant parameters.
 * The responsibilities including:
 * 1. Create all variables from files or annually written subroutines
 * 2. Initialise the required macroscopic variables and thereby the
 *    distribution functions.
 * 3. Provide some tools for accessing variables.
 */


extern std::string CASENAME;
/*!
 * SPACEDIM=2 for 2D 3 for three 3D
 */
extern int SPACEDIM;
/*!
 * Total number of blocks
 */
extern int g_BlockNum;
extern ops_block* g_Block;
/*!
 * The size of g_f in each node will be determined by the employed quadrature
 * and the model. For example, if we are simulating a two-phase flow, then the
 * size will be the product of NUMXI and NUMCOMPONENTS.
 */
extern ops_dat* g_f;
/*! might be changed to a local temporary variable
 * if we use some control routine in the main.cpp
 */
extern ops_dat* g_fStage;
extern ops_dat* g_feq;
/*!
 * Bodyforce, which is independent of the particle velocity
 */
extern ops_dat* g_Bodyforce;
/*!
 * g_MacroVars: for storing the macroscopic variables, to reduce
 * the complexity of calculating equilibrium, it will has a specific order
 */
extern ops_dat* g_MacroVars;
/*!
* Save the macroscopic variables at the previous step
* Typically used for steady flow.
*/
extern ops_dat* g_MacroVarsCopy;
/*!
 * g_dt: time step
 */
extern int g_HaloDepth;
extern Real g_dt;
/*!
 * KN: Knudsen number defined by the reference quantities
 * It must be a constant during the run time
 */
extern Real* KN;
/*!
 * Relaxation time
 * Depend on some macroscopic variables like rho,T
 */
extern ops_dat* g_Tau;

/*!
* the residual error for steady flows
* for each macroscopic variable, there are two values: the absolute
* and relative
* for each component of a vector, two values are allocated
*/
extern Real* g_ResidualError;
extern ops_reduction* g_ResidualErrorHandle;

// Boundary fitting mesh
// The following variables are introduced for
// implementing finite difference schemes
/*!
 * g_DiscreteConvectionTerm: for finite difference scheme
 */
extern ops_dat* g_DiscreteConvectionTerm;
/*! metrics structure
 * | xi_x  0 xi_y  1 |
 * | eta_x 2 eta_y 3 |
 *
 */
extern ops_dat* g_Metrics;

// Boundary fitting mesh

// Cutting cell
/*!
 * g_NodeType: boundary or fluid
 */
extern ops_dat* g_NodeType;
/*!
 * immersed solid? or the end point of the body.
 */
extern ops_dat* g_GeometryProperty;

/*!
 * Coordinate
 */
extern ops_dat* g_CoordinateXYZ;
// Cutting cell
/*
 * total number of halo series
 * Note: such as a boundary surface count as one
 */
extern int g_HaloNum;
extern ops_halo* g_Halos;
// Do we need to define many halo groups?
// note: we only need halo point for f
extern ops_halo_group g_HaloGroups;
extern int* g_BlockIterRngWhole;
extern int* g_BlockIterRngJmin;
extern int* g_BlockIterRngJmax;
extern int* g_BlockIterRngImin;
extern int* g_BlockIterRngImax;
extern int* g_BlockIterRngBulk;
extern int* g_BlockIterRngKmax;
extern int* g_BlockIterRngKmin;
/*!
 * The size of each block, i.e., each domain
 */
extern int* g_BlockSize;

/*!
 *Get the pointer pointing to the starting position of IterRng of this block
 *No NULL check for efficiency
 *Note: it looks that ops_par_loop call does not support const point.
 */
inline int* BlockIterRng(const int blockId, int* iterRng) {
    return &iterRng[blockId * 2 * SPACEDIM];
}
/*!
 * Return the starting position of memory in which we store the size of each
 * block
 */
inline const int* BlockSize(const int blockId) {
    return &g_BlockSize[blockId * SPACEDIM];
}

inline const int BlockNum() { return g_BlockNum; }
inline const int SpaceDim() { return SPACEDIM; }
inline const int HaloDepth() { return g_HaloDepth; }
/*!
 * Manually setup the flow field.
 */
void SetupFlowfield();
void SetupFlowfieldfromHdf5();
void DefineVariables();
void WriteFlowfieldToHdf5(const long timeStep);
void WriteDistributionsToHdf5(const long timeStep);
void WriteNodePropertyToHdf5(const long timeStep);
void DestroyFlowfield();
Real TotalMeshSize();
const int HaloPtNum();
#endif
