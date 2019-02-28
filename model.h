// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief Declare functions for constructing discrete velocity model
 *  @author Jianping Meng
 **/
#ifndef MODEL_H
#define MODEL_H
#include <cmath>
#include <string>
#include <vector>
#include "type.h"
/*!
 * In this module, most of variables will not change
 * when the code is running.
 */

/*!
 *NUMXI:total number of discrete velocity/lattice
 */
extern int NUMXI;

/*!
 *dimension of lattice model
 */
extern int LATTDIM;
/*!
 * Speed of Sound
 * Note: for multiple-lattice application, the lattice sound speed
 * must be same.
 */
extern Real CS;
/*!
 * The start index and end index of each component in the XI and WEIGHTS;
 */
extern int* COMPOINDEX;
/*!
 * XI: if using stream-collision scheme XI should integers actually\n
 * XI: if using finite difference scheme, it may be real and CS=1 accordingly\n
 */
extern Real* XI;
/*!
 * XIMAXVALUE: maximum value of particle speed, for calculating the CFL.
 */
extern Real XIMAXVALUE;
/*!
 * Quadrature weights
 */
extern Real* WEIGHTS;
/*!
 * Total number of components.
 * This is introduced for multiple-component fluid flow.
 */
extern int NUMCOMPONENTS;
/*!
 * OPP: the index of opposite directions of a discrete velocities\n
 * OPP: it is only useful for regular cartesian mesh\n
 * OPP: may not be initialised for complex mesh,as they needed to be
 * calculated\n
 * according to the tangential line\n
 */
extern int* OPP;
/*!
 * Number of macroscopic variables
 */
extern int NUMMACROVAR;
/*!
 * The type of macroscopic variables
 */
extern int* VARIABLETYPE;
/*!
 * Which component does the variable belong to
 */
extern int* VARIABLECOMPINDEX;
/*!
 * The start and end position of macroscopic variables of each
 * component
 */
extern int* VARIABLECOMPPOS;
/*!
 * Equilibrium function type
 */
extern int* EQUILIBRIUMTYPE;
/*!
 * Force function type
 */
extern int* FORCETYPE;
void SetLatticeName(const std::vector<std::string> latticeName);
const std::vector<std::string> LatticeName();
const std::vector<std::string> MacroVarName();
inline const int ComponentNum() { return NUMCOMPONENTS; }
inline const int MacroVarsNum() { return NUMMACROVAR; }
inline const int SizeF() { return NUMXI; }
inline const Real SoundSpeed() { return CS; }
inline const Real MaximumSpeed() { return XIMAXVALUE; }
/*!
 * In theory the size of Tau should be equal to NUMCOMPONENTS.
 */
inline const int SizeofTau() { return NUMCOMPONENTS; }
/*!
 * Specify how to calculate macroscopic variables from
 * distribution function
 */
void SetupMacroVars();
void DefineComponents(std::vector<std::string> compoNames,
                      std::vector<int> compoId,
                      std::vector<std::string> lattNames);

void DefineMacroVars(std::vector<VariableTypes> types,
                     std::vector<std::string> names, std::vector<int> varId,
                     std::vector<int> compoId);

void DefineEquilibrium(std::vector<EquilibriumType> types,
                       std::vector<int> compoId);

void DefineBodyForce(std::vector<BodyForceType> types,
                     std::vector<int> compoId);
/*
 * Define the local function for calculating the equilibrium
 * 2D BGK model
 */
Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real T = 1, const int polyOrder = 2);
/*
 * Define the local function for calculating the equilibrium
 * 3D BGK model
 */
Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real w = 0, const Real T = 1,
                const int polyOrder = 2);
/*
 * Define the local function for calculating the SWE equilibrium
 */
Real CalcSWEFeq(const int l, const Real h = 1, const Real u = 0,
                const Real v = 0, const int polyOrder = 2);

/*!
 * Free the pointer memory
 */
void DestroyModel();
// Kernel functions that will be called by ops_par_loop
/*!
 * Calculate the equilibrium function for normal fluids
 * Polynomial equilibrium function: upto the fourth order
 */
// Two-dimensional version
void KerCalcFeq(const int* nodeType, const Real* macroVars, Real* feq);
void KerCalcMacroVars(const int* nodeType, const Real* f, Real* macroVars);
void KerCalcBodyForce3D(const Real* time, const int* nodeType,
                        const Real* coordinates, const Real* macroVars,
                        Real* bodyForce);
/*!
 * @fn defining how to calculate the relaxation time
 * @param tauRef the reference relaxation time
 * @param macroVars the macroscopic variables
 * @param tau the calculated relaxation time
 */
void KerCalcTau(const int* nodeType, const Real* tauRef, const Real* macroVars,
                Real* tau);
// Three-dimensional version
// We have to create 2D and 3D version because of the difference
// of 2D and 3D OPS_ACC_MD2 macro
void KerCalcFeq3D(const int* nodeType, const Real* macroVars, Real* feq);
void KerCalcTau3D(const int* nodeType, const Real* tauRef, const Real* macroVars,
                  Real* tau);
void KerCalcMacroVars3D(const Real* dt, const int* nodeType, const Real* coordinates, const Real* f, Real* macroVars);
#endif
