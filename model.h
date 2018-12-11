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
 * In this module, most of variables will not change the value when the code
 * is running.
 */
extern std::string MODELNAME;
/*!
 *NUMXI:total number of discrete velocity/lattice
 */
extern int NUMXI;
/*!
 *Order of feq
 */
extern int FEQORDER;
/*!
 *dimension of lattice model
 */
extern int LATTDIM;
/*!
 * Speed of Sound
 */
extern Real CS;

/*!
 * Is the case a thermal problem?
 * 0 isothermal
 * 1 thermal
 */
extern int THERMALPROBLEM;
/*!
 * The start and end  of each component in the XI and WEIGHTS;
 *
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
/*!
 * Total number of components.
 * This is introduced for multiple-component fluid flow.
 */
extern int NUMCOMPONENTS;
extern Real* WEIGHTS;
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

/*
 * The type of macroscopic variables
 */

extern int* VARIABLETYPE;
/*
 * Which component does the variable belong to
 */
extern int* VARIABLECOMPINDEX;

/*!
 *Name of all macroscopic variables
 */
extern std::vector<std::string> MACROVARNAME;
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
/*!
 * Setup the D2Q9 model
 */
void SetupD2Q9Latt();
void SetupD2Q16Latt();
/*!
 * Setup a general Gauss-Hermite lattice model
 */
void SetupGenGsHermLatt(const int quadratureOrder);
/*!
 * Setup a lattice from an input file
 */
void SetupLattFromFile(const std::string filename);
/*!
 * Calculate the equilibrium function
 */
void KerCutCellCalcFeqIso(const int* nodeType, const Real* macroVars,
                          Real* feq);
void KerCutCellCalcPolyFeq(const int* polyOrder, const int* nodeType,
                           const Real* macroVars, Real* feq);
void KerCalcMacroVars(const int* nodeType, const Real* f, Real* macroVars);
void KerCalcBodyForce(const int* nodeType, const Real* f, const Real* macroVars,
                      Real* bodyForce);
/*!

*/
void KerCutCellCalcSWEFeq(const int* nodeType, const Real* macroVars,
                          Real* feq);
void KerCutCellCalcPolySWEFeq(const int* polyOrder, const int* nodeType,
                              const Real* macroVars, Real* feq);

void KerCalcSWETau(const int* nodeType, const Real* kn, const Real* macroVars,
                   Real* tau);

/*
 * Define the local function for calculating the equilibrium function based on
 * the BGK model
 */
Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real T = 1, const int polyOrder = 2);

Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real w = 0, const Real T = 1,
                const int polyOrder = 2);
/*
 * Define the local function for calculating the SWE equilibrium function based
 * on the BGK model
 */
Real CalcSWEFeq(const int l, const Real h = 1, const Real u = 0,
                const Real v = 0, const int polyOrder = 2);
/*!
 * @fn defining how to calculate the relaxation time
 * @param kn the Knudsen number
 * @param macroVars the macroscopic variables
 * @param tau the calculated relaxation time
 */
void KerCalcTau(const int* nodeType, const Real* kn, const Real* macroVars,
                Real* tau);

void SetupModel();
/*!
 * Free the pointer memory other than lattice XI
 */
void DestroyModel();
/*!
 * Free the pointer memory of XI related
 */
void DestroyLatt();
void KerCutCellCalcPolyFeq3D(const int* polyOrder, const int* nodeType,
                             const Real* macroVars, Real* feq);
void KerCalcTau3D(const int* nodeType, const Real* kn, const Real* macroVars,
                  Real* tau);
void KerCalcMacroVars3D(const int* nodeType, const Real* f, Real* macroVars);
#endif
