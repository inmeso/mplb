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

/*! @brief Declare functions for constructing discrete velocity model
 *  @author Jianping Meng
 **/
#ifndef MODEL_H
#define MODEL_H
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include "type.h"

/*!
 * Most of variables in this module will not change when the code is running.
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
 * Note: for multiple-component applications, the sound speed must be same.
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
 * Total number of components for multiple-component applications
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
 * Force function type
 */
extern int* FORCETYPE;
// Convenient functions
void SetLatticeName(const std::vector<std::string>& latticeName);
const std::vector<std::string> LatticeName();
const std::vector<std::string> MacroVarName();
/*!
 * Get collision type
 */
const std::list<std::pair<SizeType,CollisionType >> & CollisionTerms();

inline const int ComponentNum() { return NUMCOMPONENTS; }
inline const int MacroVarsNum() { return NUMMACROVAR; }
inline const int SizeF() { return NUMXI; }
inline const Real SoundSpeed() { return CS; }
inline const Real MaximumSpeed() { return XIMAXVALUE; }
/*!
 * Free the pointer memory
 */
void DestroyModel();
/*!
 * In theory the size of Tau should be equal to NUMCOMPONENTS.
 */
inline const int SizeofTau() { return NUMCOMPONENTS; }
// HiLeMMS interface, https://gitlab.com/jpmeng/hilemms
void DefineComponents(std::vector<std::string> compoNames,
                      std::vector<int> compoId,
                      std::vector<std::string> lattNames);

void DefineMacroVars(std::vector<VariableTypes> types,
                     std::vector<std::string> names, std::vector<int> varId,
                     std::vector<int> compoId);

/*!
* Define collision terms for specified components
* Must be called after DefineComponets()
* HiLeMMS interface, , https://gitlab.com/jpmeng/hilemms
*/
void DefineCollision(std::vector<CollisionType> types,
                       std::vector<SizeType> compoId);

void DefineBodyForce(std::vector<BodyForceType> types,
                     std::vector<int> compoId);
/*
 * Local function for calculating the equilibrium
 * 2D BGK model including up to fourth order terms
 * Please refer to Shan, Yuan and Chen JFM 2006(550):413-441
 */
Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real T = 1, const int polyOrder = 2);
/*
 * Local function for calculating the equilibrium
 * 3D BGK model including up to fourth order terms
 * Please refer to Shan, Yuan and Chen JFM 2006(550):413-441
 */
Real CalcBGKFeq(const int l, const Real rho = 1, const Real u = 0,
                const Real v = 0, const Real w = 0, const Real T = 1,
                const int polyOrder = 2);
/*
 * Local function for calculating the SWE equilibrium
 * Including up to fourth order terms
 * Please refer to Meng, Gu Emerson, Peng and Zhang, IJMPC 2018(29):1850080
 */
Real CalcSWEFeq(const int l, const Real h = 1, const Real u = 0,
                const Real v = 0, const int polyOrder = 2);

// Kernel functions that will be called by ops_par_loop
/*!
 * Calculate the equilibrium function for normal fluids
 * Polynomial equilibrium function: upto the fourth order
 */
// Two-dimensional version
void KerCalcFeq(const int* nodeType, const Real* macroVars, Real* feq);
void KerCalcMacroVars(const int* nodeType, const Real* f, Real* macroVars);

// Three-dimensional version
// We have to create 2D and 3D version because of the difference
// of 2D and 3D OPS_ACC_MD2 macro
void KerCalcBodyForce3D(const Real* time, const int* nodeType,
                        const Real* coordinates, const Real* macroVars,
                        Real* bodyForce);
void KerCalcBodyForce(const Real* time, const int* nodeType,
                      const Real* coordinates, const Real* macroVars,
                      Real* bodyForce);
void KerCalcFeq3D(const int* nodeType, const Real* macroVars, Real* feq);

void KerCalcMacroVars3D(const Real* dt, const int* nodeType, const Real* coordinates, const Real* f, Real* macroVars);

/**
 * @brief Implement the BGK isothermal collision model
 *
 * @param fStage the temporary variable storing distribution after collision
 * @param f distribution
 * @param macroVars macroscopic variables
 * @param nodeType
 * @param tauRef relaxation time
 * @param dt time step
 * @param componentId the component to be working on
 * Assumptions:
 * 1. fStage before collision is body force distribution
 * 2. The layout of macroVars is "rho, u, v, w"
 */
void KerCollideBGKIsothermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                  const ACC<Real>& macroVars, const ACC<int>& nodeType,
                  const Real* tauRef, const Real* dt, const int* componentId);

/**
 * @brief Implement the BGK thermal collision model
 *
 * @param fStage the temporary variable storing distribution after collision
 * @param f distribution
 * @param macroVars macroscopic variables
 * @param nodeType
 * @param tauRef relaxation time
 * @param dt time step
 * @param componentId componentId the component to be working on
 * 1. fStage before collision is body force distribution
 * 2. The layout of macroVars is "rho, u, v, w, T"
*/
void KerCollideBGKThermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                  const ACC<Real>& macroVars, const ACC<int>& nodeType,
                  const Real* tauRef, const Real* dt, const int* componentId)

#endif
