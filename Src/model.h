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

/** @brief Declare functions for constructing lattice models
 *  @author Jianping Meng
 **/
#ifndef MODEL_H
#define MODEL_H
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "type.h"

/**
 * @brief total number of discrete velocity/lattice
 */
extern int NUMXI;

/**
 * @brief LATTDIM the dimension of lattice model
 */
extern int LATTDIM;
/**
 * @brief CS a coefficient for the convience of writting unifed equilibrium
 *           function with different lattice
 * @details CS The equilibrium function can be written using \f$ \sqrt{R T_0}\f$
 *          as reference velocity where CS is the magnitude of a lattice
 *          velocity component in this reference system, see SHan, Yuan and Chen
 *          J. Fluid Mech. (2006)
 */
extern Real CS;
/**
 * @brief lattice set stored in a linear array
 * @details XI the storage looks like \f$ c^0_x \f$, \f$ c^0_y \f$, \f$ c^0_z \f$,
 * \f$ c^1_x \f$, \f$ c^1_y \f$,\f$ c^1_z \f$.... If there are multiple components
 * they are orgranised in the order when they are defined in the DefineComponent
 * call.
 */
extern Real* XI;
/**
 * @brief XIMAXVALUE maximum value of particle speed for calculating the CFL.
 */
extern Real XIMAXVALUE;
/*!
 * Quadrature weights
 */
extern Real* WEIGHTS;
/**
 * @brief NUMBCOMPONENTS Total number of components for multiple-component applications
 */
extern int NUMCOMPONENTS;
/**
 * @brief OPP The index of opposite directions of a discrete velocities
 * @details OPP It is only useful for regular cartesian mesh.
 * It may not be initialised for general structured mesh, as they needed to be
 * calculated according to the tangential line.
 */
extern int* OPP;

#include "model_host_device.h"

enum CollisionType {
    Collision_BGKIsothermal2nd = 0,
    Collision_BGKThermal4th = 1,
    Collision_BGKSWE4th = 2,
};

enum BodyForceType { BodyForce_1st = 1, BodyForce_None = 0 };

enum InitialType {Initial_BGKFeq2nd = 1};

struct MacroVariable {
    std::string name;
    int id{0};
    VariableTypes type;
};
struct Component {
    std::string name;
    int id{0};
    std::string latticeName;
    CollisionType collisionType;
    BodyForceType bodyForceType;
    InitialType initialType;
    std::map<VariableTypes, MacroVariable> macroVars;
    int index[2];
    int uId, vId;
#ifdef OPS_3D
    int wId;
#endif
    Real tauRef;
};

const std::map<int, Component>& g_Components();

inline const int ComponentNum() { return NUMCOMPONENTS; }
inline const int SizeF() { return NUMXI; }
inline const Real SoundSpeed() { return CS; }
inline const Real MaximumSpeed() { return XIMAXVALUE; }
/**
 * Free the pointer memory
 */
void DestroyModel();


void DefineComponents(const std::vector<std::string>& compoNames,
                      const std::vector<int>& compoId,
                      const std::vector<std::string>& lattNames,
                      const std::vector<Real> tauRef,
                      const SizeType timeStep = 0);


void DefineMacroVars(std::vector<VariableTypes> types,
                     std::vector<std::string> names,
                     std::vector<int> varId, std::vector<int> compoId,
                     const SizeType timeStep=0);

/*!
* Define collision terms for specified components
* Must be called after DefineComponets()
* HiLeMMS interface, , https://gitlab.com/jpmeng/hilemms
*/
void DefineCollision(std::vector<CollisionType> types,
                       std::vector<int> compoId);

void DefineBodyForce(std::vector<BodyForceType> types,
                     std::vector<SizeType> compoId);

void DefineInitialCondition(std::vector<InitialType> types,
                            std::vector<int> compoId);
#ifdef OPS_3D
void UpdateMacroVars3D();
void PreDefinedBodyForce3D();
void PreDefinedInitialCondition3D();
void PreDefinedCollision3D();
#endif
#ifdef OPS_2D
void UpdateMacroVars();
void PreDefinedBodyForce();
void PreDefinedInitialCondition();
void PreDefinedCollision();
#endif
#endif
