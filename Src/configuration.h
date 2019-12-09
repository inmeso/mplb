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

/*! @brief  Head files for configuration
 * @author  Jianping Meng
 * @details Declaring date structure for holding configuration parameters.
 * The configuration currently follows the HiLeMMS definition.
 * Usage: use Config() to get the constant pointer of the data structure.
 */
#ifndef CONFIGURATION_H
#define CONFIGURATION_H
#include <fstream>
#include <iostream>
#include <streambuf>
#include <string>
#include <vector>
#include <map>
#include "json.hpp"
#include "type.h"

/** Structure for holding various input parameters
*/

struct BlockBoundaryConfig {
    SizeType blockIndex;
    SizeType componentID;
    std::vector<Real> givenVars;
    BoundarySurface boundarySurface;
    BoundaryType boundaryType;
    std::vector<VariableTypes> macroVarTypesatBoundary;
};

struct Configuration {
    std::string caseName;
    SizeType spaceDim{3};
    std::vector<std::string> compoNames;
    std::vector<SizeType> compoIds;
    std::vector<std::string> lattNames;
//The following is optional
//TODO to write Query() method very soon.
    std::vector<VariableTypes> macroVarTypes;
    std::vector<std::string> macroVarNames;
    std::vector<SizeType> macroVarIds;
    std::vector<SizeType> macroCompoIds;
    std::vector<CollisionType> CollisionTypes;
    std::vector<SizeType> CollisionCompoIds;
    std::vector<BodyForceType> bodyForceTypes;
    std::vector<SizeType> bodyForceCompoIds;
    std::vector<InitialType> initialTypes;
    std::vector<SizeType> initialConditionCompoId;
    SchemeType schemeType{Scheme_StreamCollision};
    SizeType blockNum{1};
    std::vector<SizeType> blockSize;
    std::vector<Real> startPos;
    Real meshSize;
    std::vector<Real> tauRef;
    bool transient{true};
    Real convergenceCriteria{-1};
    SizeType timeSteps{0};
    SizeType checkPeriod{1000};
    std::vector<BlockBoundaryConfig> blockBoundaryConditions;
};

/** Reading the parameters from a input file in the json format
 *  In the MPI mode, the whole input file will be broadcasted to all nodes by
 *  the root rank 0
 */
void ReadConfiguration(std::string& configFileName);

/** Get the pointer to the data structure holding various parameters
 */
const Configuration& Config();

#endif  // CONFIGURATION_H