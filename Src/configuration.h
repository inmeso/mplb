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
 * @details Declaring data structure for holding configuration parameters.
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
#include <type_traits>
#include "json.hpp"
#include "type.h"
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "model.h"
#include "model_host_device.h"
#include "flowfield_host_device.h"
#include "boundary.h"

/**
 * Structure for holding various input parameters.
 */

struct Configuration {
    std::string caseName;
    SizeType spaceDim{3};
    std::vector<std::string> compoNames;
    std::vector<int> compoIds;
    std::vector<std::string> lattNames;
    std::vector<VariableTypes> macroVarTypes;
    std::vector<std::string> macroVarNames;
    std::vector<int> macroVarIds;
    std::vector<int> macroCompoIds;
    std::vector<CollisionType> CollisionTypes;
    std::vector<int> CollisionCompoIds;
    std::vector<BodyForceType> bodyForceTypes;
    std::vector<SizeType> bodyForceCompoIds;
    std::vector<InitialType> initialTypes;
    std::vector<int> initialConditionCompoId;
    SchemeType schemeType{Scheme_StreamCollision};
    std::vector<std::string> blockNames;
    std::vector<int> blockIds;
    std::vector<int> blockSize;
    std::vector<int> fromBlockIds;
    std::vector<int> toBlockIds;
    std::vector<BoundarySurface> fromBoundarySurface;
    std::vector<BoundarySurface> toBoundarySurface;
    std::vector<VertexType> blockConnectionType;
    std::map<int, std::vector<Real>> startPos;
    Real meshSize;
    std::vector<Real> tauRef;
    bool transient{true};
    Real convergenceCriteria{-1};
    SizeType timeStepsToRun{0};
    SizeType currentTimeStep{0};
    SizeType checkPeriod{1000};
    std::vector<BlockBoundary> blockBoundaryConfig;
};
/**
 * @brief Reading the parameters from a input file in the json format
 *
 * @details In the MPI mode, the whole input file will be broadcasted to all
 * nodes by the root rank 0
 *
 * @param configFileName configuration file name
 */
void ReadConfiguration(std::string& configFileName);

/** Get the pointer to the data structure holding various parameters
 */
const Configuration& Config();

const nlohmann::json& JsonConfig();

template <typename T>
void Query(T& value, const std::string& key) {
    const nlohmann::json& jsonConfig{JsonConfig()};
    if (jsonConfig.contains(key)) {
        if (jsonConfig[key].is_null()) {
            ops_printf(
                "Error! Please insert the %s item into the configuration!\n",
                key.c_str());
            assert(jsonConfig[key].is_null());
        };
        value = jsonConfig[key].get<T>();
    } else {
        ops_printf("Error! Please supply %s in the configuration!\n",
                   key.c_str());
        assert(jsonConfig.contains(key));
    }
}
/**
 * @brief Get non-mandatory options from JSON configuration file
 *
 * @tparam T Option type
 * @param value Option value
 * @param key Option name
 */
template <typename T>
void Check(T& value, const std::string& key) {
    const nlohmann::json& jsonConfig{JsonConfig()};
    bool haveKey{true};
    if (!jsonConfig.contains(key)) {
        ops_printf("Warning! %s is not defined in the configuration!\n",
                   key.c_str());
        haveKey = false;
    } else {
        if (jsonConfig[key].is_null()) {
            ops_printf("Warning! %s is empty in the configuration!\n",
                       key.c_str());
            haveKey = false;
        };
    }
    if (haveKey) {
        value = jsonConfig[key].get<T>();
    };
}

template <typename T>
void Query(T& value, const std::string& key0, const std::string& key1) {
    const nlohmann::json& jsonConfig{JsonConfig()};
    if (!jsonConfig.contains(key0)) {
        ops_printf("Error! Please supply %s in the configuration!\n",
                   key0.c_str());
        assert(jsonConfig.contains(key0));
    }
    if (!jsonConfig[key0].contains(key1)) {
        ops_printf("Error! Please supply %s : %s in the configuration!\n",
                   key0.c_str(), key1.c_str());
        assert(jsonConfig[key0].contains(key1));
    }
    if (jsonConfig[key0][key1].is_null()) {
        ops_printf(
            "Error! Please insert the %s->%s item into the configuration!\n",
            key0.c_str(), key1.c_str());
        assert(jsonConfig[key0][key1].is_null());
    };
    value = jsonConfig[key0][key1].get<T>();
}

/**
 * @brief Get the Config File From Cmd object.
 * @param findConfig if a configuration file is specified.
 * @param fileName the configuration file name.
 * @param argc the command line argument number.
 * @param argv the command line arguments.
 */
void GetConfigFileFromCmd(bool& findConfig, std::string& fileName,
                          const int argc, const char** argv);

#endif  // CONFIGURATION_H