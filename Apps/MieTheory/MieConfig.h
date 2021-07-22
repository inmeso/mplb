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
 * @author  Sina Haeri (s.haeri@ed.ac.uk)
 * @details Declaring date structure for holding configuration parameters
 * for Mie calculations
 */
#ifndef MIECONFIG_H
#define MIECONFIG_H
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

/** Structure for holding various input parameters
*/
// MAKE SURE YOU BUILD THE CODE IN DOUBLE PRECISION (REAL 
// should at least be a double) 
struct MieConfig {
    std::string caseName;
    int spaceDim{3};
    std::vector<int> blockSize;
    Real partRadius;
    Real partPermeability;
    Real envPermeability;
    std::complex<Real> partRefractiveIndex;
    std::complex<Real> envRefractiveIndex;
    Real vaccumWaveLength;
};

/** Reading the parameters from a input file in the json format
 *  In the MPI mode, the whole input file will be broadcasted to all nodes by
 *  the root rank 0
 */
void ReadConfiguration(std::string& configFileName);

/** Get the pointer to the data structure holding various parameters
 */
const MieConfig& Config();
const nlohmann::json& JsonConfig();

template <typename T>
void Query(T& value, std::string key) {
    const nlohmann::json& jsonConfig{JsonConfig()};
    if (jsonConfig[key].is_null()) {
        ops_printf("Error! Please insert the %s item into the configuration!\n",
                   key.c_str());
        assert(jsonConfig[key].is_null());
    };
    value = jsonConfig[key].get<T>();
}

#endif  // MIECONFIG_H