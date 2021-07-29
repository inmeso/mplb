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

/*! @brief  Function definition for configuration
 * @author  Sina Haeri (s.haeri@ed.ac.uk)
 * @details Define the functions for Json configuration input for Mie 
 * calculations.
 */
#include "MieConfig.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif

MieConfig config;
using json = nlohmann::json;
json jsonConfig;

namespace std {
    
    template< class T > void to_json(json &j, const std::complex< T > &p) {
        j = json{ {"real", p.real()}, {"imag", p.imag()} };
    }
    
    template< class T > void from_json(const json &j, std::complex< T > &p) {
        p.real(j.at("real"));
        p.imag(j.at("imag"));
    }
}

const MieConfig& Config() { return config; }

const json& JsonConfig() { return jsonConfig; }

void ParseJson() {
    if (jsonConfig["CaseName"].is_null()) {
        ops_printf(
            "Error! Please insert the CaseName item into the configuration!\n");
        assert(jsonConfig["CaseName"].is_null());
    } else {
        config.caseName = jsonConfig["CaseName"];
    }

    if (jsonConfig["spaceDim"].is_null()) {
        ops_printf(
            "Error! Please insert the spaceDim item into the configuration!\n");
        assert(jsonConfig["spaceDim"].is_null());
    } else {
        config.spaceDim = jsonConfig["spaceDim"];
    }

    if (jsonConfig["blockSize"].is_null()) {
        ops_printf(
            "Error! Please insert the blockSize item into the "
            "configuration!\n");
        assert(jsonConfig["blockSize"].is_null());
    } else {
        config.blockSize =
            jsonConfig["blockSize"].get<std::vector<int>>();
    }

    if (jsonConfig["blockExtent"].is_null()) {
        ops_printf(
            "Error! Please insert the blockExtent item into the "
            "configuration!\n");
        assert(jsonConfig["blockExtent"].is_null());
    } else {
        config.blockExtent =
            jsonConfig["blockExtent"].get<std::vector<Real>>();
    }

    if (jsonConfig["partRadius"].is_null()) {
        ops_printf(
            "Error! Please insert the partRadius item into the configuration!\n");
        assert(jsonConfig["partRadius"].is_null());
    } else {
        config.partRadius =
            jsonConfig["partRadius"].get<Real>();
    }

    if (jsonConfig["partPermeability"].is_null()) {
        ops_printf(
            "Error! Please insert the partPermeability (particle Permeability) item into the "
            "configuration!\n");
        assert(jsonConfig["partPermeability"].is_null());
    } else {
        config.partPermeability =
            jsonConfig["partPermeability"].get<Real>();
    }

    if (jsonConfig["envPermeability"].is_null()) {
        ops_printf(
            "Error! Please insert the envPermeability item into the configuration!\n");
        assert(jsonConfig["envPermeability"].is_null());
    } else {
        config.envPermeability = jsonConfig["envPermeability"].get<Real>();
    }

    if (jsonConfig["partRefractiveIndex"].is_null()) {
        ops_printf(
            "Error! Please insert the partRefractiveIndex item into the "
            "configuration!\n");
        assert(jsonConfig["partRefractiveIndex"].is_null());
    } else {
        config.partRefractiveIndex =
            jsonConfig["partRefractiveIndex"].get<std::complex<Real>>();
    }

    if (jsonConfig["envRefractiveIndex"].is_null()) {
        ops_printf(
            "Error! Please insert the envRefractiveIndex item into the "
            "configuration!\n");
        assert(jsonConfig["envRefractiveIndex"].is_null());
    } else {
        config.envRefractiveIndex = jsonConfig["envRefractiveIndex"].get<std::complex<Real>>();
    }

    if (jsonConfig["vaccumWaveLength"].is_null()) {
        ops_printf(
            "Error! Please insert the vaccumWaveLength item into the "
            "configuration!\n");
        assert(jsonConfig["vaccumWaveLength"].is_null());
    } else {
        config.vaccumWaveLength =
            jsonConfig["vaccumWaveLength"].get<Real>();
    }
}

void ReadConfiguration(std::string& configFileName) {
    std::string configString;
#ifdef OPS_MPI
    long configFileSize{0};
    if (ops_my_global_rank == MPI_ROOT) {
        std::ifstream configFile(configFileName);
        if (!configFile.is_open()) {
            ops_printf("Error! Cannot open the configuration file %s\n",
                       configFileName.c_str());
            assert(configFile.is_open());
        }
        std::string tmpStr((std::istreambuf_iterator<char>(configFile)),
                           std::istreambuf_iterator<char>());
        configString = tmpStr;
        configFileSize = configString.length();
    }
    MPI_Bcast(&configFileSize, 1, MPI_INT, 0, OPS_MPI_GLOBAL);
    char* strBuf = new char[configFileSize + 1];
    if (ops_my_global_rank == 0) {
        configString.copy(strBuf, configFileSize + 1);
        strBuf[configFileSize] = '\0';
    }
    MPI_Bcast(strBuf, configFileSize + 1, MPI_CHAR, 0, OPS_MPI_GLOBAL);
    if (ops_my_global_rank != MPI_ROOT) {
        configString = strBuf;
    }

    FreeArrayMemory(strBuf);
#else
    std::ifstream configFile(configFileName);
    if (!configFile.is_open()) {
        ops_printf("Error! Cannot open the configuration file %s\n",
                   configFileName.c_str());
        assert(configFile.is_open());
    }
    std::string tmpStr((std::istreambuf_iterator<char>(configFile)),
                       std::istreambuf_iterator<char>());
    configString = tmpStr;
#endif  // OPS_MPI
    jsonConfig = json::parse(configString);
    ParseJson();
}