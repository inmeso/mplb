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
 * @author  Jianping Meng
 * @details Define the functions for Json configuration input
 */
#include "configuration.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
Configuration config;
using json = nlohmann::json;
json jsonConfig;

// map Enum types to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM(VariableTypes,
                             {
                                 {Variable_Rho, "Variable_Rho"},
                                 {Variable_U, "Variable_U"},
                                 {Variable_V, "Variable_V"},
                                 {Variable_W, "Variable_W"},
                                 {Variable_T, "Variable_T"},
                                 {Variable_Qx, "Variable_Qx"},
                                 {Variable_Qy, "Variable_Qy"},
                                 {Variable_Qz, "Variable_Qz"},
                                 {Variable_U_Force, "Variable_U_Force"},
                                 {Variable_V_Force, "Variable_V_Force"},
                                 {Variable_W_Force, "Variable_W_Force"},
                             });

NLOHMANN_JSON_SERIALIZE_ENUM(
    VertexType, {
                    {VertexType::Inlet, "Inlet"},
                    {VertexType::OutLet, "OutLet"},
                    {VertexType::MDPeriodic, "MDPeriodic"},
                    {VertexType::FDPeriodic, "FDPeriodic"},
                    {VertexType::Symmetry, "Symmetry"},
                    {VertexType::Wall, "Wall"},
                    {VertexType::VirtualBoundary, "VirtualBoundary"},
                    {VertexType::ImmersedSolid, "ImmersedSolid"},
                    {VertexType::ImmersedBoundary, "ImmersedBoundary"},
                });

NLOHMANN_JSON_SERIALIZE_ENUM(
    BoundaryScheme,
    {
        {BoundaryScheme::KineticDiffuseWall, "KineticDiffuseWall"},
        {BoundaryScheme::ExtrapolPressure1ST, "ExtrapolPressure1ST"},
        {BoundaryScheme::ExtrapolPressure2ND, "ExtrapolPressure2ND"},
        {BoundaryScheme::MDPeriodic, "MDPeriodic"},
        {BoundaryScheme::FDPeriodic, "FDPeriodic"},
        {BoundaryScheme::BounceBack, "BounceBack"},
        {BoundaryScheme::FreeFlux, "FreeFlux"},
        {BoundaryScheme::ZouHeVelocity, "ZouHeVelocity"},
        {BoundaryScheme::EQMDiffuseRefl, "EQMDiffuseREfl"},
        {BoundaryScheme::None, "None"},
    });

NLOHMANN_JSON_SERIALIZE_ENUM(
    CollisionType, {{Collision_BGKIsothermal2nd, "Collision_BGKIsothermal2nd"},
                    {Collision_BGKThermal4th, "Collision_BGKThermal4th"},
                    {Collision_BGKSWE4th, "Collision_BGKSWE4th"}});

NLOHMANN_JSON_SERIALIZE_ENUM(BodyForceType,
                             {{BodyForce_1st, "BodyForce_1st"},
                              {BodyForce_None, "BodyForce_None"}});

NLOHMANN_JSON_SERIALIZE_ENUM(BoundarySurface,
                             {
                                 {BoundarySurface::Left, "Left"},
                                 {BoundarySurface::Right, "Right"},
                                 {BoundarySurface::Top, "Top"},
                                 {BoundarySurface::Bottom, "Bottom"},
#ifdef OPS_3D
                                 {BoundarySurface::Front, "Front"},
                                 {BoundarySurface::Back, "Back"},
                                 {BoundarySurface::LeftBack, "LeftBack"},
                                 {BoundarySurface::LeftFront, "LeftFront"},
                                 {BoundarySurface::RightBack, "RightBack"},
                                 {BoundarySurface::RightFront, "RightFront"},
                                 {BoundarySurface::TopBack, "TopBack"},
                                 {BoundarySurface::TopFront, "TopFront"},
                                 {BoundarySurface::BottomBack, "BottomBack"},
                                 {BoundarySurface::BottomFront, "BottomFront"},
#endif
                                 {BoundarySurface::LeftTop, "LeftTop"},
                                 {BoundarySurface::LeftBottom, "LeftBottom"},
                                 {BoundarySurface::RightTop, "RightTop"},
                                 {BoundarySurface::RightBottom, "RightBottom"},
                             });

NLOHMANN_JSON_SERIALIZE_ENUM(InitialType,
                             {
                                 {Initial_BGKFeq2nd, "Initial_BGKFeq2nd"},
                             });

NLOHMANN_JSON_SERIALIZE_ENUM(SchemeType,
                             {
                                 {Scheme_E1st2nd, "Scheme_E1st2nd"},
                                 {Scheme_StreamCollision,
                                  "Scheme_StreamCollision"},
                                 {Scheme_I1st2nd, " Scheme_I1st2nd"},
                             });

const Configuration& Config() { return config; }

const json& JsonConfig() { return jsonConfig; }

int GetBlockBoundaryConditionNum() {
    int num{0};
    std::string key{"BoundaryCondition" + std::to_string(num)};
    while (jsonConfig.contains(key) && (!jsonConfig[key].is_null())) {
        num++;
        key = "BoundaryCondition" + std::to_string(num);
    }
    return num;
}

void ParseJson() {
    Query(config.caseName, "CaseName");
    Query(config.spaceDim, "SpaceDim");
    Query(config.compoNames, "CompoNames");
    Query(config.lattNames, "LatticeName");
    Query(config.compoIds, "CompoIds");
    Query(config.tauRef, "TauRef");
    Query(config.macroVarNames, "MacroVarNames");
    Query(config.macroVarIds, "MacroVarIds");
    Query(config.macroCompoIds, "MacroCompoIds");
    Query(config.macroVarTypes, "MacroVarTypes");
    Query(config.CollisionTypes, "CollisionType");
    Query(config.CollisionCompoIds, "CollisionCompoIds");
    Query(config.initialTypes, "InitialType");
    Query(config.initialConditionCompoId, "InitialCompoIds");
    Query(config.bodyForceCompoIds, "BodyForceCompoId");
    Query(config.bodyForceTypes, "BodyForceType");
    Query(config.schemeType, "SchemeType");
    Query(config.blockIds, "BlockIds");
    Query(config.blockNames, "BlockNames");
    Query(config.blockSize, "BlockSize");

    Check(config.fromBlockIds, "FromBlockIds");
    Check(config.toBlockIds, "ToBlockIds");
    Check(config.fromBoundarySurface, "FromBoundarySurface");
    Check(config.toBoundarySurface, "ToBoundarySurface");
    Check(config.blockConnectionType, "BlockConnectionType");

    for (const auto id : config.blockIds) {
        std::vector<Real> pos;
        Query(pos, "StartPos", std::to_string(id));
        config.startPos.emplace(id, pos);
    }
    Query(config.checkPeriod, "CheckPeriod");
    Query(config.meshSize, "MeshSize");
    int boundaryConditionNum{GetBlockBoundaryConditionNum()};
    config.blockBoundaryConfig.resize(boundaryConditionNum);
    for (int bcIdx = 0; bcIdx < boundaryConditionNum; bcIdx++) {
        std::string bcName{"BoundaryCondition" + std::to_string(bcIdx)};
        if (jsonConfig[bcName].is_null()) {
            ops_printf(
                "Error! Please insert the %s item into the "
                "configuration!\n",
                bcName.c_str());
            assert(jsonConfig[bcName].is_null());
        } else {
            Query(config.blockBoundaryConfig[bcIdx].blockIndex, bcName,
                  "BlockIndex");
            Query(config.blockBoundaryConfig[bcIdx].componentID, bcName,
                  "ComponentId");
            Query(config.blockBoundaryConfig[bcIdx].boundarySurface, bcName,
                  "BoundarySurface");
            Query(config.blockBoundaryConfig[bcIdx].boundaryScheme, bcName,
                  "BoundaryScheme");
            Query(config.blockBoundaryConfig[bcIdx].givenVars, bcName,
                  "GivenVars");
            Query(config.blockBoundaryConfig[bcIdx].boundaryType, bcName,
                  "BoundaryType");
            Query(config.blockBoundaryConfig[bcIdx].macroVarTypesatBoundary,
                  bcName, "MacroVarTypesatBoundary");
        }
    }
    Query(config.currentTimeStep, "CurrentTimeStep");
    Query(config.transient, "Transient");

    if (config.transient) {
        Query(config.timeStepsToRun, "TimeStepsToRun");

    } else {
        Query(config.convergenceCriteria, "ConvergenceCriteria");
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

void GetConfigFileFromCmd(bool& findConfig, std::string& fileName,
                          const int argc, const char** argv) {
    findConfig = false;
    for (int i = 1; i < argc; i++) {
        const std::string arg{argv[i]};
        SizeType found{arg.find("Config=")};
        if (found == 0) {
            findConfig = true;
            fileName = arg.substr(found + 7);
            break;
        }
    }
}