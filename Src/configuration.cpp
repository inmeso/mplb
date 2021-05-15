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
    VertexType,
    {
        {VertexType::Inlet, "Inlet"},
        {VertexType::OutLet, "OutLet"},
        {VertexType::MDPeriodic, "MDPeriodic"},
        {VertexType::FDPeriodic, "FDPeriodic"},
        {VertexType::Symmetry, "Symmetry"},
        {VertexType::Wall, "Wall"},
        {VertexType::ImmersedSolid, "ImmersedSolid"},
        {VertexType::ImmersedBoundary, "ImmersedBoundary"}
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
    });

NLOHMANN_JSON_SERIALIZE_ENUM(
    CollisionType,
    {{Collision_BGKIsothermal2nd, "Collision_BGKIsothermal2nd"},
     {Collision_BGKThermal4th, "Collision_BGKThermal4th"},
     {Collision_BGKSWE4th, "Collision_BGKSWE4th"}});

NLOHMANN_JSON_SERIALIZE_ENUM(BodyForceType,
                             {{BodyForce_1st, "BodyForce_1st"},
                              {BodyForce_None, "BodyForce_None"}});

NLOHMANN_JSON_SERIALIZE_ENUM(BoundarySurface,
                             {{BoundarySurface::Left, "Left"},
                              {BoundarySurface::Right, "Right"},
                              {BoundarySurface::Top, "Top"},
                              {BoundarySurface::Bottom, "Bottom"}
#ifdef OPS_3D
                              ,
                              {BoundarySurface::Front, "Front"},
                              {BoundarySurface::Back, "Back"}
#endif
                             });

NLOHMANN_JSON_SERIALIZE_ENUM(InitialType,
                             {{Initial_BGKFeq2nd, "Initial_BGKFeq2nd"}});

NLOHMANN_JSON_SERIALIZE_ENUM(SchemeType, {{Scheme_E1st2nd, "Scheme_E1st2nd"},
                                          {Scheme_StreamCollision,
                                           "Scheme_StreamCollision"},
                                          {Scheme_I1st2nd, " Scheme_I1st2nd"}});

const Configuration& Config() { return config; }

const json& JsonConfig() { return jsonConfig; }

void ParseJson() {
    if (jsonConfig["CaseName"].is_null()) {
        ops_printf(
            "Error! Please insert the CaseName item into the configuration!\n");
        assert(jsonConfig["CaseName"].is_null());
    } else {
        config.caseName = jsonConfig["CaseName"];
    }

    if (jsonConfig["SpaceDim"].is_null()) {
        ops_printf(
            "Error! Please insert the SpaceDim item into the configuration!\n");
        assert(jsonConfig["SpaceDim"].is_null());
    } else {
        config.spaceDim = jsonConfig["SpaceDim"];
    }

    if (jsonConfig["CompoNames"].is_null()) {
        ops_printf(
            "Error! Please insert the SpaceDim item into the configuration!\n");
        assert(jsonConfig["CompoNames"].is_null());
    } else {
        config.compoNames =
            jsonConfig["CompoNames"].get<std::vector<std::string>>();
    }

    if (jsonConfig["LatticeName"].is_null()) {
        ops_printf(
            "Error! Please insert the LatticeName item into the "
            "configuration!\n");
        assert(jsonConfig["LatticeName"].is_null());
    } else {
        config.lattNames =
            jsonConfig["LatticeName"].get<std::vector<std::string>>();
    }

    if (jsonConfig["CompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the SpaceDim item into the configuration!\n");
        assert(jsonConfig["CompoIds"].is_null());
    } else {
        config.compoIds = jsonConfig["CompoIds"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["MacroVarNames"].is_null()) {
        ops_printf(
            "Error! Please insert the MacroVarNames item into the "
            "configuration!\n");
        assert(jsonConfig["MacroVarNames"].is_null());
    } else {
        config.macroVarNames =
            jsonConfig["MacroVarNames"].get<std::vector<std::string>>();
    }

    if (jsonConfig["MacroVarIds"].is_null()) {
        ops_printf(
            "Error! Please insert the MacroVarIds item into the "
            "configuration!\n");
        assert(jsonConfig["MacroVarIds"].is_null());
    } else {
        config.macroVarIds = jsonConfig["MacroVarIds"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["MacroCompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the MacroCompoIds item into the "
            "configuration!\n");
        assert(jsonConfig["MacroCompoIds"].is_null());
    } else {
        config.macroCompoIds =
            jsonConfig["MacroCompoIds"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["MacroVarTypes"].is_null()) {
        ops_printf(
            "Error! Please insert the MacroVarTypes item into the "
            "configuration!\n");
        assert(jsonConfig["MacroVarTypes"].is_null());
    } else {
        config.macroVarTypes =
            jsonConfig["MacroVarTypes"].get<std::vector<VariableTypes>>();
    }

    if (jsonConfig["CollisionType"].is_null()) {
        ops_printf(
            "Error! Please insert the CollisionType item into the "
            "configuration!\n");
        assert(jsonConfig["CollisionType"].is_null());
    } else {
        config.CollisionTypes =
            jsonConfig["CollisionType"].get<std::vector<CollisionType>>();
    }

    if (jsonConfig["CollisionCompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the EquilibriumCompoIds item into the "
            "configuration!\n");
        assert(jsonConfig["CollisionCompoIds"].is_null());
    } else {
        config.CollisionCompoIds =
            jsonConfig["CollisionCompoIds"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["InitialType"].is_null()) {
        ops_printf(
            "Error! Please insert the InitialType item into the "
            "configuration!\n");
        assert(jsonConfig["InitialType"].is_null());
    } else {
        config.initialTypes =
            jsonConfig["InitialType"].get<std::vector<InitialType>>();
    }

    if (jsonConfig["InitialCompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the InitialCompoIds item into the "
            "configuration!\n");
        assert(jsonConfig["InitialCompoIds"].is_null());
    } else {
        config.initialConditionCompoId =
            jsonConfig["InitialCompoIds"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["BodyForceType"].is_null()) {
        ops_printf(
            "Error! Please insert the BodyForceType item into the "
            "configuration!\n");
        assert(jsonConfig["BodyForceType"].is_null());
    } else {
        config.bodyForceTypes =
            jsonConfig["BodyForceType"].get<std::vector<BodyForceType>>();
    }

    if (jsonConfig["BodyForceCompoId"].is_null()) {
        ops_printf(
            "Error! Please insert the BodyForceCompoId item into the "
            "configuration!\n");
        assert(jsonConfig["BodyForceCompoId"].is_null());
    } else {
        config.bodyForceCompoIds =
            jsonConfig["BodyForceCompoId"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["SchemeType"].is_null()) {
        ops_printf(
            "Error! Please insert the SchemeType item into the "
            "configuration!\n");
        assert(jsonConfig["SchemeType"].is_null());
    } else {
        config.schemeType = jsonConfig["SchemeType"];
    }

    SizeType boundaryConditionNum{2 * config.spaceDim * config.blockNum};
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
            if (jsonConfig[bcName]["BlockIndex"].is_null()) {
                ops_printf("Error! Please insert the Block item into %s\n",
                           bcName.c_str());
                assert(jsonConfig[bcName]["BlockIndex"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].blockIndex =
                    jsonConfig[bcName]["BlockIndex"];
            }

            if (jsonConfig[bcName]["ComponentId"].is_null()) {
                ops_printf(
                    "Error! Please insert the ComponentId item into %s\n",
                    bcName.c_str());
                assert(jsonConfig[bcName]["ComponentId"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].componentID =
                    jsonConfig[bcName]["ComponentId"];
            }

            if (jsonConfig[bcName]["BoundarySurface"].is_null()) {
                ops_printf(
                    "Error! Please insert the ComponentId item into %s\n",
                    bcName.c_str());
                assert(jsonConfig[bcName]["BoundarySurface"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].boundarySurface =
                    jsonConfig[bcName]["BoundarySurface"];
            }

            if (jsonConfig[bcName]["BoundaryScheme"].is_null()) {
                ops_printf(
                    "Error! Please insert the BoundaryScheme item into %s\n",
                    bcName.c_str());
                assert(jsonConfig[bcName]["BoundaryScheme"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].boundaryScheme =
                    jsonConfig[bcName]["BoundaryScheme"].get<BoundaryScheme>();
            }

            if (jsonConfig[bcName]["BoundaryType"].is_null()) {
                ops_printf(
                    "Error! Please insert the BoundaryType item into %s\n",
                    bcName.c_str());
                assert(jsonConfig[bcName]["BoundaryType"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].boundaryType =
                    jsonConfig[bcName]["BoundaryType"].get<VertexType>();
            }

            if (jsonConfig[bcName]["GivenVars"].is_null()) {
                ops_printf("Error! Please insert the GivenVars item into %s\n",
                           bcName.c_str());
                assert(jsonConfig[bcName]["GivenVars"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].givenVars =
                    jsonConfig[bcName]["GivenVars"].get<std::vector<Real>>();
            }
            if (jsonConfig[bcName]["MacroVarTypesatBoundary"].is_null()) {
                ops_printf(
                    "Error! Please insert the MacroVarTypesatBoundary item "
                    "into %s\n",
                    bcName.c_str());
                assert(jsonConfig[bcName]["MacroVarTypesatBoundary"].is_null());
            } else {
                config.blockBoundaryConfig[bcIdx].macroVarTypesatBoundary =
                    jsonConfig[bcName]["MacroVarTypesatBoundary"]
                        .get<std::vector<VariableTypes>>();
            }
        }
    }

    if (jsonConfig["BlockNum"].is_null()) {
        ops_printf(
            "Error! Please insert the BlockNum item into the "
            "configuration!\n");
        assert(jsonConfig["BlockNum"].is_null());
    } else {
        config.blockNum = jsonConfig["BlockNum"];
    }

    if (jsonConfig["BlockSize"].is_null()) {
        ops_printf(
            "Error! Please insert the BlockSize item into the "
            "configuration!\n");
        assert(jsonConfig["BlockSize"].is_null());
    } else {
        config.blockSize = jsonConfig["BlockSize"].get<std::vector<SizeType>>();
    }

    if (jsonConfig["TauRef"].is_null()) {
        ops_printf(
            "Error! Please insert the TauRef item into the "
            "configuration!\n");
        assert(jsonConfig["TauRef"].is_null());
    } else {
        config.tauRef = jsonConfig["TauRef"].get<std::vector<Real>>();
    }

    if (jsonConfig["StartPos"].is_null()) {
        ops_printf(
            "Error! Please insert the StartPos item into the "
            "configuration!\n");
        assert(jsonConfig["StartPos"].is_null());
    } else {
        config.startPos = jsonConfig["StartPos"].get<std::vector<Real>>();
    }

    if (jsonConfig["CheckPeriod"].is_null()) {
        ops_printf(
            "Error! Please insert the CheckPeriod item into the "
            "configuration!\n");
        assert(jsonConfig["CheckPeriod"].is_null());
    } else {
        config.checkPeriod = jsonConfig["CheckPeriod"];
    }

    if (jsonConfig["MeshSize"].is_null()) {
        ops_printf(
            "Error! Please insert the MeshSize item into the "
            "configuration!\n");
        assert(jsonConfig["MeshSize"].is_null());
    } else {
        config.meshSize = jsonConfig["MeshSize"];
    }

    if (jsonConfig["Transient"].is_null()) {
        ops_printf(
            "Error! Please insert the Transient item into the "
            "configuration!\n");
        assert(jsonConfig["Transient"].is_null());
    } else {
        config.transient = jsonConfig["Transient"];
    }
    if (config.transient) {
        if (jsonConfig["TimeSteps"].is_null()) {
            ops_printf(
                "Error! Please insert the TimeSteps item into the "
                "configuration!\n");
            assert(jsonConfig["TimeSteps"].is_null());
        } else {
            config.timeSteps = jsonConfig["TimeSteps"];
        }

    } else {
        if (jsonConfig["ConvergenceCriteria"].is_null()) {
            ops_printf(
                "Error! Please insert the ConvergenceCriteria item into the "
                "configuration!\n");
            assert(jsonConfig["ConvergenceCriteria"].is_null());
        } else {
            config.convergenceCriteria = jsonConfig["ConvergenceCriteria"];
        }
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