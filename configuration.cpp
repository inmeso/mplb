#include "configuration.h"
Configuration config;
using json = nlohmann::json;

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
    BoundaryType,
    {
        {BoundaryType_KineticDiffuseWall, "Boundary_KineticDiffuseWall"},
        {BoundaryType_ExtrapolPressure1ST, "Boundary_ExtrapolPressure1ST"},
        {BoundaryType_ExtrapolPressure2ND, "Boundary_ExtrapolPressure2ND"},
        {BoundaryType_Periodic, "Boundary_Periodic"},
        {BoundaryType_BounceBackWall, "Boundary_BounceBackWall"},
        {BoundaryType_FreeFlux, "Boundary_FreeFlux"},
        {BoundaryType_ZouHeVelocity, "Boundary_ZouHeVelocity"},
        {BoundaryType_EQMDiffuseRefl, "Boundary_EQMDiffuseREfl"},
    });

NLOHMANN_JSON_SERIALIZE_ENUM(
    EquilibriumType,
    {{Equilibrium_BGKIsothermal2nd, "Equilibrium_BGKIsothermal2nd"},
     {Equilibrium_BGKThermal4th, "Equilibrium_BGKThermal4th"},
     {Equilibrium_BGKSWE4th, "Equilibrium_BGKSWE4th"}});

NLOHMANN_JSON_SERIALIZE_ENUM(BodyForceType,
                             {{BodyForce_1st, "BodyForce_1st"},
                              {BodyForce_None, "BodyForce_None"}});

const Configuration* Config() { return &config; }

void ParseJson(json& jsonConfig) {
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

    if (jsonConfig["CompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the SpaceDim item into the configuration!\n");
        assert(jsonConfig["CompoIds"].is_null());
    } else {
        config.compoIds = jsonConfig["CompoIds"].get<std::vector<int>>();
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
        config.macroVarIds = jsonConfig["MacroVarIds"].get<std::vector<int>>();
    }

    if (jsonConfig["MacroCompoIds"].is_null()) {
        ops_printf(
            "Error! Please insert the MacroCompoIds item into the "
            "configuration!\n");
        assert(jsonConfig["MacroCompoIds"].is_null());
    } else {
        config.macroCompoIds =
            jsonConfig["MacroCompoIds"].get<std::vector<int>>();
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
        config.blockSize = jsonConfig["BlockSize"].get<std::vector<int>>();
    }

    if (jsonConfig["StartPos"].is_null()) {
        ops_printf(
            "Error! Please insert the StartPos item into the "
            "configuration!\n");
        assert(jsonConfig["StartPos"].is_null());
    } else {
        config.startPos = jsonConfig["StartPos"].get<std::vector<Real>>();
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
    int boundaryConditionNum{2 * config.spaceDim * config.blockNum};
    for (int bcIdx = 0; bcIdx < boundaryConditionNum;bcIdx++){
        


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
    json jsonConfig = json::parse(configString);
    ParseJson(jsonConfig);
}
