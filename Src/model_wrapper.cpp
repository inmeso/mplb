#include <vector>
#include <map>
#include "flowfield.h"
#include "flowfield_host_device.h"
#include "model.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "model_kernel.inc"
#ifdef OPS_3D
void PreDefinedCollision3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKIsothermal3D, "KerCollideBGKIsothermal3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKThermal4th:
                    ops_par_loop(
                        KerCollideBGKThermal3D, "KerCollideBGKThermal3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_T).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified collision type is not implemented!\n");
                    break;
            }
        }
    }
#endif // OPS_3D
}

void UpdateMacroVars3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            for (auto& macroVar : compo.macroVars) {
                const int varId{macroVar.second.id};
                const VariableTypes varType{macroVar.first};
                switch (varType) {
                    case Variable_Rho:
                        ops_par_loop(
                            KerCalcDensity3D, "KerCalcDensity3D", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_U:
                        ops_par_loop(
                            KerCalcU3D, "KerCalcU3D", block.Get(), SpaceDim(),
                            iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_V:
                        ops_par_loop(
                            KerCalcV3D, "KerCalcV3D", block.Get(), SpaceDim(),
                            iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_W:
                        ops_par_loop(
                            KerCalcW3D, "KerCalcW3D", block.Get(), SpaceDim(),
                            iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_U_Force:
                        ops_par_loop(
                            KerCalcUForce3D, "KerCalcUForce3D", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_V_Force:
                        ops_par_loop(
                            KerCalcVForce3D, "KerCalcVForce3D", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_W_Force:
                        ops_par_loop(
                            KerCalcWForce3D, "KerCalcWForce3D", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    default:
                        break;
                }
            }
        }
    }
#endif // OPS_3D
}

void PreDefinedBodyForce3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const BodyForceType forceType{compo.bodyForceType};
            switch (forceType) {
                case BodyForce_1st:
                    ops_par_loop(
                        KerCalcBodyForce1ST3D, "KerCalcBodyForce1ST3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case BodyForce_None:
                    ops_par_loop(
                        KerCalcBodyForceNone3D, "KerCalcBodyForceNone",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified force type is not implemented!\n");
                    break;
            }
        }
    }
#endif // OPS_3D
}

void PreDefinedInitialCondition3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int compoId{compo.id};
            const InitialType initialType{compo.initialType};
            switch (initialType) {
                case Initial_BGKFeq2nd: {
                    ops_par_loop(
                        KerInitialiseBGK2nd3D, "KerInitialiseBGK2nd3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                default:
                    ops_printf(
                        "The specified initial type is not implemented!\n");
                    break;
            }
        }
    }
    // TODO this may be better arranged.
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
#endif // OPS_3D
}

#endif // OPS_3D outter
#ifdef OPS_2D
void PreDefinedCollision() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKIsothermal, "KerCollideBGKIsothermal",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKThermal4th:
                    ops_par_loop(
                        KerCollideBGKThermal, "KerCollideBGKThermal",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_T).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified collision type is not implemented!\n");
                    break;
            }
        }
    }
#endif // OPS_2D
}



void UpdateMacroVars() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            for (auto& macroVar : compo.macroVars) {
                const int varId{macroVar.second.id};
                const VariableTypes varType{macroVar.first};
                switch (varType) {
                    case Variable_Rho:
                        ops_par_loop(
                            KerCalcDensity, "KerCalcDensity", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_U:
                        ops_par_loop(
                            KerCalcU, "KerCalcU", block.Get(), SpaceDim(),
                            iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_V:
                        ops_par_loop(
                            KerCalcV, "KerCalcV", block.Get(), SpaceDim(),
                            iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_U_Force:
                        ops_par_loop(
                            KerCalcUForce, "KerCalcUForce", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_V_Force:
                        ops_par_loop(
                            KerCalcVForce, "KerCalcVForce", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    default:
                        break;
                }
            }
        }
    }
#endif // OPS_2D
}


void PreDefinedBodyForce() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const BodyForceType forceType{compo.bodyForceType};
            switch (forceType) {
                case BodyForce_1st:
                    ops_par_loop(
                        KerCalcBodyForce1ST, "KerCalcBodyForce1ST",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case BodyForce_None:
                    ops_par_loop(
                        KerCalcBodyForceNone, "KerCalcBodyForceNone",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified force type is not implemented!\n");
                    break;
            }
        }
    }
#endif // OPS_2D
}



void PreDefinedInitialCondition() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int compoId{compo.id};
            const InitialType initialType{compo.initialType};
            switch (initialType) {
                case Initial_BGKFeq2nd: {
                    ops_par_loop(
                        KerInitialiseBGK2nd, "KerInitialiseBGK2nd",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                default:
                    ops_printf(
                        "The specified initial type is not implemented!\n");
                    break;
            }
        }
    }
    // TODO this may be better arranged.
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
#endif // OPS_2D
}
#endif //OPS_2D outter