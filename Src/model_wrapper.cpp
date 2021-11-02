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
            const Real* pdt{pTimeStep()};
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
                case GuoForce:
                    ops_par_loop(
                        KerCalcGuoForce3D, "KerCalcGuoForce",
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
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_U_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_V_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_W_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
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

void PreDefinedInitialConditionAD3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const auto& it = g_Components().begin();
            const auto& Iter{*it};
            const Component& compoVel{Iter.second};
            const int compoId{compo.id};
            const int compoVelId{compoVel.id};
            const InitialType initialType{compo.initialType};

            const Real* pdt{pTimeStep()};
            switch (initialType) {
                case Initial_BGKFeq2ndAD: {
                    ops_par_loop(
                        KerInitialiseBGKADF3D, "KerInitialiseBGKAD3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),    
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                case Initial_BGKFeq2ndFE: {
                    ops_par_loop(
                        KerInitialiseBGKFEF3D, "KerInitialiseBGKAD3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),    
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                case Initial_BGKGeq2ndFE: {
                    ops_par_loop(
                        KerInitialiseBGKFEG3D, "KerInitialiseBGKADG3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                                    OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
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

void PreDefinedCollisionAD3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const auto& compoVel = g_Components().at(0);
            //const auto& Iter{*it};
            //const Component& compoVel{it.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKADF3D, "KerCollideBGKIsothermal",
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
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKADF:
                    ops_par_loop(
                        KerCollideBGKADF3D, "KerCollideBGKADF",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKFEF:
                    ops_par_loop(
                        KerCollideBGKFEF3D, "KerCollideBGKADF",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKFEG:
                    ops_par_loop(
                        KerCollideBGKFEG3D, "KerCollideBGKADG",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.wId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                                    OPS_READ),
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

//Subroutine to calculate central finite difference gradients at each lattice point
void Calcphi2Gradients3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{2};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        const int compoRhoId{compoRho.id};
        */
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcGradients3D, "KerCalcGradients", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIndex), 1,
                                    ONEPTREGULARSTENCIL, "Real", OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_gbl(compoRho.index, 2, "int", OPS_READ),
            ops_arg_gbl(&order, 1, "int", OPS_READ),
            ops_arg_idx());
        
    }
}

//Subroutine to update the chemical potential at each lattice point
void CalcMu3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};

        const Real A{AFE()};
        const Real Kappa{KappaFE()};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        const int compoRhoId{compoRho.id};
        */
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcMu3D, "KerCalcMu", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIndex), 1,
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(),
                        LOCALSTENCIL, "double", OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_gbl(&A, 1, "Real", OPS_READ),
            ops_arg_gbl(&Kappa, 1, "Real", OPS_READ));
            

        
    }
} 

//Subroutine to calculate the chemical potential gradients at each lattice point
void CalcmuGradients3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{1};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        */
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcGradients3D, "KerCalcGradients", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_muGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_mu()[blockIndex], 1, ONEPTREGULARSTENCIL, "double",
                        OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_gbl(compoRho.index, 2, "int", OPS_READ),
            ops_arg_gbl(&order, 1, "int", OPS_READ),
            ops_arg_idx());

        
    }
} 



//Subroutine to implement wetting boundary condition along all solid boundaries
void CalcPhiWetting3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real theta{Theta()};
        const Real A{AFE()};
        const Real Kappa{KappaFE()};
        auto compo = g_Components().at(1);

        ops_par_loop(KerUpdateRhoWetting3D, "KerUpdateRhoWetting3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIndex), 1,
                                    ONEPTLATTICESTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                    ONEPTREGULARSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_gbl(&theta, 1, "Real", OPS_READ),
                        ops_arg_gbl(&A, 1, "Real", OPS_READ),
                        ops_arg_gbl(&Kappa, 1, "Real", OPS_READ));
        
    }
}



//Subroutine to update macroscopic body force at each lattice point
void UpdateMacroscopicBodyForceFE3D(const Real time) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        /*
        auto it = g_Components().begin();
        const auto& Iter1{*it};
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        */
        auto compoVel = g_Components().at(0);
        auto compoRho = g_Components().at(1);
        ops_par_loop(KerUpdateMacroBodyForce3D, "KerUpdateMacroBodyForce",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroBodyforce().at(compoVel.id).at(blockIdx), SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIdx), 1,
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_dat(g_muGrad()[blockIdx], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_dat(g_GeometryProperty().at(blockIdx), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIdx), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_idx());

        
    }
}

//Toggle-able Free Energy implementation
void FreeEnergy3D(){

    g_mu().TransferHalos();
    const auto& compoC = g_Components().at(1);
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).TransferHalos();
    CalcPhiWetting3D();
    Calcphi2Gradients3D();
    CalcMu3D();
    g_mu().TransferHalos();
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).TransferHalos();
    CalcmuGradients3D();

    UpdateMacroscopicBodyForceFE3D(1);

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
#endif  // OPS_2D
}

void PreDefinedCollisionAD() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const auto& compoVel = g_Components().at(0);
            //const auto& Iter{*it};
            //const Component& compoVel{it.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            auto compoRho = g_Components().at(1);
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKADF, "KerCollideBGKIsothermal",
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
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKADF:
                    ops_par_loop(
                        KerCollideBGKADF, "KerCollideBGKADF",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKFEF:
                    ops_par_loop(
                        KerCollideBGKFEF, "KerCollideBGKADF",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_BGKFEG:
                    ops_par_loop(
                        KerCollideBGKFEG, "KerCollideBGKADG",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, ONEPTREGULARSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.macroVars.at(Variable_U_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.macroVars.at(Variable_V_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                                    OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;
                case Collision_MRTFEF:
                    ops_par_loop(
                        KerCollideMRTFEF, "KerCollideMRTFEF",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_RW),
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
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                    break;

                    /*
                    case Collision_BGKIsothermal2nd_diffusive:
                        ops_par_loop(
                            KerCollideBGKIsothermal, "KerCollideBGKIsothermal",
                            block.Get(), SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                        LOCALSTENCIL, "double", OPS_WRITE),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars().at(compo.uId).at(blockIndex), 1,
                                LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars().at(compo.vId).at(blockIndex), 1,
                                LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(&tau, 1, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                        */

                    default:
                        ops_printf(
                            "The specified collision type is not "
                            "implemented!\n");
                        break;
            
            }
        }
    }
#endif // OPS_2D
}

//g_MacroVars().at(varId).at(blockIndex)
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
            const Real* pdt{pTimeStep()};
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
                case GuoForce:
                    ops_par_loop(
                        KerCalcGuoForce, "KerCalcGuoForce",
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
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_U_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_V_Force).id).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
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

void PreDefinedInitialConditionAD() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const auto& it = g_Components().begin();
            const auto& Iter{*it};
            const Component& compoVel{Iter.second};
            const int compoId{compo.id};
            const int compoVelId{compoVel.id};
            const InitialType initialType{compo.initialType};

            const Real* pdt{pTimeStep()};
            switch (initialType) {
                case Initial_BGKFeq2ndAD: {
                    ops_par_loop(
                        KerInitialiseBGKADF, "KerInitialiseBGKAD",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),    
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                case Initial_BGKFeq2ndFE: {
                    ops_par_loop(
                        KerInitialiseBGKFEF, "KerInitialiseBGKAD",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),    
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                } break;
                case Initial_BGKGeq2ndFE: {
                    ops_par_loop(
                        KerInitialiseBGKFEG, "KerInitialiseBGKADG",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                    1, ONEPTREGULARSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                                    OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
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

//Subroutine to calculate central finite difference gradients at each lattice point
void Calcphi2Gradients() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{2};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        const int compoRhoId{compoRho.id};
        */
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcGradients, "KerCalcGradients", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIndex), 1,
                                    ONEPTREGULARSTENCIL, "Real", OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_gbl(compoRho.index, 2, "int", OPS_READ),
            ops_arg_gbl(&order, 1, "int", OPS_READ),
            ops_arg_idx());
        
    }
}

//Subroutine to update the chemical potential at each lattice point
void CalcMu() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        const int compoRhoId{compoRho.id};
        */

        const Real A{AFE()};
        const Real Kappa{KappaFE()};
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcMu, "KerCalcMu", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIndex), 1,
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(),
                        LOCALSTENCIL, "double", OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_gbl(&A, 1, "Real", OPS_READ),
            ops_arg_gbl(&Kappa, 1, "Real", OPS_READ));
            

        
    }
} 

//Subroutine to calculate the chemical potential gradients at each lattice point
void CalcmuGradients() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{1};
        /*
        auto it = g_Components().begin();
        std::advance(it, 1);
        const auto& Iter2{*it};
        */
        auto compoRho = g_Components().at(1);
        ops_par_loop(
            KerCalcGradients, "KerCalcGradients", block.Get(), SpaceDim(),
            iterRng.data(),
            ops_arg_dat(g_muGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                        OPS_RW),
            ops_arg_dat(g_mu()[blockIndex], 1, ONEPTREGULARSTENCIL, "double",
                        OPS_READ),
            ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                        LOCALSTENCIL, "int", OPS_READ),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
            ops_arg_gbl(compoRho.index, 2, "int", OPS_READ),
            ops_arg_gbl(&order, 1, "int", OPS_READ),
            ops_arg_idx());

        
    }
} 



//Subroutine to implement wetting boundary condition along all solid boundaries
void CalcPhiWetting() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        auto compo = g_Components().at(1);
        const Real theta{Theta()};
        const Real A{AFE()};
        const Real Kappa{KappaFE()};

        ops_par_loop(KerUpdateRhoWetting, "KerUpdateRhoWetting",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIndex), 1,
                                    ONEPTLATTICESTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                    ONEPTREGULARSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_gbl(&theta, 1, "Real", OPS_READ),
                        ops_arg_gbl(&A, 1, "Real", OPS_READ),
                        ops_arg_gbl(&Kappa, 1, "Real", OPS_READ));
        
    }
}


//Subroutine to update macroscopic body force at each lattice point
void UpdateMacroscopicBodyForceFE(const Real time) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        auto compoVel = g_Components().at(0);
        auto compoRho = g_Components().at(1);
        /*
        auto it = g_Components().begin();
        const auto& Iter1{*it};
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        */

        ops_par_loop(KerUpdateMacroBodyForce, "KerUpdateMacroBodyForce",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroBodyforce().at(compoVel.id).at(blockIdx), SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars()
                                .at(compoRho.macroVars.at(Variable_Rho).id)
                                .at(blockIdx), 1,
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_dat(g_muGrad()[blockIdx], SpaceDim(),
                        LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_GeometryProperty().at(blockIdx), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIdx), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_idx());

        
    }
}

//Toggle-able Free Energy implementation
void FreeEnergy2D(){
    g_mu().TransferHalos();
    const auto& compoC = g_Components().at(1);
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).TransferHalos();
    CalcPhiWetting();
    Calcphi2Gradients();
    CalcMu();
    g_mu().TransferHalos();
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).TransferHalos();
    CalcmuGradients();

    UpdateMacroscopicBodyForceFE(1);

}

#endif //OPS_2D outter